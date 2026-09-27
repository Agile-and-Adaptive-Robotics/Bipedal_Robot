function measure_frames_20260926()
% measure world orientations of the two hinge-source ports (quaternion)
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

% find body ports again (same ports v4 used)
bodyOf = @(h) get_param(get_param(h, 'Parent'), 'Name');
% branch sources: use pinT_KT / pinT_KB B-port lines' body-side handles
srcs = { 'pinT_KT', 'x04_01_KT_R_003_1_RIGID'; 'pinT_KB', 'x04_02_KB_R_003_1_RIGID' };
for s = 1:2
    ph = get_param([mdl '/' srcs{s, 1}], 'PortHandles');
    p = [ph.RConn ph.LConn];
    bodyPort = p(2);   % the pin frame's F output = the frame whose orientation we set
    % TS: base = World, follower = body port, SenseQ
    ts = ['qTS_' srcs{s, 1}];
    add_block('sm_lib/Frames and Transforms/Transform Sensor', [mdl '/' ts], 'Position', [1100 100+200*s 1180 180+200*s]);
    set_param([mdl '/' ts], 'SenseQ', 'on');
    w = find_system(mdl, 'SearchDepth', 1, 'Regexp', 'on', 'Name', 'World');
    wph = get_param(w{1}, 'PortHandles');
    tsph = get_param([mdl '/' ts], 'PortHandles');
    tsP = [tsph.LConn tsph.RConn];
    add_line(mdl, wph.RConn(1), tsP(1));
    add_line(mdl, bodyPort, tsP(2));
    cv = ['qC_' srcs{s, 1}];
    add_block('nesl_utility/PS-Simulink Converter', [mdl '/' cv], 'Position', [1220 100+200*s 1270 130+200*s]);
    cph = get_param([mdl '/' cv], 'PortHandles');
    add_line(mdl, tsP(3), cph.LConn(1));   % first enabled output = Q
    lg = ['qL_' srcs{s, 1}];
    add_block('simulink/Sinks/To Workspace', [mdl '/' lg], 'Position', [1310 100+200*s 1370 130+200*s]);
    set_param([mdl '/' lg], 'VariableName', ['q_' srcs{s, 1}], 'SaveFormat', 'Structure With Time');
    lph = get_param([mdl '/' lg], 'PortHandles');
    add_line(mdl, cph.Outport(1), lph.Inport(1));
end
set_param(mdl, 'StopTime', '0.01');
out = sim(mdl);
q1 = out.q_pinT_KT.signals(1).values(1, :);
q2 = out.q_pinT_KB.signals(1).values(1, :);
fprintf('KT port frame quaternion (w x y z): %.6f %.6f %.6f %.6f\n', q1);
fprintf('KB port frame quaternion (w x y z): %.6f %.6f %.6f %.6f\n', q2);
% cleanup
for s = 1:2
    for nm = {['qTS_' srcs{s, 1}], ['qC_' srcs{s, 1}], ['qL_' srcs{s, 1}]}
        try, delete_block([mdl '/' nm{1}]); catch, end
    end
end
save_system(mdl);
close_system(mdl, 0);
save(fullfile(here, 'frame_quats_20260926.mat'), 'q1', 'q2');
fprintf('=== measure_frames DONE ===\n');
end
