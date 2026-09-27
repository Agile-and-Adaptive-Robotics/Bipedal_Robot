function add_ankles_20260926()
%% add_ankles_20260926  Convert the welded block feet to REVOLUTE ankles.
%% Chain per side: KB port1 --(branch)--> ankle_T (translate to ankle point)
%%   -> ankle Revolute (sagittal, frame z = part z = mediolateral)
%%   -> foot_T (translate down half foot height) -> foot brick.
%% The old direct foot_T branch line is cut. Ankle joint sense q enabled
%% (position output for later afferent wiring).
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

% foot_T translations currently hold brick-centre offsets from the KB branch
% frames; the ankle point = brick centre + (0, +0.0175, 0)
legs = {'R', 'x04_02_KB_R_001_1_RIGID', 'x04_02_KB_R_001_1_RIGID'; ...
        'L', 'x04_04_KB_L_001_1_RIGID', 'x04_04_KB_L_001_1_RIGID'};
for k = 1:2
    S = legs{k, 1};
    footT = [sub '/' ['foot_T_' S]];
    footB = [sub '/' ['foot_' S]];
    ph = get_param(footT, 'PortHandles');
    p = [ph.RConn ph.LConn];
    % rebuild from scratch: delete the joint+offset blocks (takes their lines),
    % cut BOTH of foot_T's lines
    for nm = {['ankle_' S], ['ankle_off_' S]}
        if getSimulinkBlockHandle([sub '/' nm{1}]) > 0
            delete_block([sub '/' nm{1}]);
        end
    end
    for pp = {p(1), p(2)}
        ln = get_param(pp{1}, 'Line');
        if ln > 0, delete_line(ln); end
    end
    % find the branch SOURCE port again: the KB body port1 (KB<->TI weld line)
    kbp = get_param([sub '/' legs{k, 2}], 'PortHandles');
    kbPort = [kbp.RConn kbp.LConn];
    kbPort = kbPort(1);
    % re-purpose foot_T as the ankle-point transform
    if strcmp(S, 'R')
        set_param(footT, 'TranslationCartesianOffset', '[0.03992 -0.3148 0.00391]');
    else
        set_param(footT, 'TranslationCartesianOffset', '[0.03992 0.3148 0.00391]');
    end
    % fresh blocks
    offT = [sub '/' ['ankle_off_' S]];
    add_block('sm_lib/Frames and Transforms/Rigid Transform', offT, 'Position', [660 260+80*k 740 320+80*k]);
    set_param(offT, 'TranslationMethod', 'Cartesian', ...
        'TranslationCartesianOffset', '[0 0.0175 0]', 'RotationMethod', 'None');
    ank = [sub '/' ['ankle_' S]];
    add_block('sm_lib/Joints/Revolute Joint', ank, 'Position', [600 260+80*k 660 320+80*k]);
    % wire (all ends fresh: foot_T ports cut above, new blocks unconnected)
    fpT = get_param(getSimulinkBlockHandle(footT), 'PortHandles');  fpP = [fpT.RConn fpT.LConn];
    aPh = get_param(getSimulinkBlockHandle(ank), 'PortHandles');    aP2 = [aPh.RConn aPh.LConn];
    oPh = get_param(getSimulinkBlockHandle(offT), 'PortHandles');   oP2 = [oPh.RConn oPh.LConn];
    fPh = get_param(getSimulinkBlockHandle(footB), 'PortHandles');  fP2 = [fPh.RConn fPh.LConn];
    add_line(sub, kbPort, fpP(1));
    wire_verify(sub, fpP(2), aP2(1));
    wire_verify(sub, aP2(2), oP2(1));
    wire_verify(sub, oP2(2), fP2(1));
    fprintf('ankle_%s: revolute installed\n', S);
end

try
    set_param(mdl, 'StopTime', '0.05');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
    save_system(mdl);
    out = sim(mdl);
    fprintf('0.05 s SIM OK (ankles move under gravity)\n');
catch ME
    fprintf('UPDATE/SIM FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 6)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 170)));
    end
    save_system(mdl);
end
close_system(mdl, 0);
fprintf('=== add_ankles DONE ===\n');
end

function wire_verify(sub, src, dst)
% add_line; tolerate Simulink's throw-after-connect quirk by verifying the line
try
    add_line(sub, src, dst);
catch ME
    if get_param(dst, 'Line') > 0
        fprintf('wire threw but connected OK\n');
    else
        rethrow(ME);
    end
end
end

function wire_once(sub, src, dst)
% add_line unless the destination already has a connection (idempotent)
if get_param(dst, 'Line') > 0
    return;
end
try
    add_line(sub, src, dst);
catch ME
    % Simulink sometimes connects successfully and still throws
    if get_param(dst, 'Line') > 0
        fprintf('wire threw but connected OK (%s)\n', ME.message(1:min(end, 60)));
    else
        rethrow(ME);
    end
end
end
