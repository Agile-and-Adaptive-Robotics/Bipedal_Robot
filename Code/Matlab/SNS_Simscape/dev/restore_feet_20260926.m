function restore_feet_20260926()
% remove the ankle retrofit, restore the original welded-feet chain
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

legs = {'R', 'x04_02_KB_R_001_1_RIGID', '[0.03992 -0.3323 0.00391]'; ...
        'L', 'x04_04_KB_L_001_1_RIGID', '[0.03992 0.3323 0.00391]'};
for k = 1:2
    S = legs{k, 1};
    for nm = {['ankle_' S], ['ankle_off_' S]}
        if getSimulinkBlockHandle([sub '/' nm{1}]) > 0
            delete_block([sub '/' nm{1}]);
        end
    end
    footT = [sub '/' ['foot_T_' S]];
    ph = get_param(footT, 'PortHandles');
    p = [ph.RConn ph.LConn];
    for pp = {p(1), p(2)}
        ln = get_param(pp{1}, 'Line');
        if ln > 0, delete_line(ln); end
    end
    set_param(footT, 'TranslationCartesianOffset', legs{k, 3});
    kbp = get_param([sub '/' legs{k, 2}], 'PortHandles');
    kbPort = [kbp.RConn kbp.LConn];
    kbPort = kbPort(1);
    fpT = get_param(getSimulinkBlockHandle(footT), 'PortHandles'); fpP = [fpT.RConn fpT.LConn];
    fPh = get_param([sub '/' ['foot_' S]], 'PortHandles');         fP2 = [fPh.RConn fPh.LConn];
    add_line(sub, kbPort, fpP(1));
    wire_verify(sub, fpP(2), fP2(1));
    fprintf('foot_%s restored to welded chain\n', S);
end
try
    set_param(mdl, 'StopTime', '0.05');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK - model back to good state\n');
    save_system(mdl);
    sim(mdl);
    fprintf('SIM OK\n');
catch ME
    fprintf('FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 4)
        fprintf('  CAUSE: %s\n', ME.cause{c}.message(1:min(end, 150)));
    end
    save_system(mdl);
end
close_system(mdl, 0);
fprintf('=== restore_feet DONE ===\n');
end

function wire_verify(sub, src, dst)
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
