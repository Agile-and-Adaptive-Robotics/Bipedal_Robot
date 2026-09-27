function fall_test_20260926()
% sanity: 0.5 s sim of both models (humanoid free-falls, rig sits)
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
load_system(fullfile(sns, [mdl '.slx']));
try
    set_param(mdl, 'StopTime', '0.5');
    set_param(mdl, 'SimulationCommand', 'update');
    out = sim(mdl);
    fprintf('HUMANOID 0.5 s SIM OK (free fall expected, no sensors)\n');
catch ME
    fprintf('HUMANOID SIM FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 5), fprintf('  CAUSE: %s\n', ME.cause{c}.message(1:min(end,150))); end
end
save_system(mdl); close_system(mdl, 0);

mdl2 = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl2 '.slx']));
try
    set_param(mdl2, 'StopTime', '0.5');
    out = sim(mdl2);
    fprintf('RIG 0.5 s SIM OK\n');
catch ME
    fprintf('RIG SIM FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 5), fprintf('  CAUSE: %s\n', ME.cause{c}.message(1:min(end,150))); end
end
save_system(mdl2); close_system(mdl2, 0);
fprintf('=== fall_test DONE ===\n');
end
