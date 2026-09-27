function retry_update_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
load_system(fullfile(sns, [mdl '.slx']));
try
    set_param(mdl, 'StopTime', '0.05');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
    out = sim(mdl);
    fprintf('SIM OK\n');
    save_system(mdl);
catch ME
    fprintf('UPDATE FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 4)
        fprintf('  CAUSE: %s\n', ME.cause{c}.message(1:min(end, 150)));
    end
end
close_system(mdl, 0);
end
