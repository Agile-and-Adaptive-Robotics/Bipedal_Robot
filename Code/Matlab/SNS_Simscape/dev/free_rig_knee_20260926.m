function free_rig_knee_20260926()
% cut the KT<->KB Parallel angle constraint locking the rig four-bar knee,
% then run the SNS->BPA demo
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

try
    delete_block([mdl '/Parallel']);
    fprintf('deleted Parallel angle constraint (knee four-bar freed)\n');
catch ME
    fprintf('Parallel: %s\n', ME.message);
end

try
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
catch ME
    fprintf('UPDATE FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 6)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 170)));
    end
    save_system(mdl);
    close_system(mdl, 0);
    return;
end
save_system(mdl);
set_param(mdl, 'StopTime', '6');
out = sim(mdl);
LE = out.L_EXT.signals(1).values;
LF = out.L_FLX.signals(1).values;
FE = out.F_EXT.signals(1).values;
FF = out.F_FLX.signals(1).values;
VE = out.V_E.signals(1).values;
fprintf('EXT length: %.4f -> %.4f..%.4f m\n', LE(1), min(LE), max(LE));
fprintf('FLX length: %.4f -> %.4f..%.4f m\n', LF(1), min(LF), max(LF));
fprintf('EXT force : max %.1f N\n', max(abs(FE)));
fprintf('FLX force : max %.1f N\n', max(abs(FF)));
fprintf('V_E: %.1f..%.1f mV | knee swing dL = %.4f m\n', min(VE), max(VE), max(LE) - min(LE));
r = corrcoef(LE(:), LF(:));
fprintf('corr(L_EXT, L_FLX) = %.3f\n', r(1, 2));
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== free_rig_knee DONE ===\n');
end
