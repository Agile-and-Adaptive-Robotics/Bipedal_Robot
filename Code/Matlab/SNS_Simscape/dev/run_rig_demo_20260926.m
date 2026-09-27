function run_rig_demo_20260926()
% full SNS->BPA->knee demo on the rig: 6 s, report lengths/forces/pressures
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));
set_param(mdl, 'StopTime', '6');
out = sim(mdl);
LE = out.L_EXT.signals(1).values; tE = out.L_EXT.time;
LF = out.L_FLX.signals(1).values;
FE = out.F_EXT.signals(1).values;
FF = out.F_FLX.signals(1).values;
VE = out.V_E.signals(1).values; tV = out.V_E.time;
fprintf('EXT length: %.4f -> %.4f..%.4f m\n', LE(1), min(LE), max(LE));
fprintf('FLX length: %.4f -> %.4f..%.4f m\n', LF(1), min(LF), max(LF));
fprintf('EXT force : max %.1f N (final %.1f)\n', max(abs(FE)), FE(end));
fprintf('FLX force : max %.1f N (final %.1f)\n', max(abs(FF)), FF(end));
fprintf('V_E: %.1f..%.1f mV | knee swing = dL = %.4f m\n', min(VE), max(VE), max(LE) - min(LE));
% antiphase check between EXT and FLX lengths
r = corrcoef(LE(:), LF(:)); r = r(1, 2);
fprintf('corr(L_EXT, L_FLX) = %.3f (negative = antagonistic motion)\n', r);
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== run_rig_demo DONE ===\n');
end
