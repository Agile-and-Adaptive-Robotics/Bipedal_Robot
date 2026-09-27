function fix_hinge2_20260926()
% Final rig-knee fix:
%   pinT_KT / pinT_KB / fl2T sit on x-flipped source frames (smiData ang=pi
%   about [1,0,0]). Give each an added pi-about-x rotation (net identity) and
%   re-express their translations in source axes.
%   Then place the muscles analytically on OPPOSITE sides of the pin
%   (EXT +x, FLX -x; origins 7 cm above, insertions 5.8 cm below the pin) and
%   skip calibration. BPA Rest/Kmax updated to the new path length.
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

% counter-rotation so F axes = part axes (source frames are Rx(pi))
for nm = {'pinT_KT', 'pinT_KB', 'fl2T'}
    set_param([mdl '/' nm{1}], 'RotationMethod', 'ArbitraryAxis', ...
        'RotationAngleUnits', 'rad', 'RotationAngle', '3.14159265358979', ...
        'RotationArbitraryAxis', '[1 0 0]');
end
% translations re-expressed in SOURCE axes (Rx(pi) applied to part-frame deltas)
set_param([mdl '/pinT_KT'], 'TranslationCartesianOffset', '[-0.0061 0.0077 0.0745]');
set_param([mdl '/pinT_KB'], 'TranslationCartesianOffset', '[0.0281 -0.0242 -0.0005]');
set_param([mdl '/fl2T'],    'TranslationCartesianOffset', '[-0.0178 -0.0152 0.0740]');

% muscles on opposite sides of the pin (offsets in pin-frame axes = part axes)
set_param([mdl '/EXT_origT'], 'TranslationCartesianOffset', '[0.028 0.070 0.0]', 'RotationMethod', 'None');
set_param([mdl '/EXT_insT'],  'TranslationCartesianOffset', '[0.033 -0.058 0.0]', 'RotationMethod', 'None');
set_param([mdl '/FLX_origT'], 'TranslationCartesianOffset', '[-0.028 0.070 0.0]', 'RotationMethod', 'None');
set_param([mdl '/FLX_insT'],  'TranslationCartesianOffset', '[-0.033 -0.058 0.0]', 'RotationMethod', 'None');

% BPA resting geometry: path = sqrt(0.005^2 + 0.128^2) = 0.1281 m
set_param([mdl '/EXT_BPA'], 'Rest', '0.1281', 'Kmax', '0.1057');
set_param([mdl '/FLX_BPA'], 'Rest', '0.1281', 'Kmax', '0.1057');
fprintf('pin frames counter-rotated; muscles placed on opposite sides\n');

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

% ---- demo ---------------------------------------------------------------------
set_param(mdl, 'StopTime', '6');
out = sim(mdl);
LE = out.L_EXT.signals(1).values;
LF = out.L_FLX.signals(1).values;
FE = out.F_EXT.signals(1).values;
FF = out.F_FLX.signals(1).values;
VE = out.V_E.signals(1).values;
fprintf('EXT length: %.4f -> %.4f..%.4f m | FLX: %.4f..%.4f m\n', LE(1), min(LE), max(LE), min(LF), max(LF));
fprintf('forces: EXT max %.1f N | FLX max %.1f N\n', max(abs(FE)), max(abs(FF)));
r = corrcoef(LE(:), LF(:));
fprintf('corr(L_EXT, L_FLX) = %.3f (want NEGATIVE)\n', r(1, 2));
fprintf('V_E: %.1f..%.1f mV\n', min(VE), max(VE));
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== fix_hinge2 DONE ===\n');
end
