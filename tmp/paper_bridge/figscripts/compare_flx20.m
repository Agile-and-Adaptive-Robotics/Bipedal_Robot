%compare_flx20.m — paper FlxGrp panel-B corrected curve (minimizeFlx at
%row-77, restored mats) vs the refreshed Knee_Flexor_Data_20mm class
%pipeline (original two-point route, XiUsed = row-77 + ext Xi3, zero bend).
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(fullfile(root, 'Testing_Data/2022_02_Festo'));
warning('off', 'all'); setenv('FLX_T5_YMM', '');
S = load('KneeFlx_20mm_42cm_X3redesign_20261019.mat');   % refreshed pipeline
P = load('Plot_KneeFlx_20mm_42cm.mat', 'Angle', 'Torque');
g = [3.937370e-3, 3.997873e4, 1.473375e4];
[~, ~, ~, bpa] = minimizeFlx(g(1), g(2), g(3));          % paper pipeline
cls = S.Bifemsh_Pam0.Torque_p(:,3);
mf = bpa(2).M_p(:,3);
ph = S.phiD(:);
ang = P.Angle(:); Me = P.Torque(:);
Gc = griddedInterpolant(ph, cls); a = Gc(ang);
Gm = griddedInterpolant(ph, mf);  b = Gm(ang);
v = isfinite(a) & isfinite(b) & isfinite(Me);
fprintf('CMP finite triples %d/%d\n', sum(v), numel(v));
fprintf('CMP max|class - minimizeFlx| at measured angles = %.4f N*m\n', max(abs(a(v)-b(v))));
fprintf('CMP class-corrected     RMSE vs measured = %.4f\n', sqrt(mean((a(v)-Me(v)).^2)));
fprintf('CMP minimizeFlx-correct RMSE vs measured = %.4f\n', sqrt(mean((b(v)-Me(v)).^2)));
fprintf('CMP class mat Xi: %.6g %.6g %.6g %.6g ; rest0=%.4f tendon0=%.4f\n', ...
    S.Xi0, S.Xi1, S.Xi2, S.Xi3, S.rest0, S.tendon0);
disp('DONE_cmp');
