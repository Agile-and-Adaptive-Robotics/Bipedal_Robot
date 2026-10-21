%verify_ExtBiomimetic.m — clean re-verification of the ext-10mm-52cm GoF
%straight from the regenerated mat (no Plot-script workspace reuse).
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
S = load(fullfile(festodir, 'KneeExt_10mm_52cm.mat'));
fprintf('mat fields: %s\n', strjoin(fieldnames(S)', ', '));

ang = S.K_ang(:) * 180/pi;
Me = S.TorqueZ(:);
ph = S.phiD(:);
fprintf('measured angle: %s\n', mat2str(ang', 4));
fprintf('measured torque: %s\n', mat2str(Me', 4));

orig = S.TorqueR(:,3);
corr = S.Vas_Pam_adj.Torque_p(:,3);
fprintf('orig finite: %d/%d ; corr finite: %d/%d\n', sum(isfinite(orig)), numel(orig), sum(isfinite(corr)), numel(corr));
fprintf('orig at phi=0/90: %.4f / %.4f\n', orig(find(ph==0,1)), orig(min(sum(ph<=-90),numel(ph))));

G1 = griddedInterpolant(ph, orig); Mo = G1(ang);
G2 = griddedInterpolant(ph, corr); Mc2 = G2(ang);
fprintf('orig@meas: %s\n', mat2str(Mo', 4));
fprintf('corr@meas: %s\n', mat2str(Mc2', 4));

% hybrid: recompute per Plot_KneeExt10mm_52cm
InflatedLength = [518 520 510 495 485 475 460 450 450]/1000;
ICRtoMuscle = [38 35 40 40 45 45 40 40 41]/1000;
addpath(fullfile(festodir, 'Ext_10mm'));
pres = zeros(1,9);
for j = 1:9
    fn = sprintf('ExtTest%0.0f_%0.0f.mat', 1, j);
    T = load(fn, 'Stats'); pres(1,j) = T.Stats{'Mean',2};
end
rmpath(fullfile(festodir, 'Ext_10mm'));
KMAX = (S.rest - S.kmax)/S.rest;
rel = ((S.rest - InflatedLength)/S.rest)/KMAX;
Fn = bpaForce10(S.rest, rel, pres);
Th = ICRtoMuscle(:) .* Fn(:);
fprintf('pres: %s\n', mat2str(pres, 2));
fprintf('hybrid: %s\n', mat2str(Th', 4));

gof3 = @(p, m) [sqrt(mean((p(:)-m(:)).^2, 'omitnan')), ...
    nansum((p(:)-m(:)).^2)/nansum((m(:)-mean(m(:),'omitnan')).^2), ...
    max(abs(p(:)-m(:)), [], 'omitnan')];
fprintf('Original  [RMSE FVU Max]: %s\n', mat2str(gof3(Mo, Me), 5));
fprintf('Corrected [RMSE FVU Max]: %s\n', mat2str(gof3(Mc2, Me), 5));
fprintf('Hybrid    [RMSE FVU Max]: %s\n', mat2str(gof3(Th, Me), 5));
disp('DONE_verify');
