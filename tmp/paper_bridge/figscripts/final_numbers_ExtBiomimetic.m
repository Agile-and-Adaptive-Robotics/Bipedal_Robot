%final_numbers_ExtBiomimetic.m — GoF of record for the 10mm/52cm biomimetic
%extensor at the restored routing, via the family-standard Go_OfF:
%original rigid, corrected at adopted row-32 terms, corrected at initial
%terms (-10.9mm, initial pair, 0.123), and hybrid (back-calculated).
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
festo = fullfile(root, 'Testing_Data/2022_02_Festo');
S = load(fullfile(festo, 'KneeExt_10mm_52cm.mat'));

ang = S.K_ang(:) * 180/pi;
Me = S.TorqueZ(:);
ph = S.phiD(:);

% corrected at INITIAL extensor terms (initial identification era)
init = MonoPamDataExplicit_balanceX3(S.Name, S.Location, S.CrossPoint, S.Dia, ...
    S.T_Pam, S.rest, S.kmax, S.tendon, S.fitting, S.pres, ...
    -0.0109, 1.02e5, 1.06e4, 0.123, S.wraps, ph, S.BPAcount);

ser = {S.TorqueR(:,3), 'Original (rigid)';
       S.Vas_Pam_adj.Torque_p(:,3), 'Corrected adopted (row 32)';
       init.Torque_p(:,3), 'Corrected initial (-10.9, 0.123)'};
for k = 1:size(ser,1)
    v = ser{k,1}(:);
    G = griddedInterpolant(ph, v);
    at = G(ang);
    [r, f, m] = Go_OfF(Me, at);
    fprintf('GOF %-32s RMSE %7.4f  FVU %7.4f  Max %7.4f  (finite-at-meas %d/9)\n', ...
        ser{k,2}, r, f, m, sum(isfinite(at)));
end

% hybrid with its own pressure trace (Plot-script recipe)
addpath(fullfile(festo, 'Ext_10mm'));
pres = zeros(1,9);
for j = 1:9
    T = load(sprintf('ExtTest1_%d.mat', j), 'Stats'); pres(j) = T.Stats{'Mean',2};
end
rmpath(fullfile(festo, 'Ext_10mm'));
InflatedLength = [518 520 510 495 485 475 460 450 450]/1000;
ICRtoMuscle = [38 35 40 40 45 45 40 40 41]/1000;
KMAX = (S.rest - S.kmax)/S.rest;
rel = ((S.rest - InflatedLength)/S.rest)/KMAX;
Th = ICRtoMuscle(:) .* bpaForce10(S.rest, rel, pres);
[r, f, m] = Go_OfF(Me, Th);
fprintf('GOF %-32s RMSE %7.4f  FVU %7.4f  Max %7.4f\n', 'Hybrid (back-calculated)', r, f, m);
disp('DONE_final');
