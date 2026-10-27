%gen_Ext20mm_cam_20261026.m -- P1-a rebuild of the extensor 20 mm cam-route
% physical-validation figure from the DECLARED nine-point dataset
% (ExtTest20mm_2 tab; Knee_Extensor_20mm.m lines 590-595 arrays), with the
% audited series identities of audit_Ext20_nine_20261026.m:
%   measured    = indigo filled circles, adjoint reduction with the
%                 readserialnumbers2 extensor rig constants prf=(-0.050,
%                 -0.255) m, tab LC entries negated (force at -LC in the
%                 tibia frame), NO sign flip; matches the live torque
%                 guard, Knee_Extensor_20mm.m's compare section, and the
%                 xlsx "Torque actual" row (max diff 1.24 N*m).
%   back-calc   = gold circles, corrected cam model evaluated at each
%                 measured angle and pressure (pressure-specific
%                 force-model prediction; NOT the hybrid, which needs
%                 measured muscle length/moment arm not recorded here).
%   corrected   = magenta line, cam design of record at 620 kPa.
%   human       = black dotted vasti target.
% Supersedes gen_Ext20mm_cam_20261023.m (six points, flexor-template
% sign-flip with prf=-0.065 -- that series tracked the prediction curve,
% rms 5.1 N*m, which is why the plotted blue sat where the yellow
% back-calculated series belongs).
% Outputs: Knee_Torque_Test/Figures/Ext20mm_cam_20261026.eps (paper),
% Ext20mm_physical_20261026.pdf/.png (dissertation, return package), QC
% png, and the series-level audit CSV.
root = 'D:/GitHub/Bipedal_Robot';
resDir = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data/2022_02_Festo'), '-end');
warning('off', 'all');
cdeg = pi/180;

%% cam mat of record
C = load(fullfile(resDir, 'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), ...
    'predBest', 'ctxCamBest', 'XiUsed');
ctx = C.ctxCamBest;
pb = C.predBest;
TorqueCam620 = pb.TorqueZ(:);
phiD = ctx.phiD(:);
phi = ctx.phi(:);
T_t1_ICR = ctx.T_t1_ICR;
xiCtx = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
assert(isequal(xiCtx, C.XiUsed), 'Xi mismatch between ctxCamBest and XiUsed');

%% declared nine-point dataset (ExtTest20mm_2; tab entries as recorded)
Load9 = [180.2 181.5 168 181 116 100 126 173.4 140];   % N
K9    = [-26.2 -30 -61 -43 -119 5 -1 -22 -86.5];       % knee angle, deg
LC9   = [28 26 23 27 10.5 21 27 28 22];                % LC angle, deg (tab, positive)
P9    = [437 420 425 435 256 480 488 435 350];         % gauge pressure, kPa

%% audited adjoint reduction (guard convention)
p_rf = [-0.050, -0.255, 0]';
d = hypot(p_rf(1), p_rf(2));
ang = atan2d(p_rf(2), p_rf(1));
T_t1_rf = RpToTrans(eye(3), p_rf);
Trk = pagemtimes(TransInv(T_t1_rf), T_t1_ICR);
fcn15 = fit(phi, squeeze(Trk(1,4,:)), 'cubicspline');
fcn16 = fit(phi, squeeze(Trk(2,4,:)), 'cubicspline');
K_ang = K9 * cdeg;
LC_ang = -LC9 * cdeg;    % extensor tab entries negated (compare convention)
Tz9 = zeros(1, numel(Load9));
for i = 1:numel(Load9)
    Trk_i = RpToTrans(eye(3), [fcn15(K_ang(i)), fcn16(K_ang(i)), 0]');
    Fr = [0; 0; 0; Load9(i)*cos(LC_ang(i)); Load9(i)*sin(LC_ang(i)); 0];
    Fk = Adjoint(Trk_i)' * Fr;
    Tz9(i) = Fk(3);
end
fprintf('Reaction point (%.3f, %.3f) m: d = %.4f m, ang = %.2f deg\n', ...
    p_rf(1), p_rf(2), d, ang);
fprintf('Nine data measurements (N*m, positive extension):\n');
for i = 1:9
    fprintf('  ang %7.1f  load %6.1f N  LC %5.1f  P %5.1f kPa -> %+8.3f\n', ...
        K9(i), Load9(i), LC9(i), P9(i), Tz9(i));
end

%% predictions (cam design of record, rebuilt pipeline)
Name = 'Vastus Medialis Proximal Ring 20mm BPA';
rest = pb.rest; tendon = pb.tendon;
kmax = rest * (1 - pb.KMAX);
LocationAdj = pb.Location;
bendMeasureAdj = pb.bendMeasure;
CrossPoint = ctx.CrossPoint; Dia = ctx.Dia; fitting = ctx.fitting;
wraps = ctx.wraps; BPAcount = pb.BPAcount; T_Pam = ctx.T_Pam;

V620 = MonoPamDataExplicit_balanceX3(Name, LocationAdj, CrossPoint, Dia, ...
    T_Pam, rest, kmax, tendon, fitting, 620, ...
    ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, wraps, phiD, BPAcount, bendMeasureAdj);
dm = max(abs(V620.Torque_p(:,3) - TorqueCam620));
fprintf('AUDIT rebuild: max |MonoPam@620 - predBest.TorqueZ| = %.6f N*m\n', dm);

predAP = zeros(size(P9));
uP = unique(P9);
for u = 1:numel(uP)
    idx = P9 == uP(u);
    tmpP = MonoPamDataExplicit_balanceX3(Name, LocationAdj, CrossPoint, Dia, ...
        T_Pam, rest, kmax, tendon, fitting, uP(u), ...
        ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, wraps, phiD, BPAcount, bendMeasureAdj);
    predAP(idx) = interp1(phiD, tmpP.Torque_p(:,3), K9(idx), 'pchip');
end
pred620 = interp1(phiD, TorqueCam620, K9, 'pchip');
human9 = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), K9, 'pchip');

fprintf('\n=== audit table (N*m; meas - pred@P; meas - human) ===\n');
for i = 1:9
    fprintf('%7.1f %7.1f %6.1f %6.0f  meas %8.3f  pred@P %8.3f  pred620 %8.3f  human %8.3f  %+7.3f  %+7.3f\n', ...
        K9(i), Load9(i), LC9(i), P9(i), Tz9(i), predAP(i), pred620(i), ...
        human9(i), Tz9(i)-predAP(i), Tz9(i)-human9(i));
end
fprintf('meet-or-exceed human: %d of %d; margins %.2f to %.2f N*m\n', ...
    sum(Tz9 >= human9), 9, min(Tz9-human9), max(Tz9-human9));
fprintf('measured exceeds pred@P at %d of %d (max +%.2f; min %+.2f)\n', ...
    sum(Tz9 > predAP), 9, max(Tz9-predAP), min(Tz9-predAP));

%% CSV export (series level)
csvFile = fullfile(qcdir, 'Ext20mm_series_audit_20261026.csv');
fid = fopen(csvFile, 'w');
fprintf(fid, ['series_id,idx,K_deg,Load_N,LC_tab_deg,LC_convention,Pressure_kPa,', ...
    'reaction_x_m,reaction_y_m,reaction_d_m,reaction_ang_deg,measured_torque_Nm,', ...
    'pred_at_measured_P_Nm,pred_620kPa_Nm,human_target_Nm,meas_minus_predP,', ...
    'meas_minus_human,classification\n']);
for i = 1:9
    fprintf(fid, 'ExtTest20mm_2,%d,%.1f,%.1f,%.1f,tibia-axis(tab positive),%.0f,%.3f,%.3f,%.4f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%s\n', ...
        i, K9(i), Load9(i), LC9(i), P9(i), p_rf(1), p_rf(2), d, ang, ...
        Tz9(i), predAP(i), pred620(i), human9(i), Tz9(i)-predAP(i), ...
        Tz9(i)-human9(i), 'measured(adjoint,guard convention)');
end
fclose(fid);
fprintf('WROTE %s\n', csvFile);

%% figure
cMeas = [0 0 1];            % indigo
cGold = [1 0.843 0];        % gold (back-calculated)
cCorr = [0.804 0.204 0.710];% magenta
cHum  = [0 0 0];            % black
figT = figure('Name','Ext20mm_cam_20261026','Color','w');
figT.Position = [90 90 640 480];
ax = gca; hold(ax, 'on');
plot(ax, phiD, TorqueCam620, '-', 'Color', cCorr, 'LineWidth', 2.5, ...
    'DisplayName', 'Cam design of record, 620 kPa');
plot(ax, ctx.humanAngleD, ctx.humanTorque, ':', 'Color', cHum, 'LineWidth', 3, ...
    'DisplayName', 'Human target (vasti)');
scatter(ax, K9, predAP, 60*P9/620, 'filled', ...
    'MarkerFaceColor', cGold, 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Back-calculated (model at measured P)');
scatter(ax, K9, Tz9, 60*P9/620, 'filled', ...
    'MarkerFaceColor', cMeas, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Measured');
set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
    'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XLim', [-125 15], 'XTick', [-120 -90 -60 -30 0], ...
    'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
xlabel(ax, '\bf \theta_{\bf k} , \circ', 'Interpreter', 'tex');
ylabel(ax, '\bf M_{\bf k} , N\cdotm', 'Interpreter', 'tex');
lg = legend(ax, 'Location', 'northeast'); lg.FontSize = 8;
exportgraphics(figT, fullfile(outdir, 'Ext20mm_cam_20261026.eps'), 'ContentType','vector');
exportgraphics(figT, fullfile(qcdir, 'Ext20mm_physical_20261026.pdf'), 'ContentType','vector');
exportgraphics(figT, fullfile(qcdir, 'Ext20mm_physical_20261026.png'), 'Resolution', 300);
exportgraphics(figT, fullfile(qcdir, 'QC_Ext20mm_cam_20261026.png'), 'Resolution',150);
fprintf('WROTE Ext20mm_cam_20261026.eps + diss PDF/PNG + QC png\n'); close(figT);
disp('DONE_Ext20mm_cam_20261026');
