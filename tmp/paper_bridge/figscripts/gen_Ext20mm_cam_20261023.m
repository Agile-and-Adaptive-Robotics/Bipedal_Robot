%gen_Ext20mm_cam_20261023.m -- extensor 20 mm cam-route hardware validation
% figure (Ben's display-item restructure, 2026-10-23): the six measured
% isometric holds of the ExtTest20mm_2 tab, reduced to knee torque with the
% CORRECTED reaction point (readserialnumbers2: x=-0.065 m, y=-0.255 m,
% d=0.2632 m, ang=-104.3 deg), against the cam design of record's corrected
% prediction (Vas_Pam_20mm_Result_cam_20261014_0416.mat, predBest at 620 kPa)
% and the vasti human target (ctxCamBest, the RouteVariants basis).
% The compare section of Knee_Extensor_20mm.m is READ-ONLY; its machinery is
% reproduced here with the corrected reaction point. The xlsx row-17 pasted
% torques (50.58/32.8/28.8) disagree with this recomputation and are NOT used.
% Output: Knee_Torque_Test/Figures/Ext20mm_cam_20261023.eps + QC png.
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
xiRun = C.XiUsed;
xiCtx = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
assert(isequal(xiCtx, xiRun), 'Xi mismatch between ctxCamBest and XiUsed');
TorqueCam620 = pb.TorqueZ(:);
phiD = ctx.phiD(:);
phi = ctx.phi(:);
T_t1_ICR = ctx.T_t1_ICR;
fprintf('Cam design of record: rest=%.4f m tendon=%.5f m BPAcount=%d\n', ...
    pb.rest, pb.tendon, pb.BPAcount);
fprintf('predBest.TorqueZ (620 kPa corrected): [%.2f %.2f] N*m over %d poses\n', ...
    min(TorqueCam620), max(TorqueCam620), numel(TorqueCam620));

%% six measured holds (ExtTest20mm_2 tab rows 6/7/8/10, tests 1-6)
Load6 = [180.2 181.5 168 181 116 100];
K6    = [-26.2 -30 -61 -43 -119 5];
P6    = [437 420 425 435 256 480];
LC6   = [28 26 23 27 10.5 21];
xlsxTabTorque = [44.8 45.6 42.7 45.2 30.1 25.6];   % tab row 16 (old basis)
xlsxPasted    = [50.58 32.8 28.8 NaN NaN NaN];     % tab row 17 (flagged, unused)

%% adjoint reduction: FLEXOR-TEMPLATE machinery (the measured t1-ICR spline
% table of the physical leg, from Knee_Flexor_Data_20mm.m / gen_Flx20mm_
% 3panel_20261023.m), corrected reaction point (-0.065,-0.255), LC entries
% positive as in the tab, force at LC+pi; the flexor-template sign is then
% flipped so extension torque reads positive.
positions = 100;
knee_angle_tbl = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	22.32	21.99	21.2	22.33]')/1000;
fcn13 = fit(knee_angle_tbl, t1_ICR_x, 'cubicspline');
fcn14 = fit(knee_angle_tbl, t1_ICR_y, 'cubicspline');
kneeMin = -2.0943951;
kneeMax = 0.17453293;
phiFlx = linspace(kneeMin, kneeMax, positions);
[valF, posF] = min(abs(phiFlx));
phiFlx(posF) = 0;
T_t1_ICR_flx = zeros(4, 4, positions);
for i = 1:positions
    t1toICR = [fcn13(phiFlx(i)), fcn14(phiFlx(i)), 0];
    T_t1_ICR_flx(:, :, i) = RpToTrans(eye(3), t1toICR');
end
phiFit = phiFlx(:);

tzCorr = extFlxTemplateTorques([-0.065, -0.255, 0]', Load6, K6, LC6, phiFit, T_t1_ICR_flx);
tzOld  = extFlxTemplateTorques([-0.050, -0.255, 0]', Load6, K6, LC6, phiFit, T_t1_ICR_flx); %#ok<NASGU>

% task-stated targets, angle-ordered: -119 -61 -43 -26.2 -30 +5
[~, iOrd] = sort(K6);
taskT = [37.23 38.85 37.99 38.31 30.13 22.86];
fprintf('\n=== six adjoint torques, flexor template + corrected point, raw sign ===\n');
for k = 1:6
    fprintf('ang %7.1f  load %6.1f N  P %5.1f kPa  LC %5.1f -> tz %+8.4f N*m\n', ...
        K6(k), Load6(k), P6(k), LC6(k), tzCorr(k));
end
fprintf('sorted by angle: '); fprintf('%8.4f', tzCorr(iOrd)); fprintf('\n');
fprintf('task values    : '); fprintf('%8.2f ', taskT); fprintf('\n');
fprintf('max |flipped - task| = %.4f N*m\n', max(abs(-tzCorr(iOrd) - taskT)));

% sign convention: extension torque positive (flip the flexor-template sign)
Tz6 = -tzCorr(:);
fprintf('\n=== ADOPTED (sign-flipped to positive extension) vs tab row 16 ===\n');
for k = 1:6
    fprintf('ang %7.1f  meas %8.4f   tab16 %6.1f   delta %+.2f', ...
        K6(k), Tz6(k), xlsxTabTorque(k), Tz6(k) - xlsxTabTorque(k));
    if ~isnan(xlsxPasted(k))
        fprintf('   pasted %6.2f delta %+.2f', xlsxPasted(k), Tz6(k) - xlsxPasted(k));
    end
    fprintf('\n');
end

%% predictions at the measured angles and pressures (rebuilt via the
%% established pipeline on the saved design-of-record route arrays)
Name = 'Vastus Medialis Proximal Ring 20mm BPA';
rest = pb.rest; tendon = pb.tendon;
kmax = rest * (1 - pb.KMAX);
LocationAdj = pb.Location;
bendMeasureAdj = pb.bendMeasure;
CrossPoint = ctx.CrossPoint; Dia = ctx.Dia; fitting = ctx.fitting;
wraps = ctx.wraps; BPAcount = pb.BPAcount;
T_Pam = ctx.T_Pam;

V620 = MonoPamDataExplicit_balanceX3(Name, LocationAdj, CrossPoint, Dia, ...
    T_Pam, rest, kmax, tendon, fitting, 620, ...
    ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, wraps, phiD, BPAcount, bendMeasureAdj);
dm = max(abs(V620.Torque_p(:,3) - TorqueCam620));
fprintf('\nAUDIT rebuild: max |MonoPam@620 - predBest.TorqueZ| = %.6f N*m\n', dm);

predAP = zeros(size(P6));
uP = unique(P6);
for u = 1:numel(uP)
    idx = P6 == uP(u);
    tmpP = MonoPamDataExplicit_balanceX3(Name, LocationAdj, CrossPoint, Dia, ...
        T_Pam, rest, kmax, tendon, fitting, uP(u), ...
        ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, wraps, phiD, BPAcount, bendMeasureAdj);
    predAP(idx) = interp1(phiD, tmpP.Torque_p(:,3), K6(idx), 'pchip');
end
pred620 = interp1(phiD, TorqueCam620, K6, 'pchip');
human6 = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), K6, 'pchip');

fprintf('\n=== comparison (all N*m, positive extension) ===\n');
fprintf('angle   load_N   P_kPa    meas   pred@P  pred@620  human  meas-pred@P  meas-human\n');
for k = 1:6
    fprintf('%7.1f %7.1f %7.1f %8.2f %8.2f %8.2f %8.2f %11.2f %11.2f\n', ...
        K6(k), Load6(k), P6(k), Tz6(k), predAP(k), pred620(k), human6(k), ...
        Tz6(k) - predAP(k), Tz6(k) - human6(k));
end
fprintf('meet-or-exceed vs human target at the held angles: %d of %d\n', ...
    sum(Tz6 >= human6), numel(Tz6));
fprintf('measured range: [%.2f %.2f] N*m; predicted@620 at held angles: [%.2f %.2f]\n', ...
    min(Tz6), max(Tz6), min(pred620), max(pred620));

%% figure (A): torque vs angle, measured + predicted + human dotted
c7 = [0 0 1]; c8 = [0 0 0]; c5 = [0.804 0.204 0.710];
figT = figure('Name','Ext20mm_cam','Color','w');
figT.Position = [90 90 640 480];
ax = gca; hold(ax, 'on');
plot(ax, phiD, TorqueCam620, '-', 'Color', c5, 'LineWidth', 2.5, ...
    'DisplayName', 'Cam design of record, 620 kPa');
plot(ax, ctx.humanAngleD, ctx.humanTorque, ':', 'Color', c8, 'LineWidth', 3, ...
    'DisplayName', 'Human target (vasti)');
scatter(ax, K6, Tz6, 60*P6/620, 'filled', ...
    'MarkerFaceColor', c7, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Measured');
set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
    'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XLim', [-125 15], 'XTick', [-120 -90 -60 -30 0], ...
    'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
xlabel(ax, '\bf \theta_{k} , \circ', 'Interpreter', 'tex');
ylabel(ax, '\bf Torque, N\cdotm', 'Interpreter', 'tex');
lg = legend(ax, 'Location', 'northeast'); lg.FontSize = 8;
exportgraphics(figT, fullfile(outdir, 'Ext20mm_cam_20261023.eps'), 'ContentType','vector');
exportgraphics(figT, fullfile(qcdir, 'QC_Ext20mm_cam_20261023.png'), 'Resolution',150);
fprintf('WROTE Ext20mm_cam_20261023.eps + QC png\n'); close(figT);
disp('DONE_Ext20mm_cam');

%% local: flexor-template adjoint (LC positive as entered, force at LC+pi)
function tz = extFlxTemplateTorques(prf, Load, Kdeg, LCdeg, phi, T_t1_ICR)
K_ang = Kdeg * pi/180;
LC_ang = LCdeg * pi/180;
p_rf = prf(:);
T_t1_rf = RpToTrans(eye(3), p_rf);
Trk = pagemtimes(TransInv(T_t1_rf), T_t1_ICR);
fcn15 = fit(phi, squeeze(Trk(1,4,:)), 'cubicspline');
fcn16 = fit(phi, squeeze(Trk(2,4,:)), 'cubicspline');
tz = zeros(1, numel(Load));
for i = 1:numel(Load)
    Trk_i = RpToTrans(eye(3), [fcn15(K_ang(i)), fcn16(K_ang(i)), 0]');
    Fr = [0; 0; 0; Load(i)*cos(LC_ang(i)+pi); Load(i)*sin(LC_ang(i)+pi); 0];
    Fk = Adjoint(Trk_i)' * Fr;
    tz(i) = Fk(3);
end
end
