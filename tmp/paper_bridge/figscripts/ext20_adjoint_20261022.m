% ExtTest20mm_2 adjoint measured-torque: runner that reproduces the
% Knee_Extensor_20mm.m _2 section without the bone/GIF pipeline.
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization');
addpath('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo', '-end');
resDir = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results';
figOut = 'D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Figures';
qcDir  = 'D:/GitHub/Bipedal_Robot/tmp/paper_bridge/figscripts';

cdeg = pi/180;
S = load(fullfile(resDir, 'Vas_Pam_20mm_Result.mat'), 'xBest', 'XiUsed', 'ctx');
if isfield(S, 'XiUsed'), xiRun = S.XiUsed;
else, xiRun = [S.ctx.Xi0, S.ctx.Xi1, S.ctx.Xi2, S.ctx.Xi3]; end
ctx = buildKneeExtContext20mm();
xiBuilt = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
assert(isequal(xiBuilt, xiRun), 'Xi mismatch');
xBest = reshape(S.xBest, 1, []);
p1 = xBest(1:3); pEnd = xBest(4:6);
rest = xBest(7); tendon = xBest(8);
KMAX = ctx.KMAX; kmax = rest * (1 - KMAX);

% adj (cam) route of record -- constants verbatim from Knee_Extensor_20mm.m
tendonAdj = 0.092490;
p6Adj = [0.116323, 0.015040];
p7Adj = [0.131323, -0.018974];
p8Adj = [0.127039, -0.047873];
ctxAdj = ctx;
ctxAdj.geo.seedP6XY = p6Adj;
ctxAdj.geo.seedP7XY = p7Adj;
ctxAdj.geo.seedP8XY = p8Adj;
ctxAdj.geo.tibiaWrapCircumcircle = true;
[LocationAdj, bendMeasureAdj] = buildDistalRingLocation20mm(p1, pEnd, tendonAdj, ctxAdj);

Name = 'Vastus Medialis Proximal Ring 20mm BPA';
Vas_adj = MonoPamDataExplicit_balanceX3(Name, LocationAdj, ctx.CrossPoint, ...
    ctx.Dia, ctx.T_Pam, rest, kmax, tendonAdj, ctx.fitting, ctx.targetPressure, ...
    ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, ctx.wraps, ctx.phiD, ctx.BPAcount, bendMeasureAdj);
TorqueAdj620 = Vas_adj.Torque_p(:, 3);

phi = ctx.phi; phiD = ctx.phiD;
T_t1_ICR = ctx.T_t1_ICR;

% ExtTest20mm_2 data (sheet rows 6/7/8/10)
LoadE2 = [180.2 181.5 168 181 116 100];
KEd2   = [-26.2 -30 -61 -43 -119 5];
PE2    = [437 420 425 435 256 480];
LCEd2  = [28 26 23 27 10.5 21];
pasted = [50.58 32.8 28.8 NaN NaN NaN];   % row 17 as entered
KE2 = KEd2 * cdeg;
dExt2 = 250.9 / 1000; angExt2 = -90.83;

p_rfE2 = [dExt2 * cosd(angExt2), dExt2 * sind(angExt2), 0]';
T_t1_rfE2 = RpToTrans(eye(3), p_rfE2);
TrkE2 = pagemtimes(TransInv(T_t1_rfE2), T_t1_ICR);
fcnE2x = fit(phi, squeeze(TrkE2(1, 4, :)), 'cubicspline');
fcnE2y = fit(phi, squeeze(TrkE2(2, 4, :)), 'cubicspline');
TorqueZE2 = zeros(1, numel(LoadE2));
for i = 1:numel(LoadE2)
    LCE2 = deg2rad(90 - LCEd2(i));
    TrkE2i = RpToTrans(eye(3), [fcnE2x(KE2(i)), fcnE2y(KE2(i)), 0]');
    Fr = -[0; 0; 0; LoadE2(i) * cos(pi - LCE2); LoadE2(i) * sin(pi - LCE2); 0];
    Fk = Adjoint(TrkE2i)' * Fr;
    TorqueZE2(i) = Fk(3);
end
fprintf('ExtTest20mm_2 Torque actual (adjoint, d=%.1f mm, ang=%.2f deg):\n', ...
    1000 * dExt2, angExt2);
for i = 1:numel(LoadE2)
    fprintf('  run %d: ang %7.1f deg, %6.1f N, %5.1f kPa -> %8.4f N*m', ...
        i, KEd2(i), LoadE2(i), PE2(i), TorqueZE2(i));
    if ~isnan(pasted(i))
        fprintf('  (pasted %g, delta %+.3f)', pasted(i), TorqueZE2(i) - pasted(i));
    end
    fprintf('\n');
end

% paper figure: cam-route 620 kPa + human target + measured (size ~ p)
c7 = '#0000FF'; c8 = '#000000'; adjColor = c7;
figT = figure('Name', 'Ext20mm_measured', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 16 11]);
ax = gca; hold(ax, 'on');
plot(ax, phiD, TorqueAdj620, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, ...
    'DisplayName', 'Cam design of record, 620 kPa');
plot(ax, ctx.humanAngleD, ctx.humanTorque, ':', 'Color', c8, 'LineWidth', 3, ...
    'DisplayName', 'Human target');
scatter(ax, KEd2, TorqueZE2, 60 * PE2 / 620, 'filled', ...
    'MarkerFaceColor', c7, 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'none', ...
    'DisplayName', 'Measured (adjoint)');
set(ax, 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold', ...
    'LineWidth', 2, 'XLim', [min(phiD) max(phiD)], 'TickLength', [0.025 0.05]);
xlabel(ax, '\theta_k, \circ', 'Interpreter', 'tex', 'FontWeight', 'bold');
ylabel(ax, 'Torque, N\cdotm', 'Interpreter', 'tex', 'FontWeight', 'bold');
lg = legend(ax, 'Location', 'northeast');
set(lg, 'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', 'Box', 'off');
title(ax, '\bf 20 mm extensor route: measured vs predicted', 'Interpreter', 'tex');
exportgraphics(figT, fullfile(figOut, 'Ext20mm_measured.eps'), 'ContentType', 'vector');
exportgraphics(figT, fullfile(qcDir, 'QC_Ext20mm_measured.png'), 'Resolution', 150);
fprintf('WROTE Ext20mm_measured.eps + QC png\n');

% comparison numbers for prose
env = interp1(ctx.humanAngleD, 1.05 * abs(ctx.humanTorque), KEd2, 'pchip');
fprintf('\nangle  meas_adj  envelope  ratio  pred620scaled(P)  ratio\n');
for i = 1:numel(LoadE2)
    tp = interp1(phiD, TorqueAdj620, KEd2(i)) * PE2(i) / 620;
    fprintf('%7.1f  %8.2f  %8.2f  %5.2f  %14.2f  %5.2f\n', ...
        KEd2(i), TorqueZE2(i), env(i), TorqueZE2(i) / env(i), tp, TorqueZE2(i) / tp);
end
