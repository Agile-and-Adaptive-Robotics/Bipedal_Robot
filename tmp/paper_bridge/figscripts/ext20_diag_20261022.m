% Read-only diagnostic: measured vs back-calculated at Ben's 9 points,
% with the implied pressure the model needs to reproduce each hold.
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization');
addpath('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo', '-end');
resDir = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results';
cdeg = pi/180;

% verbatim from Ben's current compare section
Load   = [180.2 181.5 168 181 116 100 126 173.4 140];
K_deg  = [-26.2 -30 -61 -43 -119 5 -1 -22 -86.5];
LC_deg = -[28 26 23 27 10.5 21 27 28 22];
Presh  = [437 420 425 435 256 480 488 435 350];
p_rf   = [-0.05, -0.255, 0]';
K_ang  = K_deg * cdeg;
LC_ang = LC_deg * cdeg;

% context + cam route of record (tendonAdj = 92.49 mm, current)
S = load(fullfile(resDir, 'Vas_Pam_20mm_Result.mat'), 'xBest', 'XiUsed', 'ctx');
if isfield(S, 'XiUsed'), xiRun = S.XiUsed;
else, xiRun = [S.ctx.Xi0, S.ctx.Xi1, S.ctx.Xi2, S.ctx.Xi3]; end
ctx = buildKneeExtContext20mm();
assert(isequal([ctx.Xi0 ctx.Xi1 ctx.Xi2 ctx.Xi3], xiRun));
xBest = reshape(S.xBest, 1, []);
p1 = xBest(1:3); pEnd = xBest(4:6); rest = xBest(7);
KMAX = ctx.KMAX; kmax = rest * (1 - KMAX);
tendonAdj = 0.092490;
ctxAdj = ctx;
ctxAdj.geo.seedP6XY = [0.116323, 0.015040];
ctxAdj.geo.seedP7XY = [0.131323, -0.018974];
ctxAdj.geo.seedP8XY = [0.127039, -0.047873];
ctxAdj.geo.tibiaWrapCircumcircle = true;
[LocationAdj, bendMeasureAdj] = buildDistalRingLocation20mm(p1, pEnd, tendonAdj, ctxAdj);
Name = 'Vastus Medialis Proximal Ring 20mm BPA';

phi = ctx.phi; phiD = ctx.phiD;
T_t1_ICR = ctx.T_t1_ICR;

% measured adjoint (Ben's fixed convention)
Trk = pagemtimes(TransInv(RpToTrans(eye(3), p_rf)), T_t1_ICR);
f15 = fit(phi, squeeze(Trk(1, 4, :)), 'cubicspline');
f16 = fit(phi, squeeze(Trk(2, 4, :)), 'cubicspline');
TorqueZ = zeros(1, numel(Load));
for i = 1:numel(Load)
    Ti = RpToTrans(eye(3), [f15(K_ang(i)); f16(K_ang(i)); 0]);
    Fr = [0; 0; 0; Load(i) * cos(LC_ang(i)); Load(i) * sin(LC_ang(i)); 0];
    Fk = Adjoint(Ti)' * Fr;
    TorqueZ(i) = Fk(3);
end

% back-calculated at (angle, pressure)
predAP = zeros(size(Presh));
uP = unique(Presh);
for u = 1:numel(uP)
    idx = Presh == uP(u);
    tmpP = MonoPamDataExplicit_balanceX3(Name, LocationAdj, ctx.CrossPoint, ...
        ctx.Dia, ctx.T_Pam, rest, kmax, tendonAdj, ctx.fitting, uP(u), ...
        ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3, ctx.wraps, phiD, ctx.BPAcount, ...
        bendMeasureAdj);
    predAP(idx) = interp1(phiD, tmpP.Torque_p(:, 3), K_deg(idx), 'pchip');
end

% human target at the angles, for reference
tgt = interp1(ctx.humanAngleD(:), ctx.humanTorque(:), K_deg, 'pchip');

fprintf(['angle   P_meas  tau_meas  tau_pred  meas/pred  P_implied  ' ...
    'tau_human  meas/human\n']);
for i = 1:numel(Load)
    r = TorqueZ(i) / predAP(i);
    fprintf('%6.1f  %6.0f  %8.2f  %8.2f  %9.2f  %9.0f  %9.2f  %9.2f\n', ...
        K_deg(i), Presh(i), TorqueZ(i), predAP(i), r, Presh(i) * r, ...
        tgt(i), TorqueZ(i) / tgt(i));
end
