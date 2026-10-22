% Required envelope + predicted torque at the six new test angles
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
S = load('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat');
pb = S.predBest;
fn = fieldnames(pb); fprintf('predBest fields: %s\n', strjoin(fn, ','));
ctx = S.ctx; fc = fieldnames(ctx); fprintf('ctx fields: %s\n', strjoin(fc, ','));
phiD = [];
if isfield(pb, 'phiD'), phiD = pb.phiD; end
if isempty(phiD) && isfield(ctx, 'phiD'), phiD = ctx.phiD; end
if isempty(phiD) && isfield(pb, 'TorqueZ'), phiD = linspace(-120, 10, numel(pb.TorqueZ)); end
if ~isfield(pb, 'TororqueZ') && isfield(pb, 'TorqueZ'), tau0 = pb.TorqueZ; end
ht = [];
if isfield(ctx, 'humanTorque'), ht = ctx.humanTorque; end
if isfield(ctx, 'humanAngleD'), phiD = ctx.humanAngleD; end
env = 1.05 * abs(ht);
tau = abs(tau0);
angTest = [-26.2 -30 -61 -43 -119 5];
measT1 = [43.46 44.19 41.38 43.87 29.09 24.78];
measP  = [437 420 425 435 256 480];
fprintf('angle  envelope_Nm  predicted620_Nm  meas_t1_Nm  meas_kPa\n');
for i = 1:numel(angTest)
    e = interp1(phiD, env, angTest(i));
    t = interp1(phiD, tau, angTest(i));
    fprintf('%7.1f  %9.2f  %13.2f  %11.2f  %8g\n', angTest(i), e, t, measT1(i), measP(i));
end
fprintf('\nsingle-BPA-at-pressure = predicted620/2 * (P/620):\n');
for i = 1:numel(angTest)
    t = interp1(phiD, tau, angTest(i)) / 2 * measP(i) / 620;
    fprintf('%7.1f  %8.2f  (measured %8.2f, ratio %.2f)\n', angTest(i), t, measT1(i), measT1(i)/t);
end
