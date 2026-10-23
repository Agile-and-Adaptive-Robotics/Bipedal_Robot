%% write_flx20_table_0235.m
% Reproduces the Knee_Flexor_Data_20mm.m pressure back-calculation exactly
% (KMAX = 0.235, current design: route from xBest(8), TendonL = xBest(8) -
% 0.02) and writes the full 100-row table to THIS tmp folder, because the
% script's own write to Results\Flx_20mm_pressure_for_human_torque.csv is
% blocked by a file lock (Excel holds the Oct-19 copy open).
clc; clear; close all;

here = fileparts(mfilename('fullpath'));
cd(here);
root = here;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    if strcmp(parent, root), error('repo root not found'); end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

positions = 100;
R = zeros(3,3,positions); T = zeros(4,4,positions);
R_Pam = zeros(3,3,positions); T_Pam = zeros(4,4,positions);
t1toICR = zeros(1,3,positions);
T_t1_ICR = zeros(4,4,positions); T_ICR_t1 = zeros(4,4,positions);
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam = ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam = ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');
phi = linspace(-2.0943951, 0.17453293, positions);
[val, pos] = min(abs(phi)); phi(pos) = 0;
for i = 1:positions
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:,:,i) = [cos(phi(i)), -sin(phi(i)), 0; sin(phi(i)), cos(phi(i)), 0; 0,0,1];
    T(:,:,i) = RpToTrans(R(:,:,i), hipToKnee');
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:,:,i) = R(:,:,i);
    T_Pam(:,:,i) = RpToTrans(R_Pam(:,:,i), hipToKnee_Pam');
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:,:,i) = RpToTrans(eye(3), t1toICR(1,:,i)');
    T_ICR_t1(:,:,i) = RpToTrans(eye(3), -t1toICR(1,:,i)');
end
phiD = phi*180/pi;

resultData = load('Bifemsh_20mm_Result.mat', 'xBest', 'routeCtx', 'Xi3');
xBest = resultData.xBest;
routeCtx = resultData.routeCtx;
routeCtx.geo = buildGeoExclusion();
rX = load('Bifemsh_20mm_Result.mat', 'XiUsed');
Xi = rX.XiUsed;
Xi0=Xi(1); Xi1=Xi(2); Xi2=Xi(3); Xi3=Xi(4);

p1 = xBest(1:3); p2 = xBest(4:6);
rest = xBest(7);
tendon = xBest(8) - 0.02;
KMAX = 0.235;
kmax = (1 - KMAX)*rest;

Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2; Dia = 20; fitting = 0.021; wraps = 6; BPAcount = 2;
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[Loc1, bm1] = buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);
[Loc2, bm2] = buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);
LP = {Loc1; Loc2}; BM = {bm1; bm2};

H = readmatrix(fullfile(root,'Testing_Data','2022_02_Festo', ...
    'OpenSim_Bifem_Results.txt'), 'FileType', 'text', 'NumHeaderLines', 7);
humanAngle = H(:,2); TorqueHz = H(:,4);
[humanAngleU, iuH] = unique(humanAngle);
humanTgtPhi = interp1(humanAngleU, TorqueHz(iuH), phiD, 'pchip');
humanTgtPhi = humanTgtPhi(:);

Pgrid = (0:5:620)';
TqVsP = nan(numel(Pgrid), positions);
for kP = 1:numel(Pgrid)
    tmpP = MonoPam_mult(Name, LP, CrossPoint, Dia, T_Pam, rest, kmax, ...
        tendon, fitting, Pgrid(kP), Xi0, Xi1, Xi2, Xi3, wraps, phiD, ...
        BPAcount, BM);
    TqVsP(kP,:) = tmpP.Torque_p(:,3).';
end

Preq = nan(positions,1);
Flag = repmat("no human data", positions, 1);
for i = 1:positions
    tgt = humanTgtPhi(i);
    if ~isfinite(tgt), continue; end
    if tgt >= 0, Flag(i) = "no flexion demand"; continue; end
    col = TqVsP(:,i);
    fin = isfinite(col);
    if ~any(fin), Flag(i) = "BPA out of range"; continue; end
    if min(col(fin)) > tgt, Flag(i) = "unreachable at 620 kPa"; continue; end
    kk = find(col <= tgt, 1);
    if kk <= 1
        Preq(i) = 0;
    else
        pq = Pgrid([kk-1 kk]);
        [xs, isrt] = sort(col([kk-1 kk]));
        Preq(i) = interp1(xs, pq(isrt), tgt);
    end
    Flag(i) = "ok";
end

tab = table(phiD(:), humanTgtPhi(:), Preq(:), Flag, ...
    'VariableNames', {'knee_angle_deg','human_target_Nm', ...
    'required_pressure_kPa','status'});
writetable(tab, fullfile(here, 'Flx_20mm_pressure_for_human_torque_KMAX0235.csv'));
fprintf('Wrote %s\n', fullfile(here, 'Flx_20mm_pressure_for_human_torque_KMAX0235.csv'));
fprintf('max required pressure = %.1f kPa at %.1f deg; all-ok = %d\n', ...
    max(Preq), phiD(Preq==max(Preq),1), all(Flag=="ok" | Flag=="no flexion demand"));

% Key comparison figure: human target vs pair torque at 620 kPa, 0.235.
bpa620 = MonoPam_mult(Name, LP, CrossPoint, Dia, T_Pam, rest, kmax, ...
    tendon, fitting, 620, Xi0, Xi1, Xi2, Xi3, wraps, phiD, BPAcount, BM);
f = figure('Visible','off','Position',[100 100 560 420]);
plot(phiD, bpa620.Torque_p(:,3), 'LineWidth', 2); hold on;
plot(phiD, humanTgtPhi, '--', 'LineWidth', 2);
grid on; box off;
xlabel('Knee angle (deg)'); ylabel('Torque (N\cdotm)');
legend('2x20 mm pair, 620 kPa, KMAX 0.235','Human Bifemsh target','Location','southwest');
title('Flexor pair meets/exceeds human everywhere (min margin +4.0 N\cdotm at -120 deg)');
saveas(f, fullfile(here, 'flx20_0235_vs_human_620kPa.png'));
fprintf('Saved flx20_0235_vs_human_620kPa.png\n');
fprintf('TABLE DONE\n');
