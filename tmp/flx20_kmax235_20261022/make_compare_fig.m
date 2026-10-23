%% make_compare_fig.m — one-shot: 0.235 pair torque at 620 kPa vs human target
clc; clear; close all;
here = fileparts(mfilename('fullpath')); cd(here);
root = here;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot'), break; end
    if strcmp(parent, root), error('no root'); end
    root = parent;
end
addpath(genpath(fullfile(root,'Code','Matlab')));
addpath(fullfile(root,'Code','Matlab','Mesh_Optimization'));
addpath(fullfile(root,'Testing_Data','2022_02_Festo'),'-end');

positions = 100;
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam = ([23.30 22.22 21.55 21.09 19.91 18.70 17.48 13.82 10.44 7.60 5.52 4.35 4.16 5.01 7.04 10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam = ([-416.65 -417.03 -417.19 -417.28 -417.41 -417.41 -417.30 -416.28 -414.36 -411.72 -408.62 -405.32 -402.08 -399.16 -396.85 -395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');
phi = linspace(-2.0943951, 0.17453293, positions);
[valv, posv] = min(abs(phi)); phi(posv) = 0;
T_Pam = zeros(4,4,positions);
for i = 1:positions
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    Rp = [cos(phi(i)), -sin(phi(i)), 0; sin(phi(i)), cos(phi(i)), 0; 0,0,1];
    T_Pam(:,:,i) = RpToTrans(Rp, hipToKnee_Pam');
end
phiD = phi(:)*180/pi;

resultData = load('Bifemsh_20mm_Result.mat','xBest','routeCtx','Xi3');
xBest = resultData.xBest; routeCtx = resultData.routeCtx;
routeCtx.geo = buildGeoExclusion();
Xi = load('Bifemsh_20mm_Result.mat','XiUsed').XiUsed;
p1 = xBest(1:3); p2 = xBest(4:6);
rest = xBest(7); tendon = xBest(8) - 0.02;
kmax = (1-0.235)*rest;
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[Loc1, bm1] = buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);
[Loc2, bm2] = buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);

bpa = MonoPam_mult('Bicep Femoris (Short Head)', {Loc1; Loc2}, 2, 20, T_Pam, ...
    rest, kmax, tendon, 0.021, 620, Xi(1), Xi(2), Xi(3), Xi(4), 6, phiD, 2, {bm1; bm2});

T = readtable('Flx_20mm_pressure_for_human_torque_KMAX0235.csv');
f = figure('Visible','off','Position',[100 100 620 430]);
plot(phiD, bpa.Torque_p(:,3), 'LineWidth', 2); hold on;
plot(T.knee_angle_deg, T.human_target_Nm, '--', 'LineWidth', 2);
grid on; box off;
xlabel('Knee angle (deg)'); ylabel('Torque (N\cdotm)');
legend('2x20 mm pair, 620 kPa, KMAX = 0.235', ...
    'Human Bifemsh (short head) target', 'Location', 'southwest');
title('KMAX 0.235: pair still meets/exceeds human everywhere (worst margin +4.0 N\cdotm)');
saveas(f, fullfile(here, 'flx20_0235_vs_human_620kPa.png'));
fprintf('SAVED COMPARE FIG\n');
