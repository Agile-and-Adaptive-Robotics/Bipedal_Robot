%% gen_Flx20mm_3panel_20261023.m
% Three separate panels for the knee-torque paper (dissertation reuses the
% same files), assembled in LaTeX with the subcaption package -- NOT one
% MATLAB figure.
%   A  Flx20mm_humanTarget       OpenSim maximum human isometric flexion
%                                torque (biceps femoris short head)
%   B  Flx20mm_rigidVsCorrected  original rigid-body prediction vs the
%                                chi0-chi3 corrected prediction vs measured
%   C  Flx20mm_meetsHuman        human target + corrected prediction +
%                                measured (meet-or-exceed story)
% Machinery mirrors Knee_Flexor_Data_20mm.m; the measured points are the
% 56 cm flexor-test arrays Ben entered there, with the same adjoint
% transform to the moving knee center (reaction point p_rf measured
% 2026-10-23; load-cell angle entered in the bracket frame, muscle force
% pi from that angle).
%
% Style of record: measured = indigo filled circles (size ~ pressure,
% alpha 0.75); corrected prediction = light magenta dashed; original
% rigid-body prediction = gray dash-dot; human target = black dotted.
% Knee angle negative in flexion; flexion torque negative.

clc; clear; close all;

set(groot, ...
    'defaultAxesFontName','Arial', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesFontWeight','bold', ...
    'defaultAxesLabelFontSizeMultiplier',1, ...
    'defaultAxesTitleFontSizeMultiplier',1.2, ...
    'defaultAxesLineWidth',2, ...
    'defaultAxesBox','off', ...
    'defaultAxesXMinorTick','on', ...
    'defaultAxesYMinorTick','on', ...
    'defaultAxesTickLength',[0.025 0.05], ...
    'defaultAxesXGrid','off', ...
    'defaultAxesYGrid','off', ...
    'defaultLineLineWidth',2, ...
    'defaultLegendFontName','Arial', ...
    'defaultLegendFontSize',8, ...
    'defaultLegendFontWeight','bold', ...
    'defaultLegendBox','off')

%% Paths (same block as Knee_Flexor_Data_20mm.m)
here = fileparts(mfilename('fullpath'));
root = here;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', here)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

outdir = fullfile(fileparts(here),'assets');

%% Joint kinematics (copied from Knee_Flexor_Data_20mm.m)
positions = 100;

knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	22.32	21.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);
[val, pos] = min(abs(phi));
phi(pos) = 0;

T_Pam = zeros(4, 4, positions);
T_t1_ICR = zeros(4, 4, positions);
for i = 1:positions
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam = [cos(phi(i)), -sin(phi(i)), 0;
             sin(phi(i)),  cos(phi(i)), 0;
             0, 0, 1];
    T_Pam(:, :, i) = RpToTrans(R_Pam, hipToKnee_Pam');
    t1toICR = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR');
end
phiD = phi*180/pi;

%% Design of record: flexor pair route (copied from Knee_Flexor_Data_20mm.m)
Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;

resultData = load('Bifemsh_20mm_Result.mat', ...
    'xBest', 'routeCtx', 'Xi3', 'XiUsed');
xBest = resultData.xBest;
routeCtx = resultData.routeCtx;
routeCtx.geo = buildGeoExclusion();

p1 = xBest(1:3);
p2 = xBest(4:6);

[Location, bendMeasure] = buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[LocationB, bendMeasureB] = ...
    buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);
LocationPair = {Location; LocationB};
bendMeasurePair = {bendMeasure; bendMeasureB};

Dia = 20;
rest = xBest(7);
tendon = xBest(8) - 0.02;
KMAX = 0.235;                      % maximum contraction fraction of record
kmax = (1 - KMAX)*rest;
fitting = 0.021;
pres3 = 620;

XiUsed = resultData.XiUsed;
Xi0 = XiUsed(1); Xi1 = XiUsed(2); Xi2 = XiUsed(3); Xi3 = XiUsed(4);
fprintf('Design of record: rest = %.4f m, tendon = %.4f m, Xi = [%.6g %.6g %.6g %.6g]\n', ...
    rest, tendon, Xi0, Xi1, Xi2, Xi3)

wraps = 6;
BPAcount = 2;

Bifemsh_Pam3 = MonoPam_mult(Name, LocationPair, CrossPoint, Dia, T_Pam, ...
    rest, kmax, tendon, fitting, pres3, Xi0, Xi1, Xi2, Xi3, wraps, phiD, ...
    BPAcount, bendMeasurePair);
TqCorr  = Bifemsh_Pam3.Torque_p(:,3);   % corrected (chi0-chi3) at 620 kPa
TqRigid = Bifemsh_Pam3.Torque(:,3);     % original rigid-body at 620 kPa
fprintf('Pair torque at 620 kPa (N*m): rigid [min max] = [%.2f %.2f], corrected [min max] = [%.2f %.2f]\n', ...
    min(TqRigid), max(TqRigid), min(TqCorr), max(TqCorr))
fprintf('NaN count: rigid %d, corrected %d\n', sum(isnan(TqRigid)), sum(isnan(TqCorr)))

%% Human target (OpenSim, biceps femoris short head)
H = readmatrix('OpenSim_Bifem_Results.txt', ...
    'FileType', 'text', ...
    'NumHeaderLines', 7);
humanAngle = H(:,2);
TorqueHz = H(:,4);

%% Measured points: 56 cm flexor test (Ben's arrays, entered 2026-10-23)
Load   = [84	51	77.8	69	67.3	79	96	20	11.1	77	60	69.5	28	93.7	106];   % Load (N)
K_deg  = [-61.5	-83	5	-5	-52	-45	-42	-114	-119	-62	-78	-87	-105	-22.5	-14];   % Knee angle, deg (flexion negative)
LC_deg = [35	38	30	25	25	25	25	45	45	28.5	34	35	43	17	16];   % Load-cell angle, deg (bracket frame)
Presh  = [517	628	226.5	242.5	406.7	410.3	412	634	633	400	489	635	638	275.3	273.7];   % Measured pressure, kPa

K_ang  = K_deg*pi/180;
LC_ang = LC_deg*pi/180;     % bracket-frame entry, muscle force +pi

p_rf = [0.095, -0.325, 0]'; % reaction force point, m (measured 2026-10-23)
d = hypot(p_rf(1), p_rf(2));
ang = atan2d(p_rf(2), p_rf(1));
T_t1_rf = RpToTrans(eye(3), p_rf);
Trk = pagemtimes(TransInv(T_t1_rf), T_t1_ICR);
fcn15 = fit(phi', squeeze(Trk(1,4,:)), 'cubicspline');
fcn16 = fit(phi', squeeze(Trk(2,4,:)), 'cubicspline');

Trk = zeros(4,4,length(Load));
Fr = zeros(6,1,length(Load));
AdTrk = zeros(6,6,length(Load));
Fk = zeros(6,1,length(Load));

for i = 1:length(Load)
    Trk(:,:,i) = RpToTrans(eye(3), [fcn15(K_ang(i)), fcn16(K_ang(i)), 0]');
    Fr(:,:,i) = [0; 0; 0; Load(i)*cos(LC_ang(i)+pi); ...
        Load(i)*sin(LC_ang(i)+pi); 0];
    AdTrk(:,:,i) = Adjoint(Trk(:,:,i));
    Fk(:,:,i) = AdTrk(:,:,i)'*Fr(:,:,i);
end

TorqueZ56 = squeeze(Fk(3,1,:));

fprintf('\nFlxTest20mm_56cm adjoint torques (d = %.4f m, ang = %.2f deg):\n', d, ang);
TorqueZ56 = TorqueZ56(:);
tgt56 = interp1(humanAngle, TorqueHz, K_deg, 'pchip');
tgt56 = tgt56(:);
margin56 = abs(TorqueZ56) - abs(tgt56);
fprintf('angle_deg   load_N   press_kPa   meas_Nm   human_Nm   margin_Nm\n');
for i = 1:length(Load)
    fprintf('%8.1f %9.1f %10.1f %10.2f %10.2f %11.2f\n', ...
        K_deg(i), Load(i), Presh(i), TorqueZ56(i), tgt56(i), margin56(i));
end
[mnM, imnM] = min(margin56);
fprintf('Meet-or-exceed count: %d of %d points; smallest margin %+.2f N*m at %.1f deg\n', ...
    sum(margin56 >= 0), numel(margin56), mnM, K_deg(imnM));
if any(Presh > 620)
    fprintf('NOTE: %d points recorded above 620 kPa (max %.1f kPa).\n', ...
        sum(Presh > 620), max(Presh));
end

%% Defense layout only: retain all calculations above, enlarge plots and move legends.
outdir = fullfile(fileparts(here),'assets');
set(groot,'defaultAxesFontSize',16,'defaultLegendFontSize',12);
for panel=1:3
 f=figure('Visible','off','Units','centimeters','Position',[2 2 25 12],'Color','w'); hold on
 if panel==1
  plot(humanAngle,TorqueHz,':','Color','k','LineWidth',2.5); names={'Human target'}; key='human_final';
 elseif panel==2
  plot(phiD,TqRigid,'-.','Color','#8C8C8C','LineWidth',2.5);
  plot(phiD,TqCorr,'--','Color','#EA5F94','LineWidth',2.5);
  scatter(K_deg,TorqueZ56,60*Presh/620,'filled','MarkerFaceColor','#0000FF','MarkerFaceAlpha',0.75);
  names={'Rigid-body, 620 kPa','Corrected, 620 kPa','Measured, variable pressure'};key='rigid_final';
 else
  plot(humanAngle,TorqueHz,':','Color','k','LineWidth',2.5);
  plot(phiD,TqCorr,'--','Color','#EA5F94','LineWidth',2.5);
  scatter(K_deg,TorqueZ56,60*Presh/620,'filled','MarkerFaceColor','#0000FF','MarkerFaceAlpha',0.75);
  names={'Human target','Corrected, 620 kPa','Measured, variable pressure'};key='meets_final';
 end
 xlim([-120 10]);xticks([-120 -90 -60 -30 0]);
 xlabel('Knee angle, \circ');ylabel('Torque, N \cdot m');
 legend(names,'Location','northoutside','NumColumns',1);box off
 exportgraphics(f,fullfile(outdir,[key '.png']),'Resolution',200);
 exportgraphics(f,fullfile(outdir,[key '.pdf']),'ContentType','vector');close(f)
end
T=table(K_deg(:),Load(:),Presh(:),TorqueZ56,tgt56,margin56,'VariableNames',{'angle_deg','load_N','pressure_kPa','measured_Nm','human_Nm','margin_Nm'});
writetable(T,fullfile(outdir,'flexor_final_values.csv'));

