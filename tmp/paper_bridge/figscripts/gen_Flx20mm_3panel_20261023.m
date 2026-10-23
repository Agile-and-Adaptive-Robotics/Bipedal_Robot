%% gen_Flx20mm_3panel_20261023.m
% ONE MATLAB figure, three panels across, panel letters (A)(B)(C) drawn in
% per Ben's FlxPin_group.fig margin-left convention (constants copied from
% gen_ExtPin6Panel.m). Counts as ONE figure float against the Frontiers
% combined figure+table limit; no LaTeX subfigure machinery.
%   (A) OpenSim maximum human isometric flexion torque (bifemsh short head)
%   (B) original rigid-body prediction vs chi0-chi3 corrected prediction
%       vs measured points
%   (C) human target + corrected prediction + measured (meet-or-exceed)
% Machinery mirrors Knee_Flexor_Data_20mm.m; the measured points are the
% 56 cm flexor-test arrays Ben entered there, same adjoint transform.
%
% Style of record: measured = indigo filled circles (size ~ pressure,
% alpha 0.75); corrected prediction = magenta dashed; original rigid-body
% prediction = gray dash-dot; human target = black dotted.
% Output: Knee_Torque_Test/Figures/Flx20mm_56cm.{eps,pdf} (single figure,
% paper naming convention dia + resting length) + QC png here.

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

outdir = fullfile(root, 'Documentation', 'Reports and Papers', ...
    'Knee_Torque_Test', 'Figures');

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

%% One figure, three panels, letters drawn in (FlxPin_group convention)
clrMeas  = '#0000FF';   % indigo        -- measured
clrCorr  = '#EA5F94';   % light magenta -- corrected prediction
clrRigid = '#8C8C8C';   % gray          -- original rigid-body prediction
clrHuman = '#000000';   % black         -- human target

XL = [-120 10];
XT = [-120 -90 -60 -30 0];
% Per-panel YLim: the rigid curve reaches -190 N*m and would flatten the
% human/measured detail in the other panels, so each panel fits its own
% series (captions and prose carry the numbers).
YLA = [floor(min(TorqueHz)) - 1, 2];
YLB = [floor(min([TqRigid; TqCorr; TorqueZ56])) - 2, 5];
YLC = [floor(min([TqCorr; TorqueHz; TorqueZ56])) - 2, 5];
fprintf('\nShared XLim [%d %d]; YLim A [%.0f %.0f], B [%.0f %.0f], C [%.0f %.0f]\n', ...
    XL, YLA, YLB, YLC);

szM = 60*Presh/620;      % measured marker size ~ pressure

fig = figure('Name','Flx20mm_56cm','Color','w', ...
    'Units','centimeters','Position',[2 2 17.2 5.4]);
t = tiledlayout(1,3,'TileSpacing','loose','Padding','loose');
axT = gobjects(1,3);

% ---- (A) human target ---------------------------------------------------
axT(1) = nexttile(1); hold on
plot(humanAngle, TorqueHz, ':', 'Color', clrHuman)
xlim(XL); xticks(XT); ylim(YLA);
legend('Human target','Location','northeast')

% ---- (B) rigid vs corrected vs measured ---------------------------------
axT(2) = nexttile(2); hold on
plot(phiD, TqRigid, '-.', 'Color', clrRigid)
plot(phiD, TqCorr, '--', 'Color', clrCorr)
scatter(K_deg, TorqueZ56, szM, 'filled', ...
    'MarkerFaceColor', clrMeas, 'MarkerFaceAlpha', 0.75, ...
    'MarkerEdgeColor', 'none')
xlim(XL); xticks(XT); ylim(YLB);
legend('Rigid-body prediction', 'Corrected prediction', ...
    'Measured','Location','northeast')

% ---- (C) human vs corrected vs measured ---------------------------------
axT(3) = nexttile(3); hold on
plot(humanAngle, TorqueHz, ':', 'Color', clrHuman)
plot(phiD, TqCorr, '--', 'Color', clrCorr)
scatter(K_deg, TorqueZ56, szM, 'filled', ...
    'MarkerFaceColor', clrMeas, 'MarkerFaceAlpha', 0.75, ...
    'MarkerEdgeColor', 'none')
xlim(XL); xticks(XT); ylim(YLC);

ylabel(t, '\bf Torque, N\cdotm');
xlabel(t, '\bf \theta_{k} , \circ');

% panel letters per Ben's FlxPin_group.fig convention (constants from
% gen_ExtPin6Panel.m): margin labels LEFT of each tile's outer box,
% straddling the plot-box top; (A) sits at the opened left margin
% clear of the rotated y-label band.
ann = {'(A)','(B)','(C)'};
drawnow;
ylL = 0.008; ylR = 0.024;          % rotated y-label band (figure-normalized)
xL1 = 0.040;                       % (A) target left edge
if (xL1 < ylR) && (xL1 + 0.032 > ylL)   % open the margin if y-label is in the way
    S = (xL1 + 0.032 + 0.008) - ylL;
    tp = t.Position; t.Position = [tp(1) + S, tp(2), tp(3) - S, tp(4)];
    drawnow;
    ylL = ylL + S; ylR = ylR + S;
end
fprintf('AUDIT margin: yLabel band x[%.4f %.4f]  (A) left=%.4f  gap=%.4f\n', ...
    ylL, ylR, xL1, xL1 - ylR);

apos = zeros(3,2); outL = zeros(3,1);
for j = 1:3
    pos0 = axT(j).Position;
    outL(j) = axT(j).OuterPosition(1);
    tcY = pos0(2) + pos0(4);           % plot-box top: letter straddles it
    if j == 1, xj = xL1; else, xj = outL(j) - 0.037; end
    a = annotation(fig,'textbox',[xj, tcY - 0.020, 0.032, 0.04], ...
        'String',['\bf ' ann{j}], 'FontSize',12,'FontName','Arial', ...
        'FontWeight','bold','EdgeColor','none','HorizontalAlignment','left', ...
        'VerticalAlignment','top');
    try
        drawnow;
        ae = a.Extent;                 % center the letter on the box top
        a.Position(2) = a.Position(2) + (tcY - (ae(2) + ae(4)/2));
    catch
    end
    apos(j,:) = a.Position(1:2);
end

% ---- audit of the drawn state (verify, never assume) --------------------
lgAll = findall(fig, 'Type', 'legend');
fprintf('AUDIT legends: count=%d locations=%s\n', numel(lgAll), ...
    strjoin(arrayfun(@(L) L.Location, lgAll, 'UniformOutput', false), ','));
for j = 1:3
    if j == 1
        anchored = abs(apos(j,1) - xL1) < 1e-9;
    else
        anchored = abs(apos(j,1) + 0.032 - (outL(j) - 0.005)) < 5e-3;
    end
    fprintf('AUDIT tile %d %s: letter@(%.3f,%.3f) outerL=%.3f anchored=%d\n', ...
        j, ann{j}, apos(j,1), apos(j,2), outL(j), anchored);
end
% letter boxes clear of data (data coords per tile)
hAx = humanAngle(:).';  hTq = TorqueHz(:).';          % rows throughout
pR = TqRigid.';  pC = TqCorr.';  mT = TorqueZ56.';
series = {hAx, hTq; ...
    [phiD, phiD, K_deg], [pR, pC, mT]; ...
    [hAx, phiD, K_deg], [hTq, pC, mT]};
for j = 1:3
    pos0 = axT(j).Position;
    bx = axT(j).XLim(1) + ([apos(j,1), apos(j,1)+0.032] - pos0(1)) / pos0(3) * diff(axT(j).XLim);
    by = axT(j).YLim(1) + ([apos(j,2), apos(j,2)+0.04] - pos0(2)) / pos0(4) * diff(axT(j).YLim);
    xs = series{j,1}; ys = series{j,2};
    n_ov = sum(xs > bx(1) & xs < bx(2) & ys > by(1) & ys < by(2));
    fprintf('AUDIT tile %d letter data-box x[%.0f %.0f] y[%.1f %.1f] points-under-letter=%d\n', ...
        j, bx(1), bx(2), by(1), by(2), n_ov);
end
% legends in-figure and clear of data
for L = 1:numel(lgAll)
    p = lgAll(L).Position;
    inside = p(1) >= -0.001 && p(2) >= -0.001 && p(1)+p(3) <= 1.001 && p(2)+p(4) <= 1.001;
    fprintf('AUDIT legend %d in-figure=%d pos=[%.3f %.3f %.3f %.3f]\n', L, inside, p);
end

%% Exports: single figure, EPS (paper) + PDF (dissertation) + QC png
exportgraphics(fig, fullfile(outdir, 'Flx20mm_56cm.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(outdir, 'Flx20mm_56cm.pdf'), 'ContentType','vector');
exportgraphics(fig, fullfile(here, 'QC_Flx20mm_56cm.png'), 'Resolution', 300);
fprintf('Wrote Flx20mm_56cm.{eps,pdf} + QC png\n');
fprintf('\n[DONE gen_Flx20mm_3panel]\n');
