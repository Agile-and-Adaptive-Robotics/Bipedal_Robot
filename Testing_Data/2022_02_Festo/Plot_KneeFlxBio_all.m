% Biomimetic-knee flexor comparison: 10 mm, 20 mm, and redesigned two-BPA route.
% Figure 4.11. Run this script to update the dissertation PDF, PNG, and FIG.
festodir = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(festodir));
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Dissertation/Figures/30-torque');
qcdir = outdir;
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Mesh_Optimization'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');

c3 = [0.980 0.529 0.459]; % light orange - Unoptimized
c4 = [0.918 0.373 0.580]; % pink - reduced-pressure measured points
c5 = [0.804 0.204 0.710]; % magenta - Corrected
c7 = [0 0 1];             % indigo - Measured
c8 = [0 0 0];             % black - human model
sz = 60;

%% adopted terms of record (row 77)
S = load('minimizeFlxPin10_results.mat');
g = S.filtered_results(77, S.xCols);
ledger = [3.937370e-03, 3.997873e+04, 1.473375e+04];
if any(abs(g(1:3) - ledger) > 1e-6 * ledger)
    error('Row 77 mismatch vs ledger');
end

%% panels A/B machinery (P4)
[~, ~, ~, bpab] = minimizeFlx(g(1), g(2), g(3));
[~, ~, ~, bpau] = minimizeFlx(0, Inf, Inf);
Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)'; Bifemsh_T = Tab(:,4)';

angle_data = readmatrix('Results_table_FullSize.xlsx','Sheet','FlxTest20mm_42cm (3)','Range','C7:X7');
pressure_data = readmatrix('Results_table_FullSize.xlsx','Sheet','FlxTest20mm_42cm (3)','Range','C8:X8');
torque_data = readmatrix('Results_table_FullSize.xlsx','Sheet','FlxTest20mm_42cm (3)','Range','C17:X17');
angle281 = angle_data(pressure_data == 281);
mask = ismember(pressure_data, [560, 421, 281]);
near325 = abs(pressure_data - 325) <= 5;
idx325 = find(near325 & angle_data < angle281, 1, 'first');
if ~isempty(idx325), mask(idx325) = true; end

%% panel C machinery (56 cm redesigned two-BPA route)
positions = 100;
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30 22.22 21.55 21.09 19.91 18.70 17.48 13.82 10.44 7.60 5.52 4.35 4.16 5.01 7.04 10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65 -417.03 -417.19 -417.28 -417.41 -417.41 -417.30 -416.28 -414.36 -411.72 -408.62 -405.32 -402.08 -399.16 -396.85 -395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

t1_ICR_x = ([29.66 28.54 27.86 27.40 26.23 25.03 23.81 20.03 16.17 12.34 8.67 5.24 2.04 -1.01 -4.1 -7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97 25.74 25.61 25.53 25.35 25.19 25.03 24.57 24.04 23.39 22.66 21.93 22.32 21.99 21.2 22.33]')/1000;
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

Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;
resultData = load(fullfile(root, 'Code/Matlab/Mesh_Optimization/Results', ...
    'Bifemsh_20mm_Result.mat'), 'xBest', 'routeCtx', 'XiUsed');
xBest = resultData.xBest;
routeCtx = resultData.routeCtx;
routeCtx.geo = buildGeoExclusion();
p1 = xBest(1:3);
p2 = xBest(4:6);
[Location, bendMeasure] = buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[LocationB, bendMeasureB] = buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);
LocationPair = {Location; LocationB};
bendMeasurePair = {bendMeasure; bendMeasureB};

Dia = 20;
rest = xBest(7);
tendon = xBest(8) - 0.02;
KMAX = 0.235;
kmax = (1 - KMAX)*rest;
fitting = 0.021;
pres3 = 620;

XiUsed = resultData.XiUsed;
Xi0 = XiUsed(1); Xi1 = XiUsed(2); Xi2 = XiUsed(3); Xi3 = XiUsed(4);
wraps = 6;
BPAcount = 2;

Bifemsh_Pam3 = MonoPam_mult(Name, LocationPair, CrossPoint, Dia, T_Pam, ...
    rest, kmax, tendon, fitting, pres3, Xi0, Xi1, Xi2, Xi3, wraps, phiD, ...
    BPAcount, bendMeasurePair);
TqCorr = Bifemsh_Pam3.Torque_p(:,3);
fprintf('Panel C pair torque at 620 kPa: corrected [min max] = [%.2f %.2f] N*m\n', ...
    min(TqCorr), max(TqCorr));

H = readmatrix('OpenSim_Bifem_Results.txt', 'FileType', 'text', 'NumHeaderLines', 7);
humanAngle = H(:,2);
TorqueHz = H(:,4);

% fifteen measured data measurements (56 cm flexor test, entered 2026-10-23)
Load   = [84 51 77.8 69 67.3 79 96 20 11.1 77 60 69.5 28 93.7 106];
K_deg  = [-61.5 -83 5 -5 -52 -45 -42 -114 -119 -62 -78 -87 -105 -22.5 -14];
LC_deg = [35 38 30 25 25 25 25 45 45 28.5 34 35 43 17 16];
Presh  = [517 628 226.5 242.5 406.7 410.3 412 634 633 400 489 635 638 275.3 273.7];
idx = Presh>610;

K_ang  = K_deg*pi/180;
LC_ang = LC_deg*pi/180;
p_rf = [0.095, -0.325, 0]';
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
    Fr(:,:,i) = [0; 0; 0; Load(i)*cos(LC_ang(i)+pi); Load(i)*sin(LC_ang(i)+pi); 0];
    AdTrk(:,:,i) = Adjoint(Trk(:,:,i));
    Fk(:,:,i) = AdTrk(:,:,i)'*Fr(:,:,i);
end
TorqueZ56 = squeeze(Fk(3,1,:));
TorqueZ56 = TorqueZ56(:);
tgt56 = interp1(humanAngle, TorqueHz, K_deg, 'pchip');
tgt56 = tgt56(:);
margin56 = abs(TorqueZ56) - abs(tgt56);
fprintf('Panel C: meet-or-exceed %d of %d; smallest margin %+.2f N*m at %.1f deg\n', ...
    sum(margin56 >= 0), numel(margin56), min(margin56), K_deg(margin56 == min(margin56)));

%% ---- figure: 2x2, three data panels + dedicated legend tile --------------
sty = {'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 15],'XTick',[-120 -90 -60 -30 0],'GridLineStyle','none'};

fig = figure('Name','FlxGrp_2x2_20261026','Color','w', ...
    'Units','centimeters','Position',[2 2 19.5 13.5]);
t = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
axT = gobjects(1,3);

% (A) 10 mm 41.5 cm, title '10 mm'
axT(1) = nexttile(1); hold on
i = 1;
hUnopt = plot(bpau(i).Ak, bpau(i).M_p(:,3), '--', 'Color', c3, 'LineWidth', 2, ...
    'HandleVisibility','off');
hCorr = plot(bpab(i).Ak, bpab(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, ...
    'HandleVisibility','off');
hMeas = scatter(bpab(i).Aexp, bpab(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
title('\bf 10 mm');
set(gca, sty{:});

% (B) 20 mm 41.5 cm + human + reduced-P points, title '20 mm'
axT(2) = nexttile(3); hold on
i = 2;
hHum = plot(knee_angle_rT, Bifemsh_T, ':', 'Color', c8, 'LineWidth', 4, ...
    'HandleVisibility','off');
plot(bpau(i).Ak, bpau(i).M_p(:,3), '--', 'Color', c3, 'LineWidth', 2, ...
    'HandleVisibility','off');
plot(bpab(i).Ak, bpab(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, ...
    'HandleVisibility','off');
scatter(bpab(i).Aexp, bpab(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
hRed = scatter(angle_data(mask), torque_data(mask), sz, 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'none', ...
    'HandleVisibility','off');
title('\bf 20 mm');
set(gca, sty{:});

% (C) 20 mm 56 cm redesigned two-BPA route, title '20 mm'
axT(3) = nexttile(4); hold on
plot(humanAngle, TorqueHz, ':', 'Color', c8, 'LineWidth', 4, ...
    'HandleVisibility','off');
plot(phiD, TqCorr, '-', 'Color', c5, 'LineWidth', 2.5, ...
    'HandleVisibility','off');
scatter(K_deg(idx), TorqueZ56(idx), sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
scatter(K_deg(~idx), TorqueZ56(~idx), sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c4, 'MarkerEdgeColor', 'none', 'HandleVisibility','off');
title('\bf 20 mm');
set(gca, sty{:});

% shared B/C y-range and x-range (both 20 mm vs the same human target);
% panel A keeps its own y-range (10 mm scale) -- no clipping to share.
tmpRP = torque_data(mask);
yBC = [bpab(2).M_p(:); bpau(2).M_p(:); tmpRP(:); ...
       Bifemsh_T(:); TorqueHz(:); TqCorr(:); TorqueZ56(:)];
yLo = floor(min(yBC)/10)*10 - 5;
axT(2).YLim = [yLo, 5];
axT(3).YLim = [yLo, 5];
tmpX = [phiD(:); K_deg(:); humanAngle(:); bpab(2).Ak(:); bpau(2).Ak(:)];
xLo = max(min(tmpX), -125);
xHi = min(max(tmpX), 10);
axT(2).XLim = [xLo, xHi];
axT(3).XLim = [xLo, xHi];

% tile 4: dedicated legend (dummy series with the panel styles)
axL = nexttile(2); hold on
plot(NaN, NaN, '--', 'Color', c3, 'LineWidth', 2, 'DisplayName','Unoptimized');
plot(NaN, NaN, '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Optimized');
plot(NaN, NaN, ':', 'Color', c8, 'LineWidth', 4, 'DisplayName','Human model (max)');
scatter(NaN, NaN, sz, 'filled', 'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', ...
    'DisplayName','Measured');
scatter(NaN, NaN, sz, 'filled', 'MarkerFaceColor', c4, 'MarkerEdgeColor', 'none', ...
    'DisplayName','Measured at reduced pressure');
lg = legend('Location','northwest'); lg.FontSize = 11; lg.Box = 'on';
axis(axL,'off');

axT(3).YTickLabelMode = 'auto';
ylabel(t, 'Torque, N\cdotm', 'Interpreter','tex'); xlabel(t, '\bf \theta_{\bf k} , \circ');

% Reserve a left margin for A/B outside the shared torque label.
t.OuterPosition = [0.07 0 0.93 1];
drawnow;
% A/B align vertically; C sits left of its y-axis numbers at the same height.
ann = {'(A)','(B)','(C)'};
for j = 1:3
    if j < 3, letterX = -0.42; else, letterX = -0.23; end
    text(axT(j),letterX,1.005,ann{j},'Units','normalized', ...
        'FontSize',12,'FontName','Arial','FontWeight','bold', ...
        'HorizontalAlignment','left','VerticalAlignment','bottom', ...
        'Clipping','off');
end
lg.FontSize = 12;
drawnow;

% Update the active dissertation asset and its editable MATLAB figure.
savefig(fig,fullfile(outdir,'FlxGrp_2x2_20261026_publication.fig'));
exportgraphics(fig,fullfile(outdir,'FlxGrp_2x2_20261026_publication.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outdir,'FlxGrp_2x2_20261026_publication.png'),'Resolution',300);
save(fullfile(festodir,'Plot_KneeFlxBio_all.mat'),'bpab','bpau','phiD','TqCorr','K_deg','TorqueZ56');
disp('Updated dissertation Figure 4.11');
