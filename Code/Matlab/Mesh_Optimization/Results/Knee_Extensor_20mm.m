%% Knee extensor data and optimized 20 mm BPA comparison
clc
clear
close all

%% Repo root and path setup (derived from this file's location)
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end

addpath(genpath(fullfile(root, 'Code', 'Matlab')));
% Mesh_Optimization must win any shadowing contest against data subfolders.
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));

% buildKneeExtContext20mm reads OpenSim_Vasti_Results.txt and its Xi-pick
% mat by bare name; both live in Testing_Data\2022_02_Festo. Append (do not
% prepend) so Code\Matlab keeps winning any name collisions.
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

%% Load the Opt_run_Ext result and rebuild its context
% ctx is a struct built by buildKneeExtContext20mm(), called from both
% Opt_run_Ext.m and Opt_sanity_Ext.m. It stores:
%   N, phi, phiD, pos      sample count, radians/degrees, zero-angle index
%   T                     human tibia-frame -> femur transforms
%   T_Pam, T_Pam_inv       robot ICR <-> femur transforms
%   T_ICR_t1, T_t1_ICR     robot t1 <-> ICR transforms
%   T_t1_f, T_f_t1         robot femur <-> t1 composite transforms
%   routeSeed, geo         native 9-point seed and extensor CAD geometry
%   CrossPoint, routeRows  frame split (row 6) and route size (9 rows)
%   Dia, BPAcount, fitting, wraps, targetPressure, KMAX, Xi0:Xi3
%                         BPA/model constants and identified parameters
%   humanAngleD, humanTorque, targetName  SIGNED selected OpenSim target
%   x0, lb, ub             original design and optimization bounds
%   initial*              original-route/length diagnostic values
% The extensor builder constructs geo locally and uses tendonLimit20mm,
% buildDistalRingLocation20mm and its local muscleLengthNormal helper to
% initialize the original design. It does not use buildGeoExclusion.
%
% Dated Opt_run_Ext results (2026-09-10 route-elimination rework onward)
% save xBest and the XiUsed record but NOT ctx. Rebuild the context fresh
% and require its Xi block to equal the run's XiUsed; a mismatch means the
% builder's Xi pick moved on since the run and the display would mislabel
% the design. Do not load old Location/bendMeasure arrays or expect
% separate saved endpoint variables.
resultFile = fullfile(scriptDir, 'Vas_Pam_20mm_Result.mat');
S = load(resultFile, 'xBest', 'XiUsed');
if ~isfield(S,'xBest') || ~isfield(S,'XiUsed')
    error('Knee_Extensor_20mm:MissingResult', ...
        '%s must contain xBest and XiUsed from Opt_run_Ext.', resultFile)
end

ctx = buildKneeExtContext20mm();

xiBuilt = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
if ~isequal(xiBuilt, S.XiUsed)
    error('Knee_Extensor_20mm:XiMismatch', ...
        'Builder Xi [%.6g %.6g %.6g %.6g] does not equal run XiUsed [%.6g %.6g %.6g %.6g]. Update buildKneeExtContext20mm or repoint resultFile.', ...
        xiBuilt, S.XiUsed)
end

xBest = reshape(S.xBest,1,[]);
validateattributes(xBest,{'numeric'},{'real','finite','numel',8});
p1 = xBest(1:3);
pEnd = xBest(4:6);
rest = xBest(7);
tendon = xBest(8);
KMAX = ctx.KMAX;
kmax = rest*(1-KMAX);

%% Hand-adjustment block (Ben): override the route seed rows here.
% Leave empty to use the builder's seeds (ctx.routeSeed, set by the
% 2026-09-21 fcec-ray rule in buildKneeExtContext20mm). To hand-shape a
% route, enter ALL SEVEN rows p2:p8 as [x y] pairs -- FEMUR frame for
% p2:p5, T1 frame for p6:p8 (z is redistributed along the route
% automatically). Example:
% handSeed = [ ...
%     0.0839, -0.2748; ...   % p2
%     0.0653, -0.3880; ...   % p3
%     0.0641, -0.4268; ...   % p4
%     0.0378, -0.4521; ...   % p5
%     0.0699,  0.0343; ...   % p6
%     0.0735,  0.0191; ...   % p7
%     0.0720, -0.0122];      % p8
handSeed = [];
if ~isempty(handSeed)
    ctx.routeSeed(2:8,1:2) = handSeed;
end

% Rebuild the same route as predictKneeExt20mm, for every knee position.
% The first/last NATIVE rows are p1 and pEnd from xBest. The builder supplies
% p2:p8 and repeated rows for eliminated points, then converts tibia-side
% rows to ICR coordinates. Do not overwrite Location(end,:,:) with raw pEnd.
[Location,bendMeasure,routeInfo] = ...
    buildDistalRingLocation20mm(p1,pEnd,tendon,ctx);

%% Hardcoded release schedule (Ben, 2026-09-21)
% Actual point-release angles for this design (Vas_Pam_20mm_Result.mat
% = the 2026-09-20 15:19 xBest) under Ben's restated rule: every
% optional row tested at every pose between its nearest ACTIVE
% neighbors (repeating rows + frame transforms, +90/-90 deg rotations),
% NO cascade wait, NO p7-triplet anchor, bypass gates at the verified
% 3 mm relaxations, p7 stays (real wall), p3:p5 by the fcec-ray rule.
% Each angle is the first swept pose at which the row is OFF (pose
% spacing 1.31 deg); p1/p2/p9 are always active. The route above still
% derives the schedule at runtime; the loop below reports any drift
% (hand-tuning xBest or handSeed, or moving a geo knob, will show here
% and these constants need re-recording).
% KNOWN SEED SIGNALS under this state (rule open but bypass chord
% collides = seed placement to improve, not the gate): p3 blocked only
% -4.44..+4.75 deg (worst margin +0.51 deg, near the hysteresis noise);
% p8 blocked -120..+6.06 deg (hysteresis range, worst -0.62 deg).
releaseRows   = [5;  4;  6;  7;  3;  8];
releaseAngleD = [-89.797979628815; -55.656565451176; ...
                 -47.777777564029; -20.202019959013; ...
                 6.060606331479; 7.373737646003];   % deg
for iRel = 1:numel(releaseRows)
    row = releaseRows(iRel);
    trans = find(routeInfo.active(row,1:end-1) ...
        & ~routeInfo.active(row,2:end), 1);
    if isempty(trans)
        derivedD = NaN;   % never released (inactive from the start)
    else
        derivedD = ctx.phiD(trans+1);
    end
    fprintf('p%d release: hardcoded %+.2f deg, derived %+.2f deg\n', ...
        row, releaseAngleD(iRel), derivedD);
end
fprintf('p7 active anywhere (hardcoded: yes until wall lift-off): %d\n', ...
    any(routeInfo.active(7,:)));

positions = ctx.N;
phi = ctx.phi;
phiD = ctx.phiD;
T = ctx.T;
T_Pam = ctx.T_Pam;
T_ICR_t1 = ctx.T_ICR_t1;
T_t1_ICR = ctx.T_t1_ICR;
pos = ctx.pos;
if phi(pos) ~= 0
    error('Knee_Extensor_20mm:MissingZeroPose', ...
        'ctx.pos must select the zero-angle pose used to build the transforms.')
end

fprintf('The algorithm will be calculating Torque at %d different joint positions.\n',positions)

%% Human quadriceps muscle paths
cdeg = pi/180;
rect_fem_x = [0.0156367;0.0179948;0.0274274;0.029683;0.0306;0.0366;0.0422;0.0451;0.0484;0.0533;0.0617;0.0634;0.067;0.0733];
rect_fem_xD = cdeg*[-120.118;-114.871;-90.068;-83.532;-80;-60;-40;-30;-20;-10;0;1.6;5;10];
rect_fem_y = [0.0234;0.0238;0.0251;0.0253;0.025284;0.0249;0.0243;0.0239;0.0234;0.0228;0.0210;0.0206;0.0192;0.0160];
rect_fem_yD = cdeg*[-120;-114.6;-90;-83.5;-80.01;-60;-40;-30;-20;-10;0;1.6;5;10];

fcn5 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.009811),'smoothingspline');
fcn6 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.02242),'smoothingspline');
fcn7 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.008285),'smoothingspline');
fcn8 = fit(rect_fem_yD,rect_fem_y+(0.0256239-0.02346),'smoothingspline');
fcn9 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.0142881),'smoothingspline');
fcn10 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.0215281),'smoothingspline');

Name = 'Vastus Medialis';
MIF = 1294;
OFL = 0.089;
TSL = 0.126;
Pennation = 0.08726646;
HumanLocation = zeros(5,3,positions);

for i = 1:positions
    if phiD(i) < -101
        HumanLocation(:,:,i) = [0.014,-0.210,0.019; ...
                                0.036,-0.277,0.001; ...
                                0.037,-0.405,-0.013; ...
                                0.027,-0.425,-0.013; ...
                                fcn5(phi(i)),fcn6(phi(i)),-0.0146];
    elseif phiD(i) < -69 && phiD(i) >= -101
        HumanLocation(:,:,i) = [0.014,-0.210,0.019; ...
                                0.036,-0.277,0.001; ...
                                0.037,-0.405,-0.013; ...
                                0.037,-0.405,-0.013; ...
                                fcn5(phi(i)),fcn6(phi(i)),-0.0146];
    else
        HumanLocation(:,:,i) = [0.014,-0.210,0.019; ...
                                0.036,-0.277,0.001; ...
                                0.036,-0.277,0.001; ...
                                0.036,-0.277,0.001; ...
                                fcn5(phi(i)),fcn6(phi(i)),-0.0146];
    end
end

Vas_Med = MonoMuscleData(Name,HumanLocation,5,MIF,TSL,Pennation,OFL,T);

Name = 'Vastus Intermedius';
MIF = 1365;
OFL = 0.087;
TSL = 0.136;
Pennation = 0.05235988;
HumanLocation = zeros(4,3,positions);

for i = 1:positions
    if phiD(i) < -80
        HumanLocation(:,:,i) = [0.029,-0.192,0.031; ...
                                0.034,-0.208,0.029; ...
                                0.034,-0.403,0.005; ...
                                fcn7(phi(i)),fcn8(phi(i)),0.0018];
    else
        HumanLocation(:,:,i) = [0.029,-0.192,0.031; ...
                                0.034,-0.208,0.029; ...
                                0.034,-0.208,0.029; ...
                                fcn7(phi(i)),fcn8(phi(i)),0.0018];
    end
end

Vas_Int = MonoMuscleData(Name,HumanLocation,4,MIF,TSL,Pennation,OFL,T);

Name = 'Vastus Lateralis';
MIF = 1871;
OFL = 0.084;
TSL = 0.157;
Pennation = 0.08726646;
HumanLocation = zeros(5,3,positions);

for i = 1:positions
    if phiD(i) < -110
        HumanLocation(:,:,i) = [0.005,-0.185,0.035; ...
                                0.027,-0.259,0.041; ...
                                0.036,-0.403,0.021; ...
                                0.025,-0.424,0.018; ...
                                fcn9(phi(i)),fcn10(phi(i)),0.0165];
    elseif phiD(i) < -69 && phiD(i) >= -110
        HumanLocation(:,:,i) = [0.005,-0.185,0.035; ...
                                0.027,-0.259,0.041; ...
                                0.036,-0.403,0.021; ...
                                0.036,-0.403,0.021; ...
                                fcn9(phi(i)),fcn10(phi(i)),0.0165];
    else
        HumanLocation(:,:,i) = [0.005,-0.185,0.035; ...
                                0.027,-0.259,0.041; ...
                                0.027,-0.259,0.041; ...
                                0.027,-0.259,0.041; ...
                                fcn9(phi(i)),fcn10(phi(i)),0.0165];
    end
end

Vas_Lat = MonoMuscleData(Name,HumanLocation,5,MIF,TSL,Pennation,OFL,T);

%% Original and optimized X3 BPA models
Name = 'Vastus Medialis Proximal Ring 20mm BPA';
CrossPoint = ctx.CrossPoint;
Dia = ctx.Dia;
fitting = ctx.fitting;
pres = ctx.targetPressure;
Xi0 = ctx.Xi0;
Xi1 = ctx.Xi1;
Xi2 = ctx.Xi2;
Xi3 = ctx.Xi3;
wraps = ctx.wraps;
BPAcount = ctx.BPAcount;

p10 = ctx.x0(1:3);
pEnd0 = ctx.x0(4:6);
rest0 = ctx.x0(7);
tendon0 = ctx.x0(8);
kmax0 = (1-KMAX)*rest0;
[Location0,bendMeasure0] = buildDistalRingLocation20mm(p10,pEnd0,tendon0,ctx);

Vas_Pam0 = MonoPamDataExplicit_balanceX3(Name,Location0,CrossPoint,Dia,T_Pam,rest0,kmax0,tendon0,fitting,pres,Xi0,Xi1,Xi2,Xi3,wraps,phiD,BPAcount,bendMeasure0);
Vas_Pam2 = MonoPamDataExplicit_balanceX3(Name,Location,CrossPoint,Dia,T_Pam,rest,kmax,tendon,fitting,200,Xi0,Xi1,Xi2,Xi3,wraps,phiD,BPAcount,bendMeasure);
Vas_Pam3 = MonoPamDataExplicit_balanceX3(Name,Location,CrossPoint,Dia,T_Pam,rest,kmax,tendon,fitting,pres,Xi0,Xi1,Xi2,Xi3,wraps,phiD,BPAcount,bendMeasure);

Torque0 = Vas_Pam0.Torque_p(:,3);
Torque2 = Vas_Pam2.Torque_p(:,3);
Torque3 = Vas_Pam3.Torque_p(:,3);
Lm0 = rest0.*(1-Vas_Pam0.strain_p(:));
Lm2 = rest.*(1-Vas_Pam2.strain_p(:));
Lm3 = rest.*(1-Vas_Pam3.strain_p(:));
relstrain0 = Vas_Pam0.strain_f(:)/KMAX;
relstrain2 = Vas_Pam2.strain_f(:)/KMAX;
relstrain3 = Vas_Pam3.strain_f(:)/KMAX;
G0 = hypot(Vas_Pam0.mA_p(:,1),Vas_Pam0.mA_p(:,2));
G2 = hypot(Vas_Pam2.mA_p(:,1),Vas_Pam2.mA_p(:,2));
G3 = hypot(Vas_Pam3.mA_p(:,1),Vas_Pam3.mA_p(:,2));

% Preserve the signed OpenSim torque; do not replace it by its magnitude.
humanTorque = interp1(ctx.humanAngleD,ctx.humanTorque,phiD,'pchip','extrap');
humanTorque = humanTorque(:);
validHumanTorque = isfinite(humanTorque) & humanTorque ~= 0;
torqueMarginFraction = nan(size(Torque3));
torqueMarginFraction(validHumanTorque) = ...
    Torque3(validHumanTorque)./humanTorque(validHumanTorque)-1;

fprintf('\nLoaded optimized extensor design:\n')
fprintf('p1     = [%.6f %.6f %.6f] m\n',p1)
fprintf('pEnd   = [%.6f %.6f %.6f] m\n',pEnd)
fprintf('rest   = %.6f m\n',rest)
fprintf('tendon = %.6f m\n',tendon)
fprintf('kmax   = %.6f m\n',kmax)
fprintf('minimum torque margin = %+.6f (%+.2f%%)\n',min(torqueMarginFraction(validHumanTorque)),100*min(torqueMarginFraction(validHumanTorque)))

%% Plot settings
run('Colors.m')
originalColor = [0.4 0.4 0.4];
optimizedColor = c{6};
optClr2 = c{5};
humanColor = '#000000';
fontName = 'Arial';
axesFontSize = 10;
titleFontSize = 12;
legendFontSize = 8;
tickLength = [0.025 0.05];
figurePosition = [2 2 14 10.5];
xLimits = [min(phiD),max(phiD)];

%% Extensor torque
figure('Name','Extensor Torque', 'Color','w', 'Units','centimeters', 'Position',figurePosition)
ax = gca;
hold(ax,'on')
plot(ax,phiD,Torque0,'--','Color',originalColor,'LineWidth',2,'DisplayName','Original BPA')
plot(ax,phiD,Torque2,'--','Color',optClr2,'LineWidth',2,'DisplayName','Optimized BPA, 200 kPa')
plot(ax,phiD,Torque3,'-','Color',optimizedColor,'LineWidth',2.5,'DisplayName','Optimized BPA, 620 kPa')
plot(ax,ctx.humanAngleD,ctx.humanTorque,':','Color',humanColor,'LineWidth',4,'DisplayName','Human target')
formatAxes(ax,fontName,axesFontSize,tickLength,xLimits)
% Fit the signed torque data; the old [0 15] limit clipped these predictions.
torqueLimits = [min([0;Torque0(:);Torque2(:);Torque3(:);ctx.humanTorque(:)]), ...
               max([0;Torque0(:);Torque2(:);Torque3(:);ctx.humanTorque(:)])];
torquePadding = 0.05*diff(torqueLimits);
if torquePadding == 0
    torquePadding = 1;
end
ylim(ax,torqueLimits + [-torquePadding,torquePadding])
xlabel(ax,'\theta_k, °','Interpreter','tex','FontWeight','bold')
ylabel(ax,'Torque, N\cdotm','Interpreter','tex','FontWeight','bold')
title(ax,'Extensor Torque','FontName',fontName,'FontSize',titleFontSize,'FontWeight','bold')
formatLegend(ax,fontName,legendFontSize)

%% Muscle length
figure('Name','Muscle Length', 'Color','w', 'Units','centimeters', 'Position',figurePosition)
ax = gca;
hold(ax,'on')
plot(ax,phiD,Lm0,'--','Color',originalColor,'LineWidth',2,'DisplayName','Original BPA')
plot(ax,phiD,Lm2,'-','Color',optClr2,'LineWidth',2.5,'DisplayName','Optimized BPA, 200 kPa')
plot(ax,phiD,Lm3,'-','Color',optimizedColor,'LineWidth',2.5,'DisplayName','Optimized BPA, 620 kPa')
formatAxes(ax,fontName,axesFontSize,tickLength,xLimits)
xlabel(ax,'\theta_k, °','Interpreter','tex','FontWeight','bold')
ylabel(ax,'Muscle Length, m','FontWeight','bold')
title(ax,'Muscle Length, L_m','Interpreter','tex','FontName',fontName,'FontSize',titleFontSize,'FontWeight','bold')
formatLegend(ax,fontName,legendFontSize)

%% Effective relative strain including Xi3
figure('Name','Relative Strain', 'Color','w', 'Units','centimeters', 'Position',figurePosition)
ax = gca;
hold(ax,'on')
plot(ax,phiD,relstrain0,'--','Color',originalColor,'LineWidth',2,'DisplayName','Original BPA')
plot(ax,phiD,relstrain2,'-','Color',optClr2,'LineWidth',2.5,'DisplayName','Optimized BPA')
plot(ax,phiD,relstrain3,'-','Color',optimizedColor,'LineWidth',2.5,'DisplayName','Optimized BPA')
formatAxes(ax,fontName,axesFontSize,tickLength,xLimits)
xlabel(ax,'\theta_k, °','Interpreter','tex','FontWeight','bold')
ylabel(ax,[char(949) '^{*}'], 'Interpreter','tex','FontName','Arial','FontWeight','bold');
title(ax,'Relative Strain','FontName',fontName,'FontSize',titleFontSize,'FontWeight','bold')
formatLegend(ax,fontName,legendFontSize)

%% Moment arm
figure('Name','Moment Arm', 'Color','w', 'Units','centimeters', 'Position',figurePosition)
ax = gca;
hold(ax,'on')
plot(ax,phiD,G0,'--','Color',originalColor,'LineWidth',2,'DisplayName','Original BPA')
plot(ax,phiD,G2,'-','Color',optClr2,'LineWidth',2.5,'DisplayName','Optimized BPA, 200 kPa')
plot(ax,phiD,G3,'-','Color',optimizedColor,'LineWidth',2.5,'DisplayName','Optimized BPA, 620 kPa')
formatAxes(ax,fontName,axesFontSize,tickLength,xLimits)
xlabel(ax,'\theta_k, °','Interpreter','tex','FontWeight','bold')
ylabel(ax,'Moment Arm, m','FontWeight','bold')
title(ax,'Moment Arm','FontName',fontName,'FontSize',titleFontSize,'FontWeight','bold')
formatLegend(ax,fontName,legendFontSize)

%% Torque margin
figure('Name','Torque Margin Fraction', 'Color','w', 'Units','centimeters', 'Position',figurePosition)
ax = gca;
hold(ax,'on')
plot(ax,phiD,100*torqueMarginFraction,'-','Color',optimizedColor,'LineWidth',2.5,'DisplayName','Optimized BPA')
% The extensor objective originally has no extra torque-margin field.
requiredMargin = 0;
if isfield(ctx,'requiredTorqueMargin')
    requiredMargin = ctx.requiredTorqueMargin;
elseif isfield(ctx,'targetMargin')
    requiredMargin = ctx.targetMargin;
end
requiredLine = yline(ax,100*requiredMargin,':','Required margin','Color',humanColor,'LineWidth',2,'HandleVisibility','off');
requiredLine.FontName = fontName;
requiredLine.FontSize = legendFontSize;
requiredLine.FontWeight = 'bold';
formatAxes(ax,fontName,axesFontSize,tickLength,xLimits)
xlabel(ax,'\theta_k, °','FontWeight','bold')
ylabel(ax,'Torque Margin, %','FontWeight','bold')
title(ax,'BPA Torque Margin Relative to Human','FontName',fontName,'FontSize',titleFontSize,'FontWeight','bold')
formatLegend(ax,fontName,legendFontSize)

%% Plot full optimized geometry and p1:p9 route
% Port of Opt_run_Ext's route-geometry figure (Ben, 2026-09-21: "so I
% can see if it passes the smell test"). Draws full flexion, the pose
% immediately before each unique elimination event, and full extension,
% with the clearance geometry and the active route in one frame.
plt = plotStyleR();
c = plt.hexclr;

transitionIdx = find(any( ...
    routeInfo.active(:,1:end-1) & ~routeInfo.active(:,2:end), 1));
plotIdx = unique([1, transitionIdx, numel(phiD)], 'stable');
nPoseTiles = numel(plotIdx);

% Pose tiles flow four per row; the legend gets its own fifth column
% spanning all tile rows instead of consuming a pose tile.
nPosesPerRow = 4;
nTileRows = ceil(nPoseTiles/nPosesPerRow);
nTileCols = nPosesPerRow + 1;

if nPoseTiles > nTileRows*nPosesPerRow
    error('Too many route poses for the requested tiled layout.')
end

thPlot = linspace(0,2*pi,200).';

figure( ...
    'Name','Optimized 9-point extensor route geometry', ...
    'Color','w', ...
    'Position', [40, 40, 1900, 250+560*nTileRows])

tGeo = tiledlayout( ...
    nTileRows, nTileCols, ...
    'TileSpacing','compact', ...
    'Padding','compact');
hGeoLegend = gobjects(9,1);

for qPlot = 1:nPoseTiles

    ii = plotIdx(qPlot);

    % Native route; p6:p9 are native t1 coordinates -- convert
    % t1 -> ICR -> femur so everything is drawn in one frame.
    Praw = routeInfo.raw(:,:,ii);
    P = Praw;
    for j = 6:9
        qICR = RowVecTrans(T_ICR_t1(:,:,ii), P(j,:));
        P(j,:) = RowVecTrans(T_Pam(:,:,ii), qICR);
    end

    ax = nexttile(tGeo,qPlot);
    hold(ax, 'on')
    colororder(ax, plt.rgbclr)
    hGeo = gobjects(9,1);

    % Femur cylinder clearance
    C = ctx.geo.femurCylCenter;
    R = ctx.geo.femurCylClearRadius;
    hGeo(1) = plot(ax, C(1)+R*cos(thPlot), C(2)+R*sin(thPlot), ...
        '-', 'Color', c{1}, 'LineWidth', plt.lineW);

    % Femur straight-wall clearance
    hGeo(2) = plot(ax, [ctx.geo.femurLineX ctx.geo.femurLineX], ...
        ctx.geo.femurLineY, '-', 'Color', c{2}, 'LineWidth', plt.lineW);

    % True normal-offset condyle clearance
    Q = ctx.geo.femurOffsetBoundary;
    hGeo(3) = plot(ax, [Q(:,1);Q(1,1)], [Q(:,2);Q(1,2)], ...
        '-', 'Color', c{3}, 'LineWidth', plt.lineW);
    if isfield(ctx.geo, 'femurCondyleClipY')
        plot(ax, [ctx.geo.femurCondyleClipX ctx.geo.femurCondyleClipX], ...
            ctx.geo.femurCondyleClipY, '-', 'Color', hGeo(3).Color, ...
            'LineWidth', plt.lineW, 'HandleVisibility', 'off')
    end

    % Lower tibia clearance circle -> femur frame
    L = [ ...
        ctx.geo.tibiaLowerCenter(1) + ...
            ctx.geo.tibiaLowerClearRadius*cos(thPlot), ...
        ctx.geo.tibiaLowerCenter(2) + ...
            ctx.geo.tibiaLowerClearRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];
    Lf = zeros(size(L));
    for kk = 1:size(L,1)
        qICR = RowVecTrans(T_ICR_t1(:,:,ii), L(kk,:));
        Lf(kk,:) = RowVecTrans(T_Pam(:,:,ii), qICR);
    end
    hGeo(4) = plot(ax, Lf(:,1), Lf(:,2), '-', ...
        'Color', c{4}, 'LineWidth', plt.lineW);

    % Upper tibia clearance circle -> femur frame
    U = [ ...
        ctx.geo.tibiaUpperCenter(1) + ...
            ctx.geo.tibiaUpperClearRadius*cos(thPlot), ...
        ctx.geo.tibiaUpperCenter(2) + ...
            ctx.geo.tibiaUpperClearRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];
    Uf = zeros(size(U));
    for kk = 1:size(U,1)
        qICR = RowVecTrans(T_ICR_t1(:,:,ii), U(kk,:));
        Uf(kk,:) = RowVecTrans(T_Pam(:,:,ii), qICR);
    end
    hGeo(5) = plot(ax, Uf(:,1), Uf(:,2), '-', ...
        'Color', c{5}, 'LineWidth', plt.lineW);

    % Local p2/p8 bend radius lines
    tibiaLowerCenter = [ctx.geo.tibiaLowerCenter, 0];
    tibiaLowerCenter = RowVecTrans(T_ICR_t1(:,:,ii), tibiaLowerCenter);
    tibiaLowerCenter = RowVecTrans(T_Pam(:,:,ii), tibiaLowerCenter);

    radiusLineX = NaN;
    radiusLineY = NaN;

    if routeInfo.active(2,ii)
        radiusLineX = [radiusLineX, ...
            ctx.geo.femurCylCenter(1), P(2,1), NaN];
        radiusLineY = [radiusLineY, ...
            ctx.geo.femurCylCenter(2), P(2,2), NaN];
    end
    if routeInfo.active(8,ii)
        radiusLineX = [radiusLineX, tibiaLowerCenter(1), P(8,1)];
        radiusLineY = [radiusLineY, tibiaLowerCenter(2), P(8,2)];
    end
    hGeo(6) = plot(ax, radiusLineX, radiusLineY, '--', ...
        'Color', c{6}, 'LineWidth', plt.lineW);

    % Optimized route
    hGeo(7) = plot(ax, P(:,1), P(:,2), 'o-', ...
        'Color', c{5}, 'LineWidth', plt.lineW, ...
        'MarkerSize', plt.markersz, ...
        'MarkerFaceColor', c{5}, 'MarkerEdgeColor', 'none');

    % Label active route points
    for j = 1:9
        if routeInfo.active(j,ii)
            text(ax, P(j,1), P(j,2), sprintf(' p%d',j), ...
                'FontName', plt.fontN, 'FontSize', plt.lgdFontsz, ...
                'FontWeight', 'bold', 'Interpreter','none')
        end
    end

    % Highlight optimized design endpoints
    hGeo(8) = scatter(ax, P(1,1), P(1,2), plt.scattersz, ...
        'Marker', 's', 'MarkerFaceColor', c{1}, ...
        'MarkerEdgeColor', 'none');
    hGeo(9) = scatter(ax, P(9,1), P(9,2), plt.scattersz, ...
        'Marker', 'd', 'MarkerFaceColor', c{2}, ...
        'MarkerEdgeColor', 'none');

    if qPlot == 1
        hGeoLegend = hGeo;
    end

    axis(ax, 'equal')
    xlabel(ax, 'Femur-frame x, m')
    ylabel(ax, 'Femur-frame y, m')

    if qPlot == 1 || qPlot == nPoseTiles
        tileTitle = sprintf('\\theta_k = %.1f^\\circ', phiD(ii));
    else
        removedNext = find(routeInfo.active(:,ii) & ~routeInfo.active(:,ii+1));
        removedText = strjoin(cellstr(compose('-p%d', removedNext)), ', ');
        tileTitle = sprintf('%s, \\theta_k = %.1f^\\circ', removedText, phiD(ii));
    end
    title(ax, tileTitle, 'Interpreter', 'tex')
    styleAxisR(ax, plt)
end

% Legend occupies the fifth column, spanning every tile row.
axLeg = nexttile(tGeo, nPosesPerRow+1, [nTileRows, 1]);
makeRouteLegendR(axLeg, hGeoLegend, { ...
    'Femur cylinder clr', ...
    'Femur line clr', ...
    'Corrected condyle clr', ...
    'Tibia lower clr', ...
    'Tibia upper clr', ...
    'Local bend radii', ...
    'p1:p9 optimized route', ...
    'Optimized p1', ...
    'Optimized pEnd'}, plt);

%% Static muscle/bone plot at exactly zero knee angle
% One common input list prevents static/animated frame conventions diverging.
% Empty limits fit the bones and routes instead of clipping optimized p1.
bonePlotArgs = {T,T_ICR_t1,phi,pos,p1,pEnd, ...
    {Vas_Int,Vas_Lat,Vas_Med}, ...
    'DisplayRotation',eye(3), ...
    'DisplayAxisMap',eye(3), ...
    'Location',Location,'CrossPoint',CrossPoint,'T_Pam',T_Pam, ...
    'HumanLabels',{'Vastus Intermedius','Vastus Lateralis','Vastus Medialis'}, ...
    'CameraOrbitDeg', 0, ...
    'XLim',[-0.8 0.8], ...
    'YLim',[-1.3 1.0], ...
    'ZLim',[-0.8 0.8]};
run('MuscleBonePlotting.m')

%% Moving muscle/bone plot: zero -> full flexion -> full extension -> zero
% Every frame uses its own Location(:,:,i), including moving/repeated points.
frames = [pos:-1:1, 2:positions, positions-1:-1:pos];
AnimateKneeBoneMuscle(bonePlotArgs{:}, ...
    'FullSkeleton',true, ... % false = femur and tibia only
    'FrameIndices',frames,'PauseTime',0.02,'Loop',false, ...
    'ExportGif',true, ...
    'CameraOrbitDeg', -90, ...
    'XLim',[-0.8 0.8], ...
    'YLim',[-1.3 1.0], ...
    'ZLim',[-0.8 0.8], ...
    'GifFile','Knee_Extensor_20mm.gif', ...
    'FrameRate',20);
clear bonePlotArgs

%% Local functions
function formatAxes(ax,fontName,fontSize,tickLength,xLimits)
set(ax,'FontName',fontName,'FontSize',fontSize,'FontWeight','bold','LineWidth',2,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',tickLength,'XLim',xLimits)
grid(ax,'off')
end

function formatLegend(ax,fontName,fontSize)
lg = legend(ax,'Location','best');
set(lg,'FontName',fontName,'FontSize',fontSize,'FontWeight','bold','Box','off')
end

% Route-plot helpers (port of Opt_run_Ext's plotStyle/loadColors/
% styleAxis/styleLegend/makeRouteLegend, R-suffixed to avoid collisions).
function plt = plotStyleR()
[plt.hexclr, plt.rgbclr] = loadColorsR();
plt.lineW = 2;
plt.scattersz = 60;
plt.markersz = 6;
plt.fontN = 'Arial';
plt.axFontsz = 12;
plt.rulerFontsz = 10;
plt.lgdFontsz = 8;
plt.tickL = [0.025, 0.05];
end

function [hexclr, rgbclr] = loadColorsR()
Colors
hexclr = { ...
    c{1}; ... % gold
    c{2}; ... % orange
    c{3}; ... % light orange
    c{4}; ... % pink
    c{5}; ... % magenta
    c{6}; ... % purple (magenta 2 in Colors.m)
    c{7}};    % indigo
rgbclr = d;   % matching RGB rows, in the same color order
end

function styleAxisR(ax, plt)
grid(ax, 'off')
box(ax, 'off')
set(ax, ...
    'FontName', plt.fontN, ...
    'FontSize', plt.rulerFontsz, ...
    'FontWeight', 'bold', ...
    'LineWidth', plt.lineW, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', plt.tickL)
for k = 1:numel(ax.XAxis)
    ax.XAxis(k).LineWidth = plt.lineW;
    ax.XAxis(k).FontSize = plt.rulerFontsz;
    set(ax.XAxis(k).Label, 'FontName', plt.fontN, ...
        'FontSize', plt.axFontsz, 'FontWeight', 'bold')
end
for k = 1:numel(ax.YAxis)
    ax.YAxis(k).LineWidth = plt.lineW;
    ax.YAxis(k).FontSize = plt.rulerFontsz;
    set(ax.YAxis(k).Label, 'FontName', plt.fontN, ...
        'FontSize', plt.axFontsz, 'FontWeight', 'bold')
end
set(ax.Title, 'FontName', plt.fontN, ...
    'FontSize', plt.axFontsz, 'FontWeight', 'bold')
end

function styleLegendR(lgd, plt)
if isempty(lgd) || ~isgraphics(lgd)
    return
end
set(lgd, 'FontName', plt.fontN, 'FontSize', plt.lgdFontsz, ...
    'FontWeight', 'bold', 'Box', 'off')
end

function makeRouteLegendR(axLeg, hTemplate, labels, plt)
axis(axLeg, 'off')
hold(axLeg, 'on')
hLegend = gobjects(numel(labels),1);
for k = 1:numel(labels)
    hLegend(k) = plot(axLeg, NaN, NaN, 'MarkerEdgeColor', 'none');
    if k <= numel(hTemplate) && isgraphics(hTemplate(k))
        if isprop(hTemplate(k), 'Color')
            hLegend(k).Color = hTemplate(k).Color;
        elseif isprop(hTemplate(k), 'CData')
            hLegend(k).Color = hTemplate(k).CData(1,:);
            hLegend(k).LineStyle = 'none'; % scatter symbols have no line
        end
        if isprop(hTemplate(k), 'LineStyle')
            hLegend(k).LineStyle = hTemplate(k).LineStyle;
        end
        if isprop(hTemplate(k), 'LineWidth')
            hLegend(k).LineWidth = hTemplate(k).LineWidth;
        end
        if isprop(hTemplate(k), 'Marker')
            hLegend(k).Marker = hTemplate(k).Marker;
        end
        if isprop(hTemplate(k), 'MarkerSize')
            hLegend(k).MarkerSize = hTemplate(k).MarkerSize;
        elseif isprop(hTemplate(k), 'SizeData')
            hLegend(k).MarkerSize = sqrt(hTemplate(k).SizeData(1));
        end
        if isprop(hTemplate(k), 'MarkerFaceColor')
            faceColor = hTemplate(k).MarkerFaceColor;
            if isequal(faceColor, 'flat') || isequal(faceColor, 'auto')
                faceColor = hLegend(k).Color;
            end
            hLegend(k).MarkerFaceColor = faceColor;
        end
        hLegend(k).MarkerEdgeColor = 'none';
    end
end
lgd = legend(axLeg, hLegend, labels, ...
    'Location', 'northwest', 'NumColumns', 1, 'AutoUpdate', 'off');
styleLegendR(lgd, plt)
end
