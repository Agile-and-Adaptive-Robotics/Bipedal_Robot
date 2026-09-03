%% Knee extensor data and optimized 20 mm BPA comparison
clc
clear
close all

%% Load the Opt_run_Ext result
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
% Load ctx saved WITH xBest; rebuilding it here could change your bounds,
% model parameters or original reference. See README.md for the full map.
% Only these two variables are required from Opt_run_Ext. Do not load old
% Location/bendMeasure arrays or expect separate saved endpoint variables.
S = load('Vas_Pam_20mm_Result.mat', 'ctx', 'xBest');
if ~isfield(S,'ctx') || ~isfield(S,'xBest')
    error('Knee_Extensor_20mm:MissingResult', ...
        'Vas_Pam_20mm_Result.mat must contain ctx and xBest from the same run.')
end

ctx = S.ctx;
xBest = reshape(S.xBest,1,[]);
validateattributes(xBest,{'numeric'},{'real','finite','numel',8});
p1 = xBest(1:3);
pEnd = xBest(4:6);
rest = xBest(7);
tendon = xBest(8);
KMAX = ctx.KMAX;
kmax = rest*(1-KMAX);

% Rebuild the same route as predictKneeExt20mm, for every knee position.
% The first/last NATIVE rows are p1 and pEnd from xBest. The builder supplies
% p2:p8 and repeated rows for eliminated points, then converts tibia-side
% rows to ICR coordinates. Do not overwrite Location(end,:,:) with raw pEnd.
[Location,bendMeasure,routeInfo] = ...
    buildDistalRingLocation20mm(p1,pEnd,tendon,ctx);

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
xlabel(ax,'\theta_k, °','Interpreter','tex','FontWeight','bold')
ylabel(ax,'Torque Margin, %','FontWeight','bold')
title(ax,'BPA Torque Margin Relative to Human','FontName',fontName,'FontSize',titleFontSize,'FontWeight','bold')
formatLegend(ax,fontName,legendFontSize)

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
