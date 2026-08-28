%% Knee Extensor 20 mm - reproduce optimizer result
% Loads the saved optimizer solution and rebuilds the 20 mm extensor BPA
% object without rerunning surrogateopt or patternsearch.

clc;
clear;
close all;

%% Project paths
scriptDir = fileparts(mfilename('fullpath'));
meshDir = fileparts(scriptDir);
matlabDir = fileparts(meshDir);
codeDir = fileparts(matlabDir);
repoDir = fileparts(codeDir);

addpath(meshDir)
addpath(matlabDir)
addpath(fullfile(matlabDir, 'Functions'))
addpath(fullfile(matlabDir, 'Functions', 'ModernRobotics'))
addpath(fullfile(matlabDir, 'Robot_Data'))

resultFile = fullfile(scriptDir, 'Vas_Pam_20mm_Result.mat');
testDataDir = fullfile(repoDir, 'Testing_Data', '2022_02_Festo');

if ~isfile(resultFile)
    error('Knee_Extensor_20mm:MissingResult', ...
        'Cannot find optimizer result file: %s', resultFile)
end

%% Load optimized design vector
S = load(resultFile, 'xBest', 'fBest');

if ~isfield(S, 'xBest')
    error('Knee_Extensor_20mm:MissingXBest', ...
        'The result file does not contain xBest.')
end

xBest = S.xBest(:).';

if isfield(S, 'fBest')
    fBest = S.fBest;
else
    fBest = NaN;
end

%% Rebuild the optimizer context and BPA class object
% buildKneeExtContext20mm reads OpenSim_Vasti_Results.txt from the current
% folder, so temporarily run the rebuild from the testing-data folder.
startDir = pwd;
restoreStartDir = onCleanup(@() cd(startDir));

if isfolder(testDataDir)
    cd(testDataDir)
else
    warning('Knee_Extensor_20mm:MissingTestDataDir', ...
        'Expected testing-data folder not found: %s', testDataDir)
end

ctx = buildKneeExtContext20mm();
predBest = predictKneeExt20mm(xBest, ctx);

if ~predBest.ok
    error('Knee_Extensor_20mm:PredictionFailed', ...
        'Prediction failed: %s', predBest.failReason)
end

%% Class-script compatibility variables
% These names mirror the older Knee_Extensor_10mm / Knee_Extensor_40mm
% scripts, while keeping the optimizer route and X3 bend model intact.
Vas_Pam_20mm = predBest.bpa;
Vas_Pam = Vas_Pam_20mm;

Location = predBest.Location;
Praw = predBest.routeInfo.raw;
routeInfo = predBest.routeInfo;

CrossPoint = ctx.CrossPoint;
Dia = ctx.Dia;
T_Pam = ctx.T_Pam;
T_Pam_inv = ctx.T_Pam_inv;
T_t1_ICR = ctx.T_t1_ICR;
T_ICR_t1 = ctx.T_ICR_t1;

rest = predBest.rest;
kmax = predBest.kmax;
tendon = predBest.tendon;
fitting = ctx.fitting;
pres = ctx.targetPressure;

phi = ctx.phi;
phiD = ctx.phiD;

ForceR = predBest.bpa.F_p;
TorqueR = predBest.bpa.Torque_p;
TorqueZ = predBest.TorqueZ;

humanTargetZ = interp1( ...
    ctx.humanAngleD, ...
    ctx.humanTorque, ...
    phiD, ...
    'pchip', ...
    'extrap');

torqueMargin = (TorqueZ - humanTargetZ)./humanTargetZ;

%% Summary
fprintf('\n========== 20 mm OPTIMIZER REPRODUCTION ==========\n')
fprintf('Result file: %s\n', resultFile)
fprintf('fBest = %.6g\n', fBest)
fprintf('p1   = [% .3f % .3f % .3f] m\n', predBest.p1)
fprintf('pEnd = [% .3f % .3f % .3f] m\n', predBest.pEnd)
fprintf('rest = %.3f m, kmax = %.3f m, tendon = %.3f m\n', ...
    rest, kmax, tendon)
fprintf('TorqueZ range = %.3f to %.3f N m\n', ...
    min(TorqueZ), max(TorqueZ))
fprintf('Max undeformed path length = %.3f m at %.1f deg\n', ...
    predBest.maxPathLength0, phiD(predBest.idxMaxPath0))
fprintf('==================================================\n')

%% Quick reproduction plots
figure('Name', '20 mm optimized extensor torque', 'Color', 'w')
plot(phiD, TorqueZ, 'LineWidth', 2)
hold on
plot(ctx.humanAngleD, ctx.humanTorque, '--', 'LineWidth', 2)
xlabel('Knee angle, deg')
ylabel('Torque, N m')
title('20 mm optimized extensor torque')
legend('20 mm BPA optimizer result', char(ctx.targetName), ...
    'Location', 'best')

figure('Name', '20 mm optimized route state', 'Color', 'w')
imagesc(phiD, 1:ctx.routeRows, double(routeInfo.active))
axis xy
xlabel('Knee angle, deg')
ylabel('Route point')
yticks(1:ctx.routeRows)
yticklabels(compose('p%d', 1:ctx.routeRows))
title('20 mm optimized active route points')
colorbar
