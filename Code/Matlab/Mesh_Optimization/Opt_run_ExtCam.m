%% Opt_run_ExtCam: clip-on tibia cam design for the 20 mm extensor
%
% Ben, 2026-10-12: the 25.4 mm p6/p7 push is a first cut; let the
% solver find the cam shape (p6/p7/p8 + tendon) that maximizes moment
% arm and torque where the human target needs it, subject to the
% no-stretch rule -- at 0 pressure the BPA must sit at or below rest
% length over the whole knee RoM, so it cannot produce torque
% unpressurized.
%
% Design vector (7): xCam = [p6x p6y p7x p7y p8x p8y tendon], p6:p8 in
% the t1 frame. FIXED at the loaded Opt_run_Ext winner: p1, pEnd, rest.
%
% Objective: objective_KneeExtCam20mm (family objective + worst-angle
% margin reward + arm-step penalty + push-neighborhood regularizer).
% Constraints: nonlconExtCam20mm (family set incl. the hard no-stretch
% cRestLength + theta(p6) rule + y-ordering + no points inside the
% stock tibia clearance circles + the VTUG valve floor: the design must
% not exceed the human target at vtugFloorKPa, because below that
% pressure the manifold cannot throttle to the target).
%
% Env knobs (read BEFORE parpool): OPTCAM_SURR_N (default 3000) and
% OPTCAM_PS_N (default 15000) cap the surrogate/patternsearch budgets
% -- shrink both for a fast verification pass, e.g. 24/60.
clear
clc
close all
rehash

surrBudget = 1000;
psBudget = 1500;
if ~isempty(getenv('OPTCAM_SURR_N'))
    surrBudget = str2double(getenv('OPTCAM_SURR_N'));
end
if ~isempty(getenv('OPTCAM_PS_N'))
    psBudget = str2double(getenv('OPTCAM_PS_N'));
end

% Festo VTUG valve-manifold floor (Ben, 2026-10-12): the manifold does
% not like to open/close below ~200 kPa, so the winner must still NEED
% at least this pressure at every angle that carries a target. This,
% not the search boxes, is what caps the cam size.
vtugFloorKPa = 200;

% Marks a full optimizer pass; the (commented) result save clears it
% just before saving, mirroring Opt_run_Ext.
liveRun = true;

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
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
% buildKneeExtContext20mm reads OpenSim_Vasti_Results.txt and its
% Xi-pick mat by bare name; both live in Testing_Data\2022_02_Festo.
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

%% Context and fixed design
ctx = buildKneeExtContext20mm();

resDir = fullfile(scriptDir, 'Results');
resultFile = fullfile(resDir, 'Vas_Pam_20mm_Result.mat');
S = load(resultFile, 'xBest', 'XiUsed', 'ctx');
if ~isfield(S,'xBest')
    error('Opt_run_ExtCam:MissingResult', ...
        '%s must contain xBest from Opt_run_Ext.', resultFile)
end
if isfield(S,'XiUsed')
    xiRun = S.XiUsed;
elseif isfield(S,'ctx') && isfield(S.ctx,'Xi0')
    xiRun = [S.ctx.Xi0, S.ctx.Xi1, S.ctx.Xi2, S.ctx.Xi3];
else
    error('Opt_run_ExtCam:MissingResult', ...
        '%s must contain xBest and XiUsed (or a saved ctx).', resultFile)
end

xiBuilt = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
if ~isequal(xiBuilt, xiRun)
    error('Opt_run_ExtCam:XiMismatch', ...
        'Builder Xi [%.6g %.6g %.6g %.6g] does not equal run Xi [%.6g %.6g %.6g %.6g].', ...
        xiBuilt, xiRun)
end

xRef8 = reshape(S.xBest,1,[]);
xFixed = xRef8(1:7);      % [p1, pEnd, rest] stay fixed
tendonRef = xRef8(8);

%% Baseline constructed contacts + Ben's first-cut 25.4 mm push
[~,~,routeInfoBase] = buildDistalRingLocation20mm( ...
    xFixed(1:3), xFixed(4:6), tendonRef, ctx);
seedBase = routeInfoBase.seed20;
thetaP6BaseD = atan2d(seedBase(6,2), seedBase(6,1));

adjPush = 0.0254;   % m, radial push of p6/p7 about the t1 origin
refPts = seedBase(6:8,1:2);
refPts(1,:) = refPts(1,:) + adjPush*refPts(1,:)/norm(refPts(1,:));
refPts(2,:) = refPts(2,:) + adjPush*refPts(2,:)/norm(refPts(2,:));

xCamBase = [reshape(seedBase(6:8,1:2),1,6), tendonRef];
xCamPush = [reshape(refPts,1,6), tendonRef];

%% Bounds (t1 frame; the part grows on the proximal-anterior tibia)
lbCam = [ ...
    seedBase(6,1)-0.010, seedBase(6,2)-0.020, ...
    seedBase(7,1)-0.010, seedBase(7,2)-0.040, ...
    seedBase(8,1)-0.010, seedBase(8,2)-0.040, ...
    ctx.geo.tendonMin];
ubCam = [ ...
    seedBase(6,1)+0.060, seedBase(6,2)+0.040, ...
    seedBase(7,1)+0.060, seedBase(7,2)+0.040, ...
    seedBase(8,1)+0.050, seedBase(8,2)+0.030, ...
    0.3];

obj = @(x) objective_KneeExtCam20mm(x, ctx, xFixed, refPts);
nonlcon = @(x) nonlconExtCam20mm(x, ctx, xFixed, thetaP6BaseD, vtugFloorKPa);
objconstr = @(x) objconstrExt20mm(x, obj, nonlcon);

if isempty(gcp('nocreate'))
    parpool;
end

%% Surrogate global search
optsG = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', surrBudget, ...
    'MinSampleDistance', 0.001, ...
    'ConstraintTolerance', 1e-6, ...
    'InitialPoints', [xCamBase; xCamPush]);

[xG, fG] = surrogateopt(objconstr, lbCam, ubCam, optsG);

%% Pattern-search refinement
optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', psBudget, ...
    'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-6, ...
    'PlotFcn', {@psplotbestf, ...
                @psplotfuncount, ...
                @psplotmeshsize, ...
                @psplotmaxconstr}, ...
    'OutputFcn', @patternProgress);

% Do not throw away a better known design before patternsearch.
fPush = obj(xCamPush);
if fPush <= fG
    xStart = xCamPush;
else
    xStart = xG;
end

[xBestCam, fBestCam] = patternsearch( ...
    obj, xStart, [], [], [], [], lbCam, ubCam, nonlcon, optsP);

%% Winner evaluation and report
[xFullBest, ctxCamBest] = camToExtDesign20mm(xBestCam, ctx, xFixed);
predBest = predictKneeExt20mm(xFullBest, ctxCamBest);
[cBest, ~] = nonlconExtCam20mm(xBestCam, ctx, xFixed, thetaP6BaseD, vtugFloorKPa);

humanBest = interp1(ctx.humanAngleD, ctx.humanTorque, ctx.phiD, 'pchip', 'extrap');
marginBest = (predBest.TorqueZ(:) - humanBest(:))./humanBest(:);
requiredBest = (1 + ctx.requiredTorqueMargin)*abs(humanBest(:));
marginReqBest = (abs(predBest.TorqueZ(:)) - requiredBest(:))./ctx.torqueScale;
maxArmStepBest = max(abs(diff(predBest.momentArm)));

fprintf('\n========== CAM-OPTIMIZED DESIGN ==========\n')
fprintf('objective = %.6g (surrogate best %.6g, push first cut %.6g)\n', ...
    fBestCam, fG, fPush)
fprintf('p6 = [%.6f %.6f] m, theta %.3f deg (baseline limit %.3f deg)\n', ...
    xBestCam(1:2), atan2d(xBestCam(2),xBestCam(1)), thetaP6BaseD)
fprintf('p7 = [%.6f %.6f] m\n', xBestCam(3:4))
fprintf('p8 = [%.6f %.6f] m\n', xBestCam(5:6))
fprintf('tendon = %.6f m (baseline %.6f m)\n', xBestCam(7), tendonRef)
fprintf('max constraint = %+.6f (cRestLength %+.6f, cContraction %+.6f, ', ...
    max(cBest), cBest(5), cBest(6))
fprintf('cTendonMax %+.6f, cTheta %+.6f, cVtug@%dkPa %+.6f)\n', ...
    cBest(1), cBest(7), vtugFloorKPa, cBest(end))
fprintf('min torque margin vs required = %+.2f%%; vs plain target = %+.2f%%\n', ...
    100*min(marginReqBest), 100*min(marginBest))
fprintf('moment arm: min %.1f mm, mean %.1f mm, max %.1f mm; max pose step %.2f mm\n', ...
    1000*min(predBest.momentArm), 1000*mean(predBest.momentArm), ...
    1000*max(predBest.momentArm), 1000*maxArmStepBest)

%% No-stretch / 0-pressure check: baseline vs first cut vs winner
% The rule (Ben): at 0 pressure the BPA produces no torque, i.e. it is
% never stretched over the knee RoM. cRestLength <= 0 is the hard
% model form; maxTq0 is the direct 0-kPa check through the class.
[xFullPush, ctxCamPush] = camToExtDesign20mm(xCamPush, ctx, xFixed);
designs = { ...
    'baseline xBest ', xRef8, ctx; ...
    'push first cut', xFullPush, ctxCamPush; ...
    'cam winner    ', xFullBest, ctxCamBest};

fprintf('\nNo-stretch / 0-pressure check (cRestLength <= 0 = feasible):\n')
fprintf('%-16s %12s %12s %14s %16s\n', 'design', 'cRestLength', 'min strain', ...
    'max stretch mm', 'max |Tq| @0kPa')
for iD = 1:size(designs,1)
    xFullD = designs{iD,2};
    ctxD = designs{iD,3};
    predD = predictKneeExt20mm(xFullD, ctxD);
    if predD.ok
        cRestD = predD.cRestLength;
        minStrainD = min(predD.strain_p);
        stretchD = 1000*max(0, -minStrainD*predD.rest);
        tq0D = zeroPressureTorque(xFullD, ctxD);
    else
        cRestD = NaN; minStrainD = NaN; stretchD = NaN; tq0D = NaN;
    end
    fprintf('%-16s %+12.6f %+12.6f %14.3f %16.4f\n', ...
        designs{iD,1}, cRestD, minStrainD, stretchD, tq0D)
end

fprintf('\nPaste into the Knee_Extensor_20mm adj block:\n')
fprintf('tendonAdj = %.6f;%% m\n', xBestCam(7))
fprintf('p6Adj = [%.6f, %.6f];%% t1 frame\n', xBestCam(1:2))
fprintf('p7Adj = [%.6f, %.6f];%% t1 frame\n', xBestCam(3:4))
fprintf('p8Adj = [%.6f, %.6f];%% t1 frame\n', xBestCam(5:6))

%% Winner required-pressure envelope (VTUG floor check)
% Coarse 20 kPa scan over the angles that carry an extension target:
% the winner is valve-usable only if min required pressure stays at or
% above the VTUG floor (cVtug enforces this in the solver).
[pminW, pmaxW] = requiredPressureScanCam(xFullBest, ctxCamBest);
fprintf('\nWinner required pressure: min %.0f kPa, max %.0f kPa (VTUG floor %.0f kPa)\n', ...
    pminW, pmaxW, vtugFloorKPa)
if pminW < vtugFloorKPa - 10
    fprintf(['WARNING: min required pressure is below the VTUG floor ', ...
        '-- check the cVtug constraint.\n'])
end

% Dated FULL-WORKSPACE result capture into Results (Ben directive,
% 2026-09-20 pattern). liveRun is set only by a full Opt_run_ExtCam
% pass and cleared before saving, so rerunning display sections from a
% loaded result mat cannot mint a new dated mat. Uncomment for a real
% run.
% if exist('liveRun', 'var')
%     stamp = char(string(datetime('now'),'yyyyMMdd_HHmm'));
%     resDir = fullfile(scriptDir, 'Results');
%     resultFile = fullfile(resDir, sprintf('Vas_Pam_20mm_Result_cam_%s.mat', stamp));
%     XiUsed = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
%     clear liveRun
%     save(resultFile)
%     fprintf('Saved %s\n', resultFile)
% end

%% Winner diagnostic figures
xLim = [min(ctx.phiD), max(ctx.phiD)];

hFigW1 = figure('Name','Cam torque vs human target','Color','w')
hold on
plot(ctx.phiD, predBest.TorqueZ, 'LineWidth', 2)
plot(ctx.humanAngleD, ctx.humanTorque, 'LineWidth', 2)
plot(ctx.phiD, requiredBest, '--', 'LineWidth', 2)
xlim(xLim); grid off
xlabel('Knee angle, deg'); ylabel('Torque, N m')
legend('Cam BPA prediction', char(ctx.targetName), 'Required (+5%)', ...
    'Location', 'best')
title('Cam-optimized extensor torque')

hFigW2 = figure('Name','Cam torque margin','Color','w')
hold on
plot(ctx.phiD, 100*marginBest, 'LineWidth', 2)
yline(0, '--'); yline(100*ctx.requiredTorqueMargin, ':')
xlim(xLim); grid off
xlabel('Knee angle, deg'); ylabel('Torque margin, %')
legend('vs plain target', 'zero', 'required +5%', 'Location', 'best')
title('Cam torque margin relative to human')

hFigW3 = figure('Name','Cam BPA strain','Color','w')
hold on
plot(ctx.phiD, predBest.strain_f, 'LineWidth', 2)
plot(ctx.phiD, predBest.strain_p, '--', 'LineWidth', 2)
yline(0, ':'); yline(ctx.KMAX, ':')
xlim(xLim); grid off
xlabel('Knee angle, deg'); ylabel('Strain')
legend('strain_f (incl Xi3)', 'strain_p', '0 (no stretch)', 'KMAX', ...
    'Location', 'best')
title('Cam BPA strain over the RoM')

hFigW4 = figure('Name','Cam moment arm','Color','w')
plot(ctx.phiD, 1000*predBest.momentArm, 'LineWidth', 2)
xlim(xLim); grid off
xlabel('Knee angle, deg'); ylabel('Moment arm, mm')
title(sprintf('Cam moment arm (max pose step %.2f mm)', 1000*maxArmStepBest))

% The figures section is silent in the console (and invisible under
% matlab -batch), so say what was created and export PNGs to Results.
winnerFigs = [hFigW1, hFigW2, hFigW3, hFigW4];
for kF = 1:numel(winnerFigs)
    exportgraphics(winnerFigs(kF), fullfile(resDir, ...
        sprintf('Vas_Pam_20mm_cam_winner_%d.png', kF)), 'Resolution', 150);
end
fprintf('\nCreated %d winner diagnostic figures (exported to %s as Vas_Pam_20mm_cam_winner_*.png).\n', ...
    numel(winnerFigs), resDir)

%% Local functions
function [pmin, pmax] = requiredPressureScanCam(xFull, ctxEval)
% Coarse required-pressure scan (20 kPa grid; torque is monotone in
% pressure) over the angles that carry an extension target.

Pgrid = (0:20:620)';
humanTgt = interp1(ctxEval.humanAngleD, ctxEval.humanTorqueAbs, ...
    ctxEval.phiD, 'pchip', 'extrap');
hasTgt = isfinite(humanTgt) & humanTgt > 0;

restV = xFull(7);
kmaxV = restV*(1 - ctxEval.KMAX);
[LocationV, bendV] = buildDistalRingLocation20mm( ...
    xFull(1:3), xFull(4:6), xFull(8), ctxEval);

Tq = nan(numel(Pgrid), numel(ctxEval.phiD));
for k = 1:numel(Pgrid)
    b = MonoPamDataExplicit_balanceX3(ctxEval.Name, LocationV, ...
        ctxEval.CrossPoint, ctxEval.Dia, ctxEval.T_Pam, restV, kmaxV, ...
        xFull(8), ctxEval.fitting, Pgrid(k), ctxEval.Xi0, ctxEval.Xi1, ...
        ctxEval.Xi2, ctxEval.Xi3, ctxEval.wraps, ctxEval.phiD, ...
        ctxEval.BPAcount, bendV);
    Tq(k,:) = abs(b.Torque_p(:,3)).';
end

Preq = nan(numel(ctxEval.phiD),1);
% find() returns a column; for-loops iterate over columns, so transpose
% to iterate element-wise (see buildDistalRingLocation20mm's find(...).').
for i = find(hasTgt).'
    col = Tq(:,i);
    fin = isfinite(col);
    if ~any(fin) || max(col(fin)) < humanTgt(i)
        continue
    end
    k = find(col >= humanTgt(i), 1);
    if k <= 1
        Preq(i) = 0;
    else
        pq = Pgrid([k-1 k]);
        [xs, isrt] = sort(col([k-1 k]));
        Preq(i) = interp1(xs, pq(isrt), humanTgt(i));
    end
end

pmin = min(Preq(find(hasTgt)),[],'omitnan');
pmax = max(Preq(find(hasTgt)),[],'omitnan');

end

function tq0 = zeroPressureTorque(xFull, ctxEval)
% Max |extension torque| of the design at 0 kPa. Zero when the BPA is
% never stretched (cRestLength <= 0 guarantees strain_p >= 0).

rest = xFull(7);
kmax = rest*(1-ctxEval.KMAX);
[Location,bendMeasure] = buildDistalRingLocation20mm( ...
    xFull(1:3), xFull(4:6), xFull(8), ctxEval);
bpa0 = MonoPamDataExplicit_balanceX3(ctxEval.Name, Location, ...
    ctxEval.CrossPoint, ctxEval.Dia, ctxEval.T_Pam, rest, kmax, ...
    xFull(8), ctxEval.fitting, 0, ctxEval.Xi0, ctxEval.Xi1, ...
    ctxEval.Xi2, ctxEval.Xi3, ctxEval.wraps, ctxEval.phiD, ...
    ctxEval.BPAcount, bendMeasure);
tq0 = max(abs(bpa0.Torque_p(:,3)));

end

function [stop, options, optchanged] = patternProgress( ...
        optimValues, options, flag)
% Concise pattern-search progress (verbatim from Opt_run_Ext).

stop = false;
optchanged = false;

if strcmp(flag, 'init')
    fprintf(['\nPattern search progress:\n' ...
        ' Iteration    Function evaluations       Best f' ...
        '       Mesh size     Max constraint\n'])
    return
end

if ~strcmp(flag, 'iter') && ~strcmp(flag, 'done')
    return
end

maxConstraint = 0;
if isfield(optimValues, 'nonlinineq') && ...
        ~isempty(optimValues.nonlinineq)
    maxConstraint = max([0; optimValues.nonlinineq(:)]);
end

fprintf('%10d %23d %12.6g %15.6g %18.6g\n', ...
    optimValues.iteration, ...
    optimValues.funccount, ...
    optimValues.fval, ...
    optimValues.meshsize, ...
    maxConstraint)

if strcmp(flag, 'done')
    fprintf('Pattern search finished.\n')
end

end
