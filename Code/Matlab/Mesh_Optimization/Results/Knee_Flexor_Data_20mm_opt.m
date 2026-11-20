%% Knee flexor: reproduce the reported optimum and probe nearby geometry
% This script uses the same prediction and collision functions as Opt_run.
% The optimized values below were copied from the 1500-evaluation
% surrogateopt + 30000-evaluation patternsearch run. They were printed to
% six decimal places, so very small differences from the live xBest are
% expected at active constraints.

% If this is run immediately after Opt_run, preserve the full-precision
% optimizer result instead of discarding it. Otherwise use the rounded
% command-window values recorded below.
if exist('xBest','var') && isnumeric(xBest) && numel(xBest) == 8
    liveXBest = xBest(:).';
else
    liveXBest = [];
end
clearvars -except liveXBest
clc
close all
rehash

ctx = buildKneeFlexorContext20mm();
geo = buildGeoExclusion();
idxP2 = 4:6;
constraintTolerance = 1e-6;

% [p1_x p1_y p1_z p2_x p2_y p2_z rest tendon]
xReportedRounded = [ ...
   -0.020312,  0.113333,  0.058283, ...
    0.001872, -0.019369,  0.058246, ...
    0.481596,  0.054253];

if isempty(liveXBest)
    xReported = xReportedRounded;
    fprintf('Using the six-decimal command-window design values.\n')
else
    xReported = liveXBest;
    fprintf('Using full-precision xBest preserved from the workspace.\n')
end

reportedObjective = 0.0663064012;
reportedLmExtension = 0.432738;
reportedLmFlexion = 0.383266;

original = evaluateCandidate(ctx.x0, "Original", ctx, geo, idxP2, ...
    constraintTolerance);
baseline = evaluateCandidate(xReported, "Reported optimum", ctx, geo, ...
    idxP2, constraintTolerance);

if ~baseline.ok
    error('Reported design prediction failed: %s', baseline.failReason)
end

%% Baseline reproduction report
[~, idxFullExtension] = max(ctx.phiD);
[~, idxFullFlexion] = min(ctx.phiD);

fprintf('\n========== REPORTED OPTIMUM REPRODUCTION ==========\n')
fprintf('p1 = [%.6f, %.6f, %.6f] m\n', baseline.x(1:3))
fprintf('p2 = [%.6f, %.6f, %.6f] m\n', baseline.x(4:6))
fprintf('rest length   = %.6f m\n', baseline.x(7))
fprintf('tendon length = %.6f m\n', baseline.x(8))
fprintf('kmax          = %.6f m\n', baseline.pred.kmax)
fprintf('objective from rounded values = %.9g\n', baseline.objective)
fprintf('reported full-precision objective = %.9g\n', reportedObjective)
fprintf('objective difference = %+.9g\n', ...
    baseline.objective - reportedObjective)
fprintf('Lm at extension = %.6f m (reported %.6f m)\n', ...
    baseline.pred.activeLength(idxFullExtension), reportedLmExtension)
fprintf('Lm at flexion   = %.6f m (reported %.6f m)\n', ...
    baseline.pred.activeLength(idxFullFlexion), reportedLmFlexion)
fprintf('minimum torque margin = %+.2f%%\n', ...
    100*baseline.minTorqueMarginFraction)
fprintf('mean remaining torque shortfall = %.2f%%\n', ...
    100*baseline.meanTorqueShortfallFraction)
fprintf('rest-series reach c = %+.6f mm\n', ...
    1000*baseline.pred.cRestLength)
fprintf('minimum tibia clearance = %.6f mm\n', ...
    1000*baseline.collisionInfo.minClearanceTibia)
fprintf('minimum femur clearance = %.6f mm\n', ...
    1000*baseline.collisionInfo.minClearanceFemur)
fprintf('collision feasible = %d\n', baseline.collisionFeasible)
fprintf('strain feasible    = %d\n', baseline.strainFeasible)
fprintf('===================================================\n')

%% Publication-style baseline figures: one figure per result
style = publicationStyle();
plotBaselineResults(original, baseline, ctx, style)

%% Equal negative-x shift sweep
% p1_x and p2_x are shifted by the same amount in their own frames.
% Two versions are evaluated:
%   raw       - only p1_x and p2_x change;
%   corrected - tendon is lengthened just enough to preserve the physical
%               rest-series reach check at extension, plus 0.1 mm margin.
shiftMM = (0:-1:-30).';
reachBuffer = 0.0001;
nShift = numel(shiftMM);
rawShift = repmat(emptyCandidate(), nShift, 1);
correctedShift = repmat(emptyCandidate(), nShift, 1);
tendonAddedMM = zeros(nShift, 1);

for i = 1:nShift
    xRaw = xReported;
    xRaw([1 4]) = xRaw([1 4]) + shiftMM(i)/1000;
    rawShift(i) = evaluateCandidate(xRaw, ...
        sprintf('Both x %+.0f mm', shiftMM(i)), ...
        ctx, geo, idxP2, constraintTolerance);

    [xCorrected, tendonAddedMM(i)] = correctSeriesReach( ...
        xRaw, ctx, reachBuffer);
    correctedShift(i) = evaluateCandidate(xCorrected, ...
        sprintf('Both x %+.0f mm, reach corrected', shiftMM(i)), ...
        ctx, geo, idxP2, constraintTolerance);
end

rawTable = candidateTable(shiftMM, rawShift, zeros(nShift,1));
correctedTable = candidateTable(shiftMM, correctedShift, tendonAddedMM);

fprintf('\nEqual -x sweep, raw geometry:\n')
disp(rawTable)
fprintf('\nEqual -x sweep, tendon-corrected reach:\n')
disp(correctedTable)

bestTorqueIdx = bestFeasibleIndex(correctedShift, 'torque');
bestObjectiveIdx = bestFeasibleIndex(correctedShift, 'objective');

if ~isempty(bestTorqueIdx)
    bestShift = correctedShift(bestTorqueIdx);
    fprintf(['\nBest feasible corrected -x shift by worst-case torque ' ...
        'margin: %+.0f mm\n'], shiftMM(bestTorqueIdx))
    printCandidateSummary(bestShift)
else
    warning('No all-feasible corrected -x candidate was found.')
end

if ~isempty(bestObjectiveIdx)
    fprintf('\nBest feasible corrected -x shift by objective: %+.0f mm\n', ...
        shiftMM(bestObjectiveIdx))
    printCandidateSummary(correctedShift(bestObjectiveIdx))
end

plotShiftSweep(rawTable, correctedTable, style)

if ~isempty(bestTorqueIdx) && bestTorqueIdx ~= 1
    plotCandidateTorque(baseline, correctedShift(bestTorqueIdx), ctx, ...
        style, 'Best Equal -x Shift')
end

%% p1 y + 0.012 m test
xP1Y = xReported;
xP1Y(2) = xP1Y(2) + 0.012;
p1YRaw = evaluateCandidate(xP1Y, "p1 y +12 mm", ctx, geo, idxP2, ...
    constraintTolerance);
[xP1YCorrected, p1YTendonAddedMM] = correctSeriesReach( ...
    xP1Y, ctx, reachBuffer);
p1YCorrected = evaluateCandidate(xP1YCorrected, ...
    "p1 y +12 mm, reach corrected", ctx, geo, idxP2, ...
    constraintTolerance);

fprintf('\np1 y +12 mm, raw geometry:\n')
printCandidateSummary(p1YRaw)
fprintf('p1 y +12 mm, tendon-corrected reach (added %.3f mm):\n', ...
    p1YTendonAddedMM)
printCandidateSummary(p1YCorrected)

plotCandidateTorque(baseline, p1YRaw, ctx, style, 'p1 y +12 mm')

%% Save numeric results from this run
save('Knee_Flexor_Design_Sweep_Results.mat', ...
    'xReported', 'baseline', 'original', 'shiftMM', 'rawShift', ...
    'correctedShift', 'rawTable', 'correctedTable', ...
    'p1YRaw', 'p1YCorrected')
writetable(rawTable, 'Knee_Flexor_NegativeX_Sweep_Raw.csv')
writetable(correctedTable, ...
    'Knee_Flexor_NegativeX_Sweep_ReachCorrected.csv')

%% Local functions
function result = emptyCandidate()
result = struct( ...
    'label', "", 'x', nan(1,8), 'ok', false, 'failReason', "", ...
    'pred', struct(), 'objective', Inf, 'cCollision', nan(4,1), ...
    'collisionInfo', struct(), 'collisionFeasible', false, ...
    'reachFeasible', false, 'strainFeasible', false, ...
    'withinBounds', false, 'allFeasible', false, ...
    'torqueMarginFraction', nan(0,1), ...
    'minTorqueMarginFraction', -Inf, ...
    'meanTorqueShortfallFraction', Inf);
end

function result = evaluateCandidate(x, label, ctx, geo, idxP2, tol)
result = emptyCandidate();
result.label = string(label);
result.x = x(:).';
result.withinBounds = all(result.x >= ctx.lb - tol) && ...
    all(result.x <= ctx.ub + tol);

pred = predictKneeFlexor20mm(result.x, ctx);
result.pred = pred;
result.ok = pred.ok;
result.failReason = pred.failReason;
if ~pred.ok
    return
end

[cCollision, ~, collisionInfo] = nonlconExclusion( ...
    result.x, geo, ctx, idxP2);
result.cCollision = cCollision;
result.collisionInfo = collisionInfo;
result.collisionFeasible = all(cCollision <= tol);
result.reachFeasible = pred.cRestLength <= tol;
result.strainFeasible = ...
    max(pred.strain) <= ctx.maxRelStrain*pred.KMAX + tol && ...
    min(pred.strain) >= ctx.minStrain - tol;
result.allFeasible = result.withinBounds && result.collisionFeasible && ...
    result.reachFeasible && result.strainFeasible;

humanAbs = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, ...
    'pchip', 'extrap');
valid = humanAbs > 100*eps;
margin = nan(size(pred.TorqueZ));
margin(valid) = abs(pred.TorqueZ(valid))./humanAbs(valid) - 1;
result.torqueMarginFraction = margin;
result.minTorqueMarginFraction = min(margin(valid));
result.meanTorqueShortfallFraction = mean(max(0, -margin(valid)));
result.objective = objectiveFromPrediction(result.x, pred, ctx);
end

function J = objectiveFromPrediction(x, pred, ctx)
if ~pred.ok || any(~isfinite(pred.TorqueZ)) || ...
        any(~isfinite(pred.strain))
    J = 1e11;
    return
end

humanAbs = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, ...
    'pchip', 'extrap');
robotAbs = abs(pred.TorqueZ);
torqueDeficit = max(0, humanAbs(:) - robotAbs(:));
torqueOvershoot = max(0, robotAbs(:) - humanAbs(:));
Jtorque = mean((torqueDeficit/ctx.torqueScale).^2);
Jovershoot = 1e-3*mean((torqueOvershoot/ctx.torqueScale).^2);
JrestLength = 1e6*max(0, pred.cRestLength).^2;
maxStrainAllowed = ctx.maxRelStrain*pred.KMAX;
JstrainHi = 1e4*max(0, max(pred.strain)-maxStrainAllowed).^2;
JstrainLo = 1e4*max(0, ctx.minStrain-min(pred.strain)).^2;
dx = x(:)-ctx.x0(:);
geomScale = [0.030; 0.300; 0.030; 0.030; 0.030; 0.030];
Jgeom = 1e-2*sum((dx(1:6)./geomScale).^2);
Jlen = 1e-3*((x(7)-ctx.x0(7))/0.040).^2;
J = 100*Jtorque + Jovershoot + JrestLength + ...
    JstrainHi + JstrainLo + Jgeom + Jlen;
end

function [xCorrected, tendonAddedMM] = correctSeriesReach(x, ctx, buffer)
xCorrected = x;
pred = predictKneeFlexor20mm(xCorrected, ctx);
if ~pred.ok
    tendonAddedMM = NaN;
    return
end
needed = max(0, pred.cRestLength + buffer);
xCorrected(8) = min(ctx.ub(8), xCorrected(8) + needed);
tendonAddedMM = 1000*(xCorrected(8)-x(8));
end

function T = candidateTable(shiftMM, result, tendonAddedMM)
n = numel(result);
objective = nan(n,1);
minMarginPct = nan(n,1);
meanShortfallPct = nan(n,1);
minTibiaClearanceMM = nan(n,1);
minFemurClearanceMM = nan(n,1);
reachViolationMM = nan(n,1);
collisionFeasible = false(n,1);
reachFeasible = false(n,1);
strainFeasible = false(n,1);
allFeasible = false(n,1);

for i = 1:n
    objective(i) = result(i).objective;
    minMarginPct(i) = 100*result(i).minTorqueMarginFraction;
    meanShortfallPct(i) = 100*result(i).meanTorqueShortfallFraction;
    collisionFeasible(i) = result(i).collisionFeasible;
    reachFeasible(i) = result(i).reachFeasible;
    strainFeasible(i) = result(i).strainFeasible;
    allFeasible(i) = result(i).allFeasible;
    if result(i).ok
        reachViolationMM(i) = 1000*result(i).pred.cRestLength;
        minTibiaClearanceMM(i) = ...
            1000*result(i).collisionInfo.minClearanceTibia;
        minFemurClearanceMM(i) = ...
            1000*result(i).collisionInfo.minClearanceFemur;
    end
end

T = table(shiftMM, tendonAddedMM, objective, minMarginPct, ...
    meanShortfallPct, minTibiaClearanceMM, minFemurClearanceMM, ...
    reachViolationMM, collisionFeasible, reachFeasible, ...
    strainFeasible, allFeasible, ...
    'VariableNames', {'Shift_mm','TendonAdded_mm','Objective', ...
    'MinTorqueMargin_pct','MeanTorqueShortfall_pct', ...
    'MinTibiaClearance_mm','MinFemurClearance_mm', ...
    'ReachViolation_mm','CollisionFeasible','ReachFeasible', ...
    'StrainFeasible','AllFeasible'});
end

function idx = bestFeasibleIndex(result, criterion)
feasible = [result.allFeasible].';
if ~any(feasible)
    idx = [];
    return
end
switch lower(criterion)
    case 'torque'
        score = [result.minTorqueMarginFraction].';
        score(~feasible) = -Inf;
        [~, idx] = max(score);
    case 'objective'
        score = [result.objective].';
        score(~feasible) = Inf;
        [~, idx] = min(score);
    otherwise
        error('Unknown selection criterion: %s', criterion)
end
end

function printCandidateSummary(result)
if ~result.ok
    fprintf('prediction failed: %s\n', result.failReason)
    return
end
fprintf('p1 = [%.6f, %.6f, %.6f] m\n', result.x(1:3))
fprintf('p2 = [%.6f, %.6f, %.6f] m\n', result.x(4:6))
fprintf('rest = %.6f m; tendon = %.6f m\n', result.x(7:8))
fprintf('objective = %.9g\n', result.objective)
fprintf('minimum torque margin = %+.2f%%\n', ...
    100*result.minTorqueMarginFraction)
fprintf('mean remaining shortfall = %.2f%%\n', ...
    100*result.meanTorqueShortfallFraction)
fprintf('reach c = %+.3f mm\n', 1000*result.pred.cRestLength)
fprintf('tibia/femur minimum clearances = %.3f / %.3f mm\n', ...
    1000*result.collisionInfo.minClearanceTibia, ...
    1000*result.collisionInfo.minClearanceFemur)
fprintf('collision/reach/strain/all feasible = %d / %d / %d / %d\n', ...
    result.collisionFeasible, result.reachFeasible, ...
    result.strainFeasible, result.allFeasible)
end

function style = publicationStyle()
run('Colors.m')
style.originalColor = [0.4 0.4 0.4];
style.optimizedColor = c{5};
style.candidateColor = c{3};
style.limitColor = c{7};
style.humanColor = '#000000';
style.fontName = 'Arial';
style.axesFontSize = 10;
style.titleFontSize = 12;
style.legendFontSize = 8;
style.axesLineWidth = 2;
style.originalLineWidth = 2;
style.optimizedLineWidth = 2.5;
style.humanLineWidth = 4;
style.tickLength = [0.025 0.05];
style.figurePosition = [2 2 14 10.5];
end

function plotBaselineResults(original, optimized, ctx, style)
humanTorque = -ctx.humanTorqueAbs;

figure('Name','Flexor Torque','Color','w','Units','centimeters', ...
    'Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,ctx.phiD,original.pred.TorqueZ,'--', ...
    'Color',style.originalColor,'LineWidth',style.originalLineWidth, ...
    'DisplayName','Original BPA')
plot(ax,ctx.phiD,optimized.pred.TorqueZ,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName','Optimized BPA')
plot(ax,ctx.humanAngleD,humanTorque,':','Color',style.humanColor, ...
    'LineWidth',style.humanLineWidth,'DisplayName','Human target')
finishAxes(ax,'Flexor Torque','Torque, N\cdotm',ctx,style,true)

figure('Name','Muscle Length','Color','w','Units','centimeters', ...
    'Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,ctx.phiD,original.pred.activeLength,'--', ...
    'Color',style.originalColor,'LineWidth',style.originalLineWidth, ...
    'DisplayName','Original BPA')
plot(ax,ctx.phiD,optimized.pred.activeLength,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName','Optimized BPA')
finishAxes(ax,'Muscle Length, L_m','Muscle Length, m',ctx,style,true)

figure('Name','Relative Strain','Color','w','Units','centimeters', ...
    'Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,ctx.phiD,original.pred.relativeStrain,'--', ...
    'Color',style.originalColor,'LineWidth',style.originalLineWidth, ...
    'DisplayName','Original BPA')
plot(ax,ctx.phiD,optimized.pred.relativeStrain,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName','Optimized BPA')
yline(ax,ctx.maxRelStrain,':','Color',style.limitColor, ...
    'LineWidth',style.originalLineWidth,'HandleVisibility','off')
finishAxes(ax,'Relative Strain','Relative Strain, \epsilon/K_{MAX}', ...
    ctx,style,true)

figure('Name','Moment Arm','Color','w','Units','centimeters', ...
    'Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,ctx.phiD,original.pred.momentArm,'--', ...
    'Color',style.originalColor,'LineWidth',style.originalLineWidth, ...
    'DisplayName','Original BPA')
plot(ax,ctx.phiD,optimized.pred.momentArm,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName','Optimized BPA')
finishAxes(ax,'Moment Arm','Moment Arm, m',ctx,style,true)

figure('Name','Torque Margin Fraction','Color','w', ...
    'Units','centimeters','Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,ctx.phiD,100*optimized.torqueMarginFraction,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName','Optimized BPA')
yline(ax,0,':','Color',style.humanColor, ...
    'LineWidth',style.originalLineWidth,'HandleVisibility','off')
finishAxes(ax,'BPA Torque Margin Relative to Human', ...
    'Torque Margin, %',ctx,style,true)
end

function plotShiftSweep(rawTable, correctedTable, style)
figure('Name','Negative-x Shift Sweep','Color','w', ...
    'Units','centimeters','Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,rawTable.Shift_mm,rawTable.MinTorqueMargin_pct,'--', ...
    'Color',style.originalColor,'LineWidth',style.originalLineWidth, ...
    'DisplayName','Raw shift')
plot(ax,correctedTable.Shift_mm, ...
    correctedTable.MinTorqueMargin_pct,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName','Tendon-corrected reach')
yline(ax,0,':','Color',style.humanColor, ...
    'LineWidth',style.originalLineWidth,'HandleVisibility','off')
set(ax,'XDir','reverse')
finishAxes(ax,'Equal Negative-x Shift Sweep', ...
    'Worst-case Torque Margin, %',[],style,true)
xlabel(ax,'Equal p1 and p2 x Shift, mm','FontName',style.fontName, ...
    'FontSize',style.axesFontSize,'FontWeight','bold')
end

function plotCandidateTorque(baseline, candidate, ctx, style, name)
if ~candidate.ok
    return
end
figure('Name',name,'Color','w','Units','centimeters', ...
    'Position',style.figurePosition)
ax = gca; hold(ax,'on')
plot(ax,ctx.phiD,baseline.pred.TorqueZ,'--', ...
    'Color',style.originalColor,'LineWidth',style.originalLineWidth, ...
    'DisplayName','Reported optimum')
plot(ax,ctx.phiD,candidate.pred.TorqueZ,'-', ...
    'Color',style.optimizedColor,'LineWidth',style.optimizedLineWidth, ...
    'DisplayName',char(candidate.label))
plot(ax,ctx.humanAngleD,-ctx.humanTorqueAbs,':', ...
    'Color',style.humanColor,'LineWidth',style.humanLineWidth, ...
    'DisplayName','Human target')
finishAxes(ax,name,'Torque, N\cdotm',ctx,style,true)
end

function finishAxes(ax,titleText,yLabelText,ctx,style,showLegend)
set(ax,'FontName',style.fontName,'FontSize',style.axesFontSize, ...
    'FontWeight','bold','LineWidth',style.axesLineWidth, ...
    'Box','off','XMinorTick','on','YMinorTick','on', ...
    'TickLength',style.tickLength,'GridLineStyle','none')
if ~isempty(ctx)
    xlim(ax,[min(ctx.phiD),max(ctx.phiD)])
end
xlabel(ax,'\theta_k, °','Interpreter','tex', ...
    'FontName',style.fontName,'FontSize',style.axesFontSize, ...
    'FontWeight','bold')
ylabel(ax,yLabelText,'Interpreter','tex', ...
    'FontName',style.fontName,'FontSize',style.axesFontSize, ...
    'FontWeight','bold')
title(ax,titleText,'Interpreter','tex','FontName',style.fontName, ...
    'FontSize',style.titleFontSize,'FontWeight','bold')
if showLegend
    lg = legend(ax,'Location','best');
    set(lg,'FontName',style.fontName,'FontSize',style.legendFontSize, ...
        'FontWeight','bold','Box','off')
end
grid(ax,'off')
end
