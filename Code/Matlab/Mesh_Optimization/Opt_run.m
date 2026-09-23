clear
clear functions
clc
rehash

% Marks a full optimizer pass. The result save clears it just before
% saving, so rerunning sections from a loaded result mat (which therefore
% lacks liveRun) cannot mint a new dated mat from an old design.
liveRun = true;


ctx = buildKneeFlexorContext20mm();
geo = ctx.geo;
idxP2 = 4:6;


obj = @(x) objective_KneeFlexor20mm(x, ctx);
objconstr = @(x) objconstrExclusion(x, obj, geo, ctx, idxP2);
nonlcon = @(x) nonlconExclusion(x, geo, ctx, idxP2);


rng default


if isempty(gcp('nocreate'))
    parpool;
end


%% Global search
% Seed the surrogate construction phase with known-good designs (2026-09-17:
% identical-config runs landed at -39% vs -4.3% margin depending on which
% basin the quasi-random initial phase found). Rows are full 8-variable
% points from Results, clipped into [lb, ub] for safety.
seedFiles = { ...
    'Bifemsh_20mm_Result_20260910_1234.mat', ...  % +1.02% record (Xi3 0.294 ctx)
    'Bifemsh_20mm_Result_20260917_1927.mat', ...  % A1 run
    'Bifemsh_20mm_Result_20260917_2004.mat', ...  % B1 run, -4.29% basin
    'Bifemsh_20mm_Result_20260918_1052.mat', ...  % seeded pick-1 run, +5.01%
    'Bifemsh_20mm_Result_20260918_1138.mat'};     % seeded pick-107 run, +3.61%
initPts = ctx.x0(:).';
for kSeed = 1:numel(seedFiles)
    Sseed = load(fullfile(fileparts(mfilename('fullpath')), ...
        'Results', seedFiles{kSeed}), 'xBest');
    initPts = [initPts; Sseed.xBest(:).']; 
end

% Ben, 2026-09-21: additional corner seed targeting a LOWER-TIBIA p2 --
% p1x = lb, p1y = ub, p2x = lb, p2y = lb, tendon = lb, and rest sized
% from the 5-deg extension pose so the series-length constraint starts
% satisfied: rest = pathLength0(idxExtension) - tendon - 2*fitting.
% (With Xi0 > 0 this seed sits Xi0 inside the cRestLength boundary.)
% p1z/p2z keep the x0 values. Rest does not affect route geometry, so
% the placeholder rest (lb) only has to let the predictor return
% pathLength0.
xCorner = [ctx.lb(1), ctx.ub(2), ctx.x0(3), ...
           ctx.lb(4), ctx.lb(5), ctx.x0(6), ctx.lb(7), ctx.lb(8)];
try
    predCorner = predictKneeFlexor20mm(xCorner, ctx);
    dCorner = predCorner.pathLength0(ctx.idxExtension);
    xCorner(7) = dCorner - xCorner(8) - 2*ctx.fitting;
    initPts = [initPts; xCorner(:).']; 
    fprintf(['Corner seed (lower-tibia p2): rest = %.4f m from ' ...
        'pathLength0(5 deg) = %.4f m.\n'], xCorner(7), dCorner)
catch
    fprintf(['Corner seed SKIPPED: predictor failed at the corner ' ...
        'geometry.\n'])
end

initPts = min(max(initPts, ctx.lb(:).'), ctx.ub(:).');
fprintf('Seeding surrogateopt with %d initial points.\n', size(initPts,1))

% Ben, 2026-09-21: 7000 evals was too many -- runs plateau by ~200
% function evaluations. 1000 leaves ample room past the plateau;
% patternsearch below still refines the winner.
optsG = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 1000, ...
    'MinSampleDistance', 0.001, ...
    'ConstraintTolerance', 1e-6, ...
    'InitialPoints', initPts);


[xG, fG, exitflagG, outputG] = surrogateopt( ...
    objconstr, ctx.lb, ctx.ub, optsG);


%% Pattern-search refinement
optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 15000, ...
    'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-6, ...
    'PlotFcn', {@psplotbestf, ...
                @psplotfuncount, ...
                @psplotmeshsize, ...
                @psplotmaxconstr}, ...
    'OutputFcn', @patternProgress);


[xBest, fBest, exitflagP, outputP] = patternsearch( ...
    obj, xG, [], [], [], [], ctx.lb, ctx.ub, nonlcon, optsP);


%% Run with adjusted seed
% Section commented out if unused.


% Display-from-mat workflow: this section is often rerun after loading a
% dated result mat (xBest, ctx, fBest, predBest, ...). The run-setup
% variables below exist only after the optimizer stages of a full run;
% rebuild them from the mat's ctx so the section is self-sufficient.
% Each guard is a no-op during a full Opt_run.
if ~exist('geo', 'var')
    geo = ctx.geo;
end
if ~exist('idxP2', 'var')
    idxP2 = 4:6;
end
if ~exist('nonlcon', 'var')
    nonlcon = @(x) nonlconExclusion(x, geo, ctx, idxP2);
end
if ~exist('optsP', 'var')
    optsP = optimoptions('patternsearch', ...
        'Display', 'iter', ...
        'UseParallel', true, ...
        'MaxFunctionEvaluations', 15000, ...
        'MeshTolerance', 1e-4, ...
        'StepTolerance', 1e-4, ...
        'ConstraintTolerance', 1e-6);
end

%Listing constraints (can put this before or after the xSeed ... cSeed
%block)
constraintNames = { ...
    'p2 exclusion', ...
    'tibia collision', ...
    'femur collision', ...
    'series length', ...
    'wrap', ...
    'tendon/wrap length', ...
    'relative strain'};

xSeed = [xBest(1:3), ...
         xBest(4:6), ...
         xBest(7), xBest(8)];

% Result mats saved before 2026-09-21 carry a ctx WITHOUT the soft
% BPA-2 anchor knobs (that is what "Unrecognized field name
% bpa2AnchorWeight" below used to mean).  Backfill the builder defaults
% so this section runs from any mat; keep values in sync with
% buildKneeFlexorContext20mm.
if ~isfield(ctx, 'bpa2AnchorWeight')
    ctx.bpa2AnchorWeight = 10;
    ctx.bpa2AnchorNorm   = 0.01;
    ctx.bpa2TargetP1     = [];
    ctx.bpa2TargetP2     = [];
end

% Soft BPA-2 anchor for this refinement (Ben, 2026-09-21): pull the
% search toward keeping the second BPA's DERIVED endpoints near this
% design's values -- the Q-angle / aLDFA build line, easier to build.
% Score punishment only (ctx.bpa2AnchorWeight, 0 disables), never a
% constraint; scaled so a 1 cm drift on both ends (20) stays far below
% 1 N m of worst-angle torque shortfall (~4300).
[ctx.bpa2TargetP1, ctx.bpa2TargetP2] = ...
    flexorBpa2Endpoints20mm(xSeed(1:3), xSeed(4:6));
fprintf(['BPA-2 anchor targets (soft, weight %g): ' ...
    'p1{2} = [% .6f % .6f % .6f] m, pEnd{2} = [% .6f % .6f % .6f] m\n'], ...
    ctx.bpa2AnchorWeight, ctx.bpa2TargetP1, ctx.bpa2TargetP2)

Jseed = objective_KneeFlexor20mm(xSeed, ctx);
cSeed = nonlcon(xSeed);

fprintf('\nTwo-route constraint values (BPA 1 + BPA 2, worst of both):\n');
for i = 1:numel(cSeed)
    fprintf('%-22s = %+0.9f\n', constraintNames{i}, cSeed(i));
end


boundViolation = max([ ...
    ctx.lb(:) - xSeed(:); ...
    xSeed(:) - ctx.ub(:)]);

fprintf('Earlier-design objective = %.6g\n', Jseed);
fprintf('Maximum nonlinear constraint = %.6g\n', max(cSeed));
fprintf('Maximum bound violation = %.6g\n', boundViolation);

% Adjust only coordinates outside the existing bounds.
xStart = min(max(xSeed(:), ctx.lb(:)), ctx.ub(:)).';

% Use the current context and objective.
objRefine = @(x) objective_KneeFlexor20mm(x, ctx);
conRefine = @(x) nonlcon(x);

% Keep the new result separate from your current xBest.
[xRefined, fRefined, exitRefined] = patternsearch( ...
    objRefine, xStart, ...
    [], [], [], [], ctx.lb, ctx.ub, conRefine, optsP);

cRefined = conRefine(xRefined);

fprintf('Refined objective = %.6g\n', fRefined);
fprintf('Maximum nonlinear constraint = %.6g\n', max(cRefined));
fprintf('Maximum bound violation = %.6g\n', max([ ...
    ctx.lb(:) - xRefined(:); ...
    xRefined(:) - ctx.ub(:)]));
fprintf('Exit flag = %d\n', exitRefined);

xBest = xRefined;
fBest = fRefined;

% How far the refinement moved the second BPA's endpoints from the
% anchor targets above.
if isfield(ctx, 'bpa2TargetP1') && ~isempty(ctx.bpa2TargetP1)
    [p1Bfinal, p2Bfinal] = flexorBpa2Endpoints20mm(xBest(1:3), xBest(4:6));
    fprintf(['BPA-2 anchor drift after refinement: p1{2} moved %.1f mm, ' ...
        'pEnd{2} moved %.1f mm\n'], ...
        100*norm(p1Bfinal(:) - ctx.bpa2TargetP1(:)), ...
        100*norm(p2Bfinal(:) - ctx.bpa2TargetP2(:)));
end

predBest = predictKneeFlexor20mm(xBest, ctx);
[cBest, ~] = nonlcon(xBest);
relativeContractionBest = predBest.relativeContraction;

%% Evaluate original and optimized designs
% Display-from-mat workflow: load a dated result mat (xBest, ctx, fBest,
% predBest, ...) then run this section. The guards below supply the
% run-setup variables this section reads that the mat does not carry;
% each is a no-op during a full Opt_run.
if ~exist('geo', 'var')
    geo = ctx.geo;
end
if ~exist('idxP2', 'var')
    idxP2 = 4:6;
end


predOriginal = predictOriginalKneeFlexor20mm(ctx);
predX0 = predictKneeFlexor20mm(ctx.x0, ctx);
predBest = predictKneeFlexor20mm(xBest, ctx);
[cCollision, ~, collisionInfo] = nonlconExclusion( ...
    xBest, geo, ctx, idxP2);


if ~predOriginal.ok
    error('Original-design prediction failed: %s', predOriginal.failReason)
end


if ~predX0.ok
    error('Original optimizer-guess prediction failed: %s', predX0.failReason)
end


if ~predBest.ok
    error('Optimized-design prediction failed: %s', predBest.failReason)
end


if exist('optsP', 'var')
    constraintTolerance = optsP.ConstraintTolerance;
else
    constraintTolerance = 1e-6;  % display-from-mat run: options not loaded
end
collisionFeasible = all(cCollision <= constraintTolerance);


% Do not present an infeasible final iterate as an optimized solution.
if ~collisionFeasible
    error(['Pattern search returned an infeasible point. Maximum nonlinear ' ...
        'constraint violation = %.9g m; tolerance = %.9g m.'], ...
        max(cCollision), constraintTolerance)
end


% Three-row design matrices: p1 is in femur; pWrap and pEnd are in t1.
pInitialWrapped = [predX0.p1; predX0.pWrap; predX0.pEnd];
pOptimized = [predBest.p1; predBest.pWrap; predBest.pEnd];
pChanged = pOptimized - pInitialWrapped;


% Save every optimizer-specific input consumed by Knee_Flexor_Data_20mm.
routeCtx = struct;


routeCtx.N = ctx.N;
routeCtx.phiD = ctx.phiD;
routeCtx.T_t1_f = ctx.T_t1_f;
routeCtx.T_ICR_t1 = ctx.T_ICR_t1;
routeCtx.wrapPointT1XY = ctx.wrapPointT1XY;
routeCtx.wrapAngleToleranceD = ctx.wrapAngleToleranceD;
routeCtx.BPAcount = ctx.BPAcount;
routeCtx.bpaRadiusMode = ctx.bpaRadiusMode;


% Save the exact scalar radii or converged per-frame bpaR arrays used by
% the optimized prediction. wRap remains the independent Xi3 bend radius.
routeCtx.geo = predBest.geo;


Xi3 = ctx.Xi3;


% Dated FULL-WORKSPACE result capture into Results (Ben directive,
% 2026-09-20: bare save so any driver section reruns from the loaded
% mat); does not overwrite prior results. liveRun is set only by a full
% Opt_run pass and cleared before saving, so rerunning this section from
% a loaded result mat cannot mint a new dated mat from an old design.
% XiUsed documents the Xi the run used (redundant with ctx; kept for the
% display loaders). Uncomment for a real run.
if exist('liveRun', 'var')
    stamp = char(string(datetime('now'),'yyyyMMdd_HHmm'));
    resDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
    resultFile = fullfile(resDir, sprintf('Bifemsh_20mm_Result_%s.mat', stamp));
    XiUsed = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
    clear liveRun
    save(resultFile)
    fprintf('Saved %s\n', resultFile)
end


%% Full-extension/full-flexion muscle-length and travel calculations
[~, idxFullExtension] = max(ctx.phiD);  % +10 deg normal calculation limit
[~, idxFullFlexion] = min(ctx.phiD);    % -120 deg normal calculation limit


LmExtension = predBest.activeLength(idxFullExtension);
LmFlexion = predBest.activeLength(idxFullFlexion);
deltaLmSigned = LmFlexion - LmExtension;
deltaLmAbsolute = abs(deltaLmSigned);


% Maximum physical BPA shortening from resting to fully contracted length.
maxContractionTravel = predBest.rest - predBest.kmax;


%% Display results
fprintf('\n========== OPTIMIZED DESIGN VALUES ==========\n')
fprintf('surrogateopt exitflag    = %d\n', exitflagG)
fprintf('patternsearch exitflag   = %d\n', exitflagP)
fprintf('surrogate evaluations   = %d\n', outputG.funccount)
fprintf('pattern evaluations     = %d\n', outputP.funccount)
if exist('fG', 'var')
    fprintf('fG                       = %.9g\n', fG)
else
    fprintf('fG                       = (not in result mat)\n')
end
fprintf('fBest                    = %.9g\n', fBest)


fprintf('\np (original), m; rows = [p1; pWrap; pEnd]:\n')
fprintf('[% .6f, % .6f, % .6f;\n', pInitialWrapped(1,:))
fprintf(' % .6f, % .6f, % .6f;\n', pInitialWrapped(2,:))
fprintf(' % .6f, % .6f, % .6f]\n', pInitialWrapped(3,:))


fprintf('\np (optimized), m; rows = [p1; pWrap; pEnd]:\n')
fprintf('[% .6f, % .6f, % .6f;\n', pOptimized(1,:))
fprintf(' % .6f, % .6f, % .6f;\n', pOptimized(2,:))
fprintf(' % .6f, % .6f, % .6f]\n', pOptimized(3,:))


fprintf('\np (changed), m; optimized minus original:\n')
fprintf('[%+.6f, %+.6f, %+.6f;\n', pChanged(1,:))
fprintf(' %+.6f, %+.6f, %+.6f;\n', pChanged(2,:))
fprintf(' %+.6f, %+.6f, %+.6f]\n', pChanged(3,:))


fprintf('wrap y-z line fraction at extension  = %.6f\n', ...
    predBest.routeInfo.wrapYZFraction)
fprintf('wrap release found                  = %d\n', ...
    predBest.routeInfo.releaseFound)
if predBest.routeInfo.releaseFound
    fprintf('first inactive wrap angle           = %+.6f deg\n', ...
        predBest.routeInfo.releaseAngleD)
else
    fprintf('first inactive wrap angle           = NONE IN MODELED RANGE\n')
end

if predBest.BPAcount == 2
    fprintf('\nBPA 2 route (same-side p1, mirrored pEnd, 2026-09-21):\n')
    fprintf('p1{2}, femur frame          = [% .6f, % .6f, % .6f] m\n', ...
        predBest.p1B)
    fprintf('pEnd{2}, t1 frame           = [% .6f, % .6f, % .6f] m\n', ...
        predBest.p2B)
    fprintf('BPA 2 wrap release found    = %d\n', ...
        predBest.routeInfo2.releaseFound)
    if predBest.routeInfo2.releaseFound
        fprintf('BPA 2 first inactive wrap angle = %+.6f deg\n', ...
            predBest.routeInfo2.releaseAngleD)
    else
        fprintf('BPA 2 first inactive wrap angle = NONE IN MODELED RANGE\n')
    end
end


fprintf('\nBPA and tendon lengths:\n')
fprintf('number of parallel BPAs             = %d\n', predBest.BPAcount)
fprintf('BPA radius source                    = %s\n', predBest.bpaRadiusMode)
fprintf('BPA radius range                     = %.6f to %.6f m\n', ...
    min(predBest.bpaRadius), max(predBest.bpaRadius))
if predBest.bpaRadiusMode == "bpaR"
    fprintf('BPA radius fixed-point iterations    = %d\n', ...
        predBest.bpaRadiusIteration)
    fprintf('BPA radius update converged          = %d\n', ...
        predBest.bpaRadiusConverged)
    fprintf('final maximum radius change          = %.9g m\n', ...
        predBest.bpaRadiusChange)
end
fprintf('rest length                         = %.6f m\n', predBest.rest)
fprintf('fully contracted length, kmax       = %.6f m\n', predBest.kmax)
fprintf('maximum contraction fraction, KMAX  = %.6f\n', predBest.KMAX)
fprintf('maximum BPA contraction travel      = %.6f m\n', maxContractionTravel)
fprintf('tendon length                       = %.6f m\n', predBest.tendon)


fprintf('\nMuscle length over normal RoM:\n')
fprintf('Lm at full extension (%+.6f deg) = %.6f m\n', ...
    ctx.phiD(idxFullExtension), LmExtension)
fprintf('Lm at full flexion   (%+.6f deg) = %.6f m\n', ...
    ctx.phiD(idxFullFlexion), LmFlexion)
fprintf('Lm flexion minus extension         = %+.6f m\n', deltaLmSigned)
fprintf('absolute Lm difference             = %.6f m\n', deltaLmAbsolute)


fprintf('\nRest-length check:\n')
fprintf('restLmt           = %.6f m\n', predBest.restLmt)
fprintf('extensionDistance = %.6f m\n', predBest.extensionDistance)
fprintf('cRestLength       = %.6f m\n', predBest.cRestLength)
fprintf('max Contraction/KMAX = %.6f\n', max(predBest.relativeContraction))


offAxisExcessX = max(0, ...
    abs(predBest.TorqueX) - abs(ctx.originalTorqueX));
offAxisExcessY = max(0, ...
    abs(predBest.TorqueY) - abs(ctx.originalTorqueY));
fprintf('\nOff-axis torque relative to original BPA:\n')
fprintf('maximum |original Tx|        = %.6f N m\n', ...
    max(abs(ctx.originalTorqueX)))
fprintf('maximum |optimized Tx|       = %.6f N m\n', ...
    max(abs(predBest.TorqueX)))
fprintf('maximum pointwise |Tx| excess= %.6f N m\n', ...
    max(offAxisExcessX))
fprintf('maximum |original Ty|        = %.6f N m\n', ...
    max(abs(ctx.originalTorqueY)))
fprintf('maximum |optimized Ty|       = %.6f N m\n', ...
    max(abs(predBest.TorqueY)))
fprintf('maximum pointwise |Ty| excess= %.6f N m\n', ...
    max(offAxisExcessY))


fprintf('\nCollision check at %.6f deg:\n', collisionInfo.angleD)
fprintf('constraint tolerance       = %.9g m\n', constraintTolerance)
fprintf('all collision constraints feasible = %d\n', collisionFeasible)
fprintf('p2 exclusion constraint    = %.6f m\n', cCollision(1))
fprintf('tibia collision constraint = %.6f m\n', cCollision(2))
fprintf('femur collision constraint = %.6f m\n', cCollision(3))
fprintf('series-length constraint   = %.6f m\n', cCollision(4))
fprintf('routes checked             = %d (BPA 2 = same-side p1, mirrored pEnd)\n', ...
    collisionInfo.routeCount)
fprintf('binding route (worst)      = BPA %d\n', collisionInfo.bindingRoute)
fprintf('min tibia clearance by route = %.6f / %.6f m\n', ...
    collisionInfo.minClearanceTibiaByRoute)
fprintf('min femur clearance by route = %.6f / %.6f m\n', ...
    collisionInfo.minClearanceFemurByRoute)
fprintf('radius source              = %s\n', collisionInfo.bpaRadiusMode)
fprintf('pWrap radius bpaRb         = %.6f m\n', collisionInfo.bpaRb)
fprintf('collision radius bpaRs     = %.6f m\n', collisionInfo.bpaRs)
fprintf('Xi3 wrap radius wRap       = %.6f m\n', collisionInfo.wRap)
fprintf('tendon radius              = %.6f m\n', collisionInfo.tendonRadius)
fprintf('optimized tendon length    = %.6f m\n', collisionInfo.tendon)
fprintf('current muscle length Lm   = %.6f m\n', ...
    collisionInfo.currentMuscleLength)
fprintf('Lm + two fittings          = %.6f m\n', ...
    collisionInfo.bpaFittingsLengthChecked)
fprintf('minimum tibia clearance    = %.6f m\n', ...
    collisionInfo.minClearanceTibia)
fprintf('minimum femur clearance    = %.6f m\n', ...
    collisionInfo.minClearanceFemur)
fprintf('extra femur clearance      = %.6f m\n', ...
    collisionInfo.requiredClearance)
fprintf('worst tibia region         = %s\n', ...
    collisionInfo.worstTibiaRegion)
fprintf('worst tibia component      = %s (radius %.6f m)\n', ...
    collisionInfo.worstTibiaComponent, collisionInfo.worstTibiaRadius)
fprintf('worst tibia center, t1     = [%.6f, %.6f, %.6f] m\n', ...
    collisionInfo.worstTibiaCenterT1)
fprintf('worst femur component      = %s (radius %.6f m)\n', ...
    collisionInfo.worstFemurComponent, collisionInfo.worstFemurRadius)
fprintf('worst femur center, femur  = [%.6f, %.6f, %.6f] m\n', ...
    collisionInfo.worstFemurCenterFemur)
fprintf('=============================================\n')


fprintf('min Contraction/KMAX = %+.6f\n', ...
    min(predBest.Contraction)/predBest.KMAX);
fprintf('max Contraction/KMAX = %+.6f\n', ...
    max(predBest.Contraction)/predBest.KMAX);


%% Plot original and optimized results in separate figure windows
humanTorque = -ctx.humanTorqueAbs;  % Flexor torque remains negative.


% Use the established accessible project palette.  Colors.m supplies c for
% line colors and d for RGB scatter-marker colors; these figures use c.
run('Colors.m')
originalColor = [0.4 0.4 0.4];
optimizedColor = c{5};
humanColor = '#000000';
limitColor = c{7};


fontName = 'Arial';
axesFontSize = 10;
titleFontSize = 12;
legendFontSize = 8;
axesLineWidth = 2;
originalLineWidth = 2;
optimizedLineWidth = 2.5;
humanLineWidth = 4;
tickLength = [0.025 0.05];
figurePosition = [2 2 14 10.5];  % centimeters; 14 cm publication width
xLimits = [min(ctx.phiD), max(ctx.phiD)];


%% Flexor torque
figure('Name', 'Flexor Torque', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predOriginal.TorqueZ, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(ax, ctx.phiD, predBest.TorqueZ, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
plot(ax, ctx.humanAngleD, humanTorque, ':', ...
    'Color', humanColor, 'LineWidth', humanLineWidth, ...
    'DisplayName', 'Human target')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'Torque, N\cdotm', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'Flexor Torque', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Muscle length
figure('Name', 'Muscle Length', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predOriginal.activeLength, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(ax, ctx.phiD, predBest.activeLength, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'Muscle Length, m', 'FontName', fontName, ...
    'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'Muscle Length, L_m', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Strain definitions
figure('Name', 'Flexor Strain', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predBest.strain_f, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'strain_f, includes Xi3')
plot(ax, ctx.phiD, predBest.strain_p, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'strain_p, excludes Xi3')
plot(ax, ctx.phiD, predBest.Contraction, '-.', ...
    'Color', c{6}, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Contraction')
minimumLine = yline(ax, 0, ':', 'Minimum strain', ...
    'Color', humanColor, 'LineWidth', originalLineWidth, ...
    'HandleVisibility', 'off');
maximumLine = yline(ax, predBest.KMAX, ':', 'KMAX', ...
    'Color', limitColor, 'LineWidth', originalLineWidth, ...
    'HandleVisibility', 'off');
set([minimumLine, maximumLine], 'FontName', fontName, ...
    'FontSize', legendFontSize, 'FontWeight', 'bold')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'Strain', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'Flexor Strain Definitions', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Moment arm
figure('Name', 'Moment Arm', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predOriginal.momentArm, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(ax, ctx.phiD, predBest.momentArm, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'Moment Arm, m', 'FontName', fontName, ...
    'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'Moment Arm', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% X-axis torque relative to the original no-wrap BPA
figure('Name', 'X-axis Torque', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predOriginal.TorqueX, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(ax, ctx.phiD, predBest.TorqueX, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'T_x, N\cdotm', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'X-axis Torque', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Y-axis torque relative to the original no-wrap BPA
figure('Name', 'Y-axis Torque', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predOriginal.TorqueY, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(ax, ctx.phiD, predBest.TorqueY, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'T_y, N\cdotm', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'Y-axis Torque', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Torque margin fraction
% Positive means the BPA exceeds the required human flexor-torque
% magnitude. Negative means a remaining torque shortfall. The flexor
% torque curves themselves remain signed and negative; absolute values are
% used only here to form the magnitude ratio.
humanAbsAtRobotAngles = interp1( ...
    ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');
validHumanTorque = humanAbsAtRobotAngles > 100*eps;
torqueMarginFraction = nan(size(predBest.TorqueZ));
torqueMarginFraction(validHumanTorque) = ...
    abs(predBest.TorqueZ(validHumanTorque)) ./ ...
    humanAbsAtRobotAngles(validHumanTorque) - 1;


fprintf('\nTorque margin relative to human target:\n')
fprintf('minimum margin fraction = %+.6f (%+.2f%%)\n', ...
    min(torqueMarginFraction(validHumanTorque)), ...
    100*min(torqueMarginFraction(validHumanTorque)))
fprintf('mean remaining shortfall = %.6f (%.2f%%)\n', ...
    mean(max(0, -torqueMarginFraction(validHumanTorque))), ...
    100*mean(max(0, -torqueMarginFraction(validHumanTorque))))


figure('Name', 'Torque Margin Fraction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, 100*torqueMarginFraction, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
zeroLine = yline(ax, 0, ':', 'Human target', ...
    'Color', humanColor, 'LineWidth', originalLineWidth, ...
    'HandleVisibility', 'off');
zeroLine.FontName = fontName;
zeroLine.FontSize = legendFontSize;
zeroLine.FontWeight = 'bold';
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'Torque Margin, %', 'FontName', fontName, ...
    'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, 'BPA Torque Margin Relative to Human', ...
    'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
grid(ax, 'off')


%% Local output function: concise pattern-search progress in Command Window
function [stop, options, optchanged] = patternProgress( ...
        optimValues, options, flag)


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
