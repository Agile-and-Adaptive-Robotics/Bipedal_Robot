% Opt_run_pulley.m -- reverse-pulley (block-and-tackle) spin-off of Opt_run.m
% Author: Ben Bolen
% Date: 2026-09-24
% Description: Same driver structure as Opt_run.m (liveRun flag, context
% from buildKneeFlexorContext20mm, surrogateopt seeded from Results mats +
% corner seed, patternsearch refinement, adjusted-seed section, display
% sections, dated FULL-WORKSPACE bare save behind liveRun) with the
% evaluation pipeline swapped to MonoPam_pulley: each BPA route is wrapped
% in a MonoPam_pulley object and the objective's torque requirement is
% computed on the INSERTION torque Torque_ins so the tackle's 1/G force
% loss is priced into the search. The result mat is
% Bifemsh_20mm_Result_pulley_<stamp>.mat.
%
% Transmission config (Ben, 2026-09-24 -- run-level, NOT free scalars):
%   ctx.pulleyBPACount    BPAs in parallel feeding the tackle (the rig: 2).
%   ctx.tackleLineParts   rope parts supporting the moving block = travel
%                         gain G; a four-wheel tackle rigges to at most 4.
%                         1 = pulley disabled (straight tendon). G = 1
%                         reproduces the straight-tendon baseline EXACTLY
%                         on the sanity/smoke geometry (non-tendon rows
%                         orientation-constant; Opt_sanity_pulley.m), and
%                         only approximately on the real flexor route
%                         (see the note at the ctx.pulley* knobs below).
%   ctx.pulleyRoutingMode 'moving_via' (default; the tendon exit is
%                         body-fixed to the proximal body and u_t rotates
%                         with the joint) or 'bowden' (housing anchor on
%                         the tibia; u_t fixed in the tibia frame).
%   ctx.pulleyExitIndex   row of the tendon exit point on the proximal
%                         body (default 1: with the flexor route's
%                         CrossPoint = 2 the crossed segment IS the
%                         tendon line, as the sanity gate asserts).
%   ctx.bowdenBossDia /   Bowden packaging min-clearance envelope,
%   ctx.bowdenRunClearance standard Shimano-type parts (M7 x 1.0 mm barrel
%                         adjuster ~7 mm boss; 5 mm brake housing + 5 mm
%                         bracket plate); reported in the display block.
% These are FIXED per optimization run; the discrete comparison across
% configs is the outer for-loop below (OPT_PULLEY_CONFIGS=1) producing one
% dated mat per config -- NOT integer-constrained nonlinear programming.
%
% Spin-off differences from Opt_run.m (all deliberate):
%   * ctx.optimizePulleyGain (default false): when false the problem is
%     the IDENTICAL 8-variable search of Opt_run (bounds, constraint
%     count, and structure); when true G joins the search as a 9th design
%     variable within ctx.pulleyGainBounds.
%   * The exclusion constraint vector gains an 8th entry: the
%     route-to-route clearance constraint. The moving via sweeps an
%     envelope as the joint rotates and the nPulleyBPA parallel BPAs each
%     occupy a body diameter, so the two routes' segment geometry
%     (including the swept via segments) must stay at least
%     2*bpaRs + ctx.routeClearanceMargin apart, worst over the
%     orientation sweep -- same worst-over-sweep style as the existing
%     exclusion constraints.
%   * ctx.geo.bpaRadiusMode is forced to "scalar" for this spin-off: the
%     radius-iteration coupling lives in the X3 pipeline, and scalar
%     radii keep nonlconExclusion reusable here without modification.
%   * MonoPam_pulley carries no Xi3/BendMeasure term (its constructor
%     mirrors MonoPamDataExplicit_balance), so the route-based strains
%     here exclude the bend-loss machinery of predictKneeFlexor20mm; the
%     run still saves Xi3 (unused by the pulley class) so Knee_Data.m can
%     rebuild either model from the same mat.
%   * Parallel-BPA footprint convention (biPulleySpecsFromOpenSim): the
%     bundle is placed SYMMETRICALLY about the original OpenSim line so
%     the first-order line of action and moment arm are preserved;
%     asymmetric placement shifts the line of action.
%
% SMOKE mode: set environment variable OPT_PULLEY_SMOKE=1 for tiny eval
% budgets, UseParallel false on both solvers, the SAME sanity geometry
% Opt_sanity_pulley.m asserts (straight-tendon 92-orientation route whose
% crossed segment is the tendon line), and no result save (liveRun stays
% undefined). Completes in minutes. Full optimization runs (parpool(10),
% the 1000/15000-eval budgets below) belong on easteregg2; this box caps
% the pool small.

clear
clear functions
clc
rehash

% Marks a full optimizer pass. The result save clears it just before
% saving, so rerunning sections from a loaded result mat (which therefore
% lacks liveRun) cannot mint a new dated mat from an old design.
% SMOKE mode never saves (liveRun stays undefined). strtrim guards against
% the cmd.exe "set VAR=1 && ..." trailing-space trap.
smokeMode = strcmp(strtrim(getenv('OPT_PULLEY_SMOKE')), '1');
if smokeMode
    fprintf(['OPT_PULLEY_SMOKE=1: tiny budgets, solvers on one core ' ...
        '(UseParallel false), the Opt_sanity_pulley geometry self-check, ' ...
        'no result save.\n'])
else
    % Marks a full optimizer pass. The result save clears it just before
    % saving, so rerunning sections from a loaded result mat (which
    % therefore lacks liveRun) cannot mint a new dated mat from an old
    % design.
    liveRun = true;
end

% Repo root and path setup (same block as Results/Knee_Flexor_Data_20mm.m
% and Opt_sanity_pulley.m): makes the driver runnable via matlab -batch
% and from fresh sessions; addpath is idempotent for Ben's saved paths.
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
% Append (do not prepend) so Code\Matlab keeps winning name collisions.
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');


ctx = buildKneeFlexorContext20mm();
geo = ctx.geo;
idxP2 = 4:6;

% --- Reverse-pulley spin-off knobs (Ben, 2026-09-24) --------------------
% G = 1 with nPulleyBPA = 1 reproduces the straight-tendon baseline on the
% sanity/smoke geometry (the regression identity Opt_sanity_pulley.m
% asserts: its route's non-tendon rows are orientation-constant). On the
% REAL flexor route the moving wrap segment is inside MuscleLength but
% outside the pulley's exit-to-insertion span, so the G = 1 solve is a
% close approximation there, NOT an identity: wrap-active frames can
% deviate from MonoPamDataExplicit_balance by a few percent of peak
% torque (measured ~3.3 N*m of ~59 N*m at the seed route).
ctx.pulleyBPACount = 2;          % BPAs in parallel feeding the tackle (the rig)
ctx.tackleLineParts = 1;         % travel gain G; 1 = pulley disabled
ctx.pulleyRoutingMode = 'moving_via';   % or 'bowden'
ctx.pulleyExitIndex = 1;         % tendon exit row; crossed segment = tendon line
ctx.bowdenBossDia = 0.007;       % M7 x 1.0 mm barrel-adjuster boss, m
ctx.bowdenRunClearance = 0.005;  % 5 mm brake housing + 5 mm bracket plate, m
ctx.routeClearanceMargin = 0.002;% route-to-route clearance margin, m
ctx.optimizePulleyGain = false;  % true: G becomes the 9th design variable
ctx.pulleyGainBounds = [1, 4];   % four-wheel tackle: travel gain <= 4

% Discrete config comparison (outer for-loop below): one dated mat per
% {nPulleyBPA} x {tackleLineParts} candidate. Default OFF: the loop runs
% the single ctx config and behaves exactly like Opt_run.
ctx.pulleyConfigSweep = strcmp(strtrim(getenv('OPT_PULLEY_CONFIGS')), '1');
ctx.pulleyBPACountList = [1, 2];       % physical rig = 2
ctx.tackleLinePartsList = [1, 2, 4];   % four-wheel tackle max = 4

% The spin-off runs scalar radii (see header) so nonlconExclusion is
% reused without modification.
ctx.bpaRadiusMode = "scalar";
ctx.geo.bpaRadiusMode = "scalar";
geo = ctx.geo;

% Extend the design vector with G only when the flag asks for it.
if ctx.optimizePulleyGain
    ctx.lb = [ctx.lb(:); ctx.pulleyGainBounds(1)];
    ctx.ub = [ctx.ub(:); ctx.pulleyGainBounds(2)];
    ctx.x0 = [ctx.x0(:); ctx.tackleLineParts];
end


obj = @(x) objective_KneeFlexorPulley(x, ctx);
objconstr = @(x) objconstrExclusionPulley(x, obj, geo, ctx, idxP2);
nonlcon = @(x) nonlconExclusionPulley(x, geo, ctx, idxP2);


rng default


if ~smokeMode
    % Full runs belong on easteregg2 (parpool(10)); this box caps small.
    if isempty(gcp('nocreate'))
        parpool(6);
    end
end


%% SMOKE: sanity-geometry check through the pulley class
% The SAME straight-tendon construction Opt_sanity_pulley.m asserts
% (92-orientation route whose crossing segment is the tendon line): the
% G = 1 regression identity vs MonoPamDataExplicit_balance and the G = 2
% equilibrium closure, at the sanity gate's tolerances.
if smokeMode
    smokeSanityCheck()
end


%% Transmission-config candidate list
% One dated mat per {nPulleyBPA} x {tackleLineParts} candidate behind
% OPT_PULLEY_CONFIGS=1; the default single-candidate list keeps the run
% identical in shape to Opt_run.
cfgBPACounts = ctx.pulleyBPACount;
cfgLineParts = ctx.tackleLineParts;
if ctx.pulleyConfigSweep
    % ROW vectors on purpose: a for-loop iterates over COLUMNS, so a
    % column-vector candidate list would run ONE iteration binding the
    % whole vector (and the vectors would then leak into %d fprintf
    % fields as garbage).
    cfgBPACounts = ctx.pulleyBPACountList;
    cfgLineParts = ctx.tackleLinePartsList;
    validateattributes(cfgBPACounts, {'numeric'}, ...
        {'row', 'positive', 'integer'});
    validateattributes(cfgLineParts, {'numeric'}, ...
        {'row', 'positive', 'integer'});
    fprintf(['Pulley config sweep ON: %d x %d = %d configs, one dated ' ...
        'mat per config.\n'], numel(cfgBPACounts), numel(cfgLineParts), ...
        numel(cfgBPACounts)*numel(cfgLineParts))
end
nCfg = numel(cfgBPACounts)*numel(cfgLineParts);
kCfg = 0;
cfgBestScore = inf;
cfgBestX = [];
cfgBestF = [];
cfgBestBPA = ctx.pulleyBPACount;
cfgBestG = ctx.tackleLineParts;

%% Discrete config loop: one full optimizer pass per transmission config
for nB = cfgBPACounts
    for nG = cfgLineParts
        kCfg = kCfg + 1;
        ctx.pulleyBPACount = nB;
        ctx.tackleLineParts = nG;
        fprintf(['\n===== Transmission config %d of %d: nPulleyBPA = %d, ' ...
            'tackleLineParts (G) = %d, mode = %s =====\n'], ...
            kCfg, nCfg, nB, nG, ctx.pulleyRoutingMode)

        % Each config mints exactly ONE dated mat (sweep mode only): the
        % loop body is the single full-run pass for that config.
        if ctx.pulleyConfigSweep && ~smokeMode
            liveRun = true;
        end


        %% Global search
        % Seed the surrogate construction phase with known-good designs
        % (2026-09-17: identical-config runs landed at -39% vs -4.3% margin
        % depending on which basin the quasi-random initial phase found).
        % Rows are full 8-variable points from Results, clipped into
        % [lb, ub] for safety. The design vector is config-independent, so
        % the seeds are computed once per config loop pass.
        seedFiles = { ...
            'Bifemsh_20mm_Result_20260910_1234.mat', ...  % +1.02% record (Xi3 0.294 ctx)
            'Bifemsh_20mm_Result_20260917_1927.mat', ...  % A1 run
            'Bifemsh_20mm_Result_20260917_2004.mat', ...  % B1 run, -4.29% basin
            'Bifemsh_20mm_Result_20260918_1052.mat', ...  % seeded pick-1 run, +5.01%
            'Bifemsh_20mm_Result_20260918_1138.mat'};     % seeded pick-107 run, +3.61%
        initPts = ctx.x0(:).';
        for kSeed = 1:numel(seedFiles)
            try
                Sseed = load(fullfile(fileparts(mfilename('fullpath')), ...
                    'Results', seedFiles{kSeed}), 'xBest');
                xSeedRow = Sseed.xBest(:).';
                if ctx.optimizePulleyGain
                    if numel(xSeedRow) < 9
                        xSeedRow = [xSeedRow, ctx.tackleLineParts];
                    end
                else
                    xSeedRow = xSeedRow(1:min(8, numel(xSeedRow)));
                end
                initPts = [initPts; xSeedRow];
            catch
                fprintf('Seed %s not found; skipped.\n', seedFiles{kSeed})
            end
        end

        % Ben, 2026-09-21: additional corner seed targeting a LOWER-TIBIA
        % p2 -- p1x = lb, p1y = ub, p2x = lb, p2y = lb, tendon = lb, and
        % rest sized from the 5-deg extension pose so the series-length
        % constraint starts satisfied:
        %   rest = pathLength0(idxExtension) - tendon - 2*fitting.
        % (With Xi0 > 0 this seed sits Xi0 inside the cRestLength
        % boundary.) p1z/p2z keep the x0 values.
        xCorner = [ctx.lb(1), ctx.ub(2), ctx.x0(3), ...
                   ctx.lb(4), ctx.lb(5), ctx.x0(6), ctx.lb(7), ctx.lb(8)];
        try
            predCorner = predictKneeFlexorPulley(xCorner, ctx);
            dCorner = predCorner.pathLength0(ctx.idxExtension);
            xCorner(7) = dCorner - xCorner(8) - 2*ctx.fitting;
            if ctx.optimizePulleyGain
                xCorner = [xCorner, ctx.pulleyGainBounds(1)];
            end
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
        % SMOKE caps the surrogate at 40 evals (<= 60) and single core.
        if smokeMode
            maxEvalsG = 40;
            useParallelG = false;
        else
            maxEvalsG = 1000;
            useParallelG = true;
        end
        optsG = optimoptions('surrogateopt', ...
            'Display', 'iter', ...
            'UseParallel', useParallelG, ...
            'MaxFunctionEvaluations', maxEvalsG, ...
            'MinSampleDistance', 0.001, ...
            'ConstraintTolerance', 1e-6, ...
            'InitialPoints', initPts);


        [xG, fG, exitflagG, outputG] = surrogateopt( ...
            objconstr, ctx.lb, ctx.ub, optsG);


        %% Pattern-search refinement
        % SMOKE caps patternsearch at 60 evals (<= 100), single core, and
        % drops the plot functions (no figures in -batch).
        if smokeMode
            maxEvalsP = 60;
            useParallelP = false;
            plotFcnP = {};
        else
            maxEvalsP = 15000;
            useParallelP = true;
            plotFcnP = {@psplotbestf, ...
                        @psplotfuncount, ...
                        @psplotmeshsize, ...
                        @psplotmaxconstr};
        end
        optsP = optimoptions('patternsearch', ...
            'Display', 'iter', ...
            'UseParallel', useParallelP, ...
            'MaxFunctionEvaluations', maxEvalsP, ...
            'MeshTolerance', 1e-4, ...
            'StepTolerance', 1e-4, ...
            'ConstraintTolerance', 1e-6, ...
            'PlotFcn', plotFcnP, ...
            'OutputFcn', @patternProgress);


        [xBest, fBest, exitflagP, outputP] = patternsearch( ...
            obj, xG, [], [], [], [], ctx.lb, ctx.ub, nonlcon, optsP);


        %% One dated mat per transmission config (sweep mode only)
        % Bare FULL-WORKSPACE save behind liveRun (Ben directive). The
        % non-sweep single-config run saves once, after the adjusted-seed
        % refinement below, exactly like Opt_run.
        if ctx.pulleyConfigSweep && exist('liveRun', 'var')
            stamp = char(string(datetime('now'),'yyyyMMdd_HHmm'));
            resDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
            resultFile = fullfile(resDir, sprintf( ...
                'Bifemsh_20mm_Result_pulley_nB%d_G%d_%s.mat', ...
                ctx.pulleyBPACount, ctx.tackleLineParts, stamp));
            XiUsed = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
            clear liveRun
            save(resultFile)
            fprintf('Saved %s\n', resultFile)
        end

        % Track the overall winner across configs (the adjusted-seed
        % refinement below runs on it).
        if fBest < cfgBestScore
            cfgBestScore = fBest;
            cfgBestX = xBest;
            cfgBestF = fBest;
            cfgBestBPA = ctx.pulleyBPACount;
            cfgBestG = ctx.tackleLineParts;
        end
    end
end

% Restore the winning config for the refinement and the final mat.
xBest = cfgBestX;
fBest = cfgBestF;
ctx.pulleyBPACount = cfgBestBPA;
ctx.tackleLineParts = cfgBestG;
if nCfg > 1
    fprintf(['\nConfig sweep winner: nPulleyBPA = %d, tackleLineParts ' ...
        '(G) = %d at fBest = %.6g.\n'], ctx.pulleyBPACount, ...
        ctx.tackleLineParts, fBest)
end
% Re-arm the result capture for the final winner mat: the per-config
% saves above consumed liveRun (sweep mode), and a fresh full-run pass
% section follows.
if ~smokeMode
    liveRun = true;
end


%% Run with adjusted seed
% Section commented out if unused.


% Display-from-mat workflow: this section is often rerun after loading a
% dated result mat (xBest, ctx, fBest, predBest, ...). The run-setup
% variables below exist only after the optimizer stages of a full run;
% rebuild them from the mat's ctx so the section is self-sufficient.
% Each guard is a no-op during a full Opt_run_pulley.
if ~exist('geo', 'var')
    geo = ctx.geo;
end
if ~exist('idxP2', 'var')
    idxP2 = 4:6;
end
if ~exist('nonlcon', 'var')
    nonlcon = @(x) nonlconExclusionPulley(x, geo, ctx, idxP2);
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
%block). c(8) is this spin-off's route-to-route clearance constraint.
constraintNames = { ...
    'p2 exclusion', ...
    'tibia collision', ...
    'femur collision', ...
    'series length', ...
    'wrap', ...
    'tendon/wrap length', ...
    'relative strain', ...
    'route-route clearance'};

xSeed = [xBest(1:3), ...
         xBest(4:6), ...
         xBest(7), xBest(8)];
if numel(xBest) > 8
    xSeed = [xSeed, xBest(9)];   % pulley gain rode along as variable 9
end

% Result mats saved before the pulley spin-off carry a ctx WITHOUT the
% pulley knobs. Backfill the builder defaults so this section runs from
% any mat; keep values in sync with the spin-off block at the top.
if ~isfield(ctx, 'pulleyBPACount')
    ctx.pulleyBPACount = 2;
end
if ~isfield(ctx, 'tackleLineParts')
    ctx.tackleLineParts = 1;
end
if ~isfield(ctx, 'pulleyRoutingMode')
    ctx.pulleyRoutingMode = 'moving_via';
end
if ~isfield(ctx, 'pulleyExitIndex')
    ctx.pulleyExitIndex = 1;
end
if ~isfield(ctx, 'bowdenBossDia')
    ctx.bowdenBossDia = 0.007;
end
if ~isfield(ctx, 'bowdenRunClearance')
    ctx.bowdenRunClearance = 0.005;
end
if ~isfield(ctx, 'routeClearanceMargin')
    ctx.routeClearanceMargin = 0.002;
end
if ~isfield(ctx, 'optimizePulleyGain')
    ctx.optimizePulleyGain = false;
end
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

Jseed = objective_KneeFlexorPulley(xSeed, ctx);
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
objRefine = @(x) objective_KneeFlexorPulley(x, ctx);
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

predBest = predictKneeFlexorPulley(xBest, ctx);
[cBest, ~] = nonlcon(xBest);
relativeContractionBest = predBest.relativeContraction;

%% Evaluate original and optimized designs
% Display-from-mat workflow: load a dated result mat (xBest, ctx, fBest,
% predBest, ...) then run this section. The guards below supply the
% run-setup variables this section reads that the mat does not carry;
% each is a no-op during a full Opt_run_pulley.
if ~exist('geo', 'var')
    geo = ctx.geo;
end
if ~exist('idxP2', 'var')
    idxP2 = 4:6;
end


predOriginal = predictOriginalKneeFlexor20mm(ctx);
predX0 = predictKneeFlexorPulley(ctx.x0(1:min(8,end)), ctx);
predBest = predictKneeFlexorPulley(xBest, ctx);
[cCollision, ~, collisionInfo] = nonlconExclusionPulley( ...
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


% Save every optimizer-specific input consumed by Knee_Data.
routeCtx = struct;


routeCtx.N = ctx.N;
routeCtx.phiD = ctx.phiD;
routeCtx.T_t1_f = ctx.T_t1_f;
routeCtx.T_ICR_t1 = ctx.T_ICR_t1;
routeCtx.wrapPointT1XY = ctx.wrapPointT1XY;
routeCtx.wrapAngleToleranceD = ctx.wrapAngleToleranceD;
routeCtx.BPAcount = ctx.BPAcount;
routeCtx.bpaRadiusMode = ctx.bpaRadiusMode;
routeCtx.pulleyBPACount = ctx.pulleyBPACount;
routeCtx.tackleLineParts = ctx.tackleLineParts;
routeCtx.pulleyRoutingMode = ctx.pulleyRoutingMode;
routeCtx.pulleyExitIndex = ctx.pulleyExitIndex;
routeCtx.bowdenBossDia = ctx.bowdenBossDia;
routeCtx.bowdenRunClearance = ctx.bowdenRunClearance;


% Save the exact scalar radii or converged per-frame bpaR arrays used by
% the optimized prediction.
routeCtx.geo = predBest.geo;


Xi3 = ctx.Xi3;


% Dated FULL-WORKSPACE result capture into Results (Ben directive,
% 2026-09-20: bare save so any driver section reruns from the loaded
% mat); does not overwrite prior results. liveRun is set only by a full
% Opt_run_pulley pass and cleared before saving, so rerunning this
% section from a loaded result mat cannot mint a new dated mat from an
% old design. XiUsed documents the Xi the run used (redundant with ctx;
% kept for the display loaders).
if exist('liveRun', 'var')
    stamp = char(string(datetime('now'),'yyyyMMdd_HHmm'));
    resDir = fullfile(fileparts(mfilename('fullpath')), 'Results');
    resultFile = fullfile(resDir, sprintf('Bifemsh_20mm_Result_pulley_%s.mat', stamp));
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
fprintf('\n========== OPTIMIZED DESIGN VALUES (PULLEY) ==========\n')
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


fprintf('\nReverse-pulley transmission (MonoPam_pulley):\n')
fprintf('nPulleyBPA (parallel BPAs, rig)     = %d\n', ...
    ctx.pulleyBPACount)
fprintf('route corridors modeled             = %d (each carries 1 BPA line into the shared tackle)\n', ...
    predBest.BPAcount)
fprintf('tackleLineParts = travel gain G     = %d\n', ...
    predBest.pulleyConfig.tackleLineParts)
if ctx.optimizePulleyGain
    fprintf('gain was OPTIMIZED (variable 9)     = %.6f, bounds [%g, %g]\n', ...
        xBest(9), ctx.pulleyGainBounds(1), ctx.pulleyGainBounds(2))
end
fprintf('routing mode                        = %s\n', ...
    predBest.pulleyConfig.routingMode)
fprintf('tackle exit row index               = %d\n', ...
    predBest.pulleyConfig.pulleyExitIndex)
fprintf(['Bowden envelope (report only in %s mode): boss dia = %.1f mm, ' ...
    'run clearance = %.1f mm\n'], predBest.pulleyConfig.routingMode, ...
    1000*predBest.pulleyConfig.bowdenBossDia, ...
    1000*predBest.pulleyConfig.bowdenRunClearance)
fprintf('infeasible frames (both routes)     = %d of %d\n', ...
    predBest.infeasibleCount, ctx.N)
fprintf('slack frames (both routes)          = %d of %d\n', ...
    predBest.slackCount, ctx.N)
fprintf('|F_ins| range, routes summed        = %.3f to %.3f N\n', ...
    min(predBest.FinsMag), max(predBest.FinsMag))
fprintf('||ReactionF|| range, routes summed  = %.3f to %.3f N\n', ...
    min(predBest.ReactionMag), max(predBest.ReactionMag))


fprintf('\nBPA and tendon lengths:\n')
fprintf('number of parallel BPAs             = %d\n', predBest.BPAcount)
fprintf('BPA radius source                    = %s\n', predBest.bpaRadiusMode)
fprintf('BPA radius range                     = %.6f to %.6f m\n', ...
    min(predBest.bpaRadius), max(predBest.bpaRadius))
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
fprintf('route-route clearance      = %.6f m (min separation demanded %.6f m)\n', ...
    -cCollision(8), collisionInfo.routeClearanceRequired)
fprintf('min route-route separation = %.6f m (worst over the sweep)\n', ...
    collisionInfo.minRouteSeparation)
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


%% Flexor torque (insertion torque carries the tackle's force loss)
figure('Name', 'Flexor Torque', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predOriginal.TorqueZ, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(ax, ctx.phiD, predBest.TorqueInsZ, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA (insertion)')
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
title(ax, sprintf('Flexor Torque After the Pulley (G = %g)', ...
    predBest.pulleyGain), 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Route vs insertion torque: the transmission trade
figure('Name', 'Insertion Torque', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predBest.TorqueZ, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Route-bookkeeping torque')
plot(ax, ctx.phiD, predBest.TorqueInsZ, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Insertion torque (after tackle)')
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
title(ax, 'Flexor Torque After the Pulley', 'FontName', fontName, ...
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
plot(ax, ctx.phiD, predBest.strain_p, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'strain\_p (pulley class)')
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


%% Pulley transmission trade: insertion force and mount reaction
figure('Name', 'Pulley Forces', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, predBest.FinsMag, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', '|F\_ins|, insertion force')
plot(ax, ctx.phiD, predBest.ReactionMag, '-.', ...
    'Color', c{6}, 'LineWidth', originalLineWidth, ...
    'DisplayName', '||ReactionF||, mount reaction')
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, 'Force, N', 'FontName', fontName, ...
    'FontSize', axesFontSize, 'FontWeight', 'bold')
title(ax, sprintf(['Pulley Insertion and Reaction Forces ' ...
    '(nBPA = %d, G = %g, %s)'], predBest.pulleyConfig.nPulleyBPA, ...
    predBest.pulleyGain, predBest.pulleyConfig.routingMode), ...
    'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
lg = legend(ax, 'Location', 'best');
set(lg, 'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'Box', 'off')
grid(ax, 'off')


%% Torque margin fraction
% Positive means the BPA exceeds the required human flexor-torque
% magnitude. Negative means a remaining torque shortfall. The flexor
% torque curves themselves remain signed and negative; absolute values are
% used only here to form the magnitude ratio. The margin uses the
% INSERTION torque so the tackle's force loss is visible.
humanAbsAtRobotAngles = interp1( ...
    ctx.humanAngleD, ctx.humanTorqueAbs, ctx.phiD, 'pchip', 'extrap');
validHumanTorque = humanAbsAtRobotAngles > 100*eps;
insertionMarginFraction = nan(size(predBest.TorqueInsZ));
insertionMarginFraction(validHumanTorque) = ...
    abs(predBest.TorqueInsZ(validHumanTorque)) ./ ...
    humanAbsAtRobotAngles(validHumanTorque) - 1;


fprintf('\nTorque margin relative to human target (insertion torque):\n')
fprintf('minimum margin fraction = %+.6f (%+.2f%%)\n', ...
    min(insertionMarginFraction(validHumanTorque)), ...
    100*min(insertionMarginFraction(validHumanTorque)))
fprintf('mean remaining shortfall = %.6f (%.2f%%)\n', ...
    mean(max(0, -insertionMarginFraction(validHumanTorque))), ...
    100*mean(max(0, -insertionMarginFraction(validHumanTorque))))


figure('Name', 'Torque Margin Fraction', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', figurePosition)
ax = gca;
hold(ax, 'on')
plot(ax, ctx.phiD, 100*insertionMarginFraction, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Insertion margin (after tackle)')
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


%% ============================================================
%  Local functions: the pulley evaluation pipeline
%  ============================================================
function J = objective_KneeFlexorPulley(x, ctx)
% Objective for the reverse-pulley spin-off. Structurally objective_
% KneeFlexor20mm with one deliberate swap: the torque requirement is
% computed on the INSERTION torque Torque_ins (after the tackle's 1/G
% force loss), so a design can no longer meet the human target on
% BPA-side force alone. With the pulley disabled (ctx.optimizePulleyGain
% false and tackleLineParts == 1) the problem keeps Opt_run's 8-variable
% structure (bounds, constraint count). Torque_ins equals the route
% torque at G = 1 on the sanity/smoke geometry only; on the real flexor
% route it is close but not identical (the moving wrap segment sits
% inside MuscleLength but outside the pulley span -- see the note at the
% ctx.pulley* knobs).


pred = predictKneeFlexorPulley(x, ctx);

if ~pred.ok
    J = 1e12;
    return
end

if any(~isfinite(pred.TorqueIns(:))) || ...
        any(~isfinite(pred.strain_p)) || ...
        any(~isfinite(pred.Contraction))
    J = 1e11;
    return
end

humanAbs = interp1( ...
    ctx.humanAngleD, ...
    ctx.humanTorqueAbs, ...
    ctx.phiD, ...
    'pchip', ...
    'extrap');

requiredTorque = (1 + ctx.requiredTorqueMargin)*humanAbs(:);
robotAbs = abs(pred.TorqueInsZ(:));

% Pointwise torque requirement with worst-position priority.
torqueShortfall = max(0, requiredTorque - robotAbs);
shortfallFraction = torqueShortfall ./ ctx.torqueScale;
Jworst = max(shortfallFraction);
Jtorque = mean(shortfallFraction.^2);

% Very small shape term; meeting the requirement takes precedence.
shapeError = (robotAbs - requiredTorque) ./ ctx.torqueScale;
Jshape = mean(shapeError.^2);

% Enforce the off-axis requirement separately for x and y.  At each
% knee position, penalize only the amount by which the new component
% magnitude exceeds the corresponding original no-wrap component.
offAxisExcessX = max(0, ...
    abs(pred.TorqueX(:)) - abs(ctx.originalTorqueX(:)));
offAxisExcessY = max(0, ...
    abs(pred.TorqueY(:)) - abs(ctx.originalTorqueY(:)));

normalizedExcessX = offAxisExcessX ./ ctx.offAxisTorqueScaleX;
normalizedExcessY = offAxisExcessY ./ ctx.offAxisTorqueScaleY;

JoffAxis = ctx.offAxisPenaltyWeight * ( ...
    max(normalizedExcessX).^2 + mean(normalizedExcessX.^2) + ...
    max(normalizedExcessY).^2 + mean(normalizedExcessY.^2));

% Resting length constraint:
% rest + tendon + 2*fitting >= distance(p1, w2)
% pred.cRestLength <= 0 is feasible.
JrestLength = 1e6 * max(0, pred.cRestLength).^2;

% Strain feasibility. The pulley class carries no Xi3, so strain_f
% (used by the base objective for force) equals strain_p here.
maxStrainAllowed = ctx.maxRelStrain * pred.KMAX;

JstrainHi = 1e7 * max(0, ...
    max([pred.strain_f; pred.Contraction]) - maxStrainAllowed).^2;
JstrainLo = 1e7 * max(0, ...
    -min(pred.Contraction)/pred.KMAX).^2;

% The new contact must release during the extension-to-flexion sweep.
% The constant term makes a no-release route unacceptable.
if ~any(pred.routeInfo.active) || pred.routeInfo.releaseFound
    JwrapRelease = 0;
else
    remainingTurn = max(0, -pred.routeInfo.finalSignedTurnD)/180;
    JwrapRelease = ctx.wrapReleasePenaltyWeight*(1 + remainingTurn).^2;
end

% Keep solution near practical geometry unless torque requires otherwise.
x0 = ctx.x0(:);
dx = x(:) - x0;

geomScale = [ ...
0.030;   % p1 x
0.300;   % p1 y -- large movement allowed
0.030;   % p1 z
0.030;   % p2 x
0.030;   % p2 y
0.030];  % p2 z

Jgeom = 1e-2 * sum((dx(1:6)./geomScale).^2);
Jlen  = 1e-3 * ((x(7) - x0(7))/0.040).^2;

% Soft BPA-2 anchor (Ben, 2026-09-21): prefer keeping the second
% BPA's DERIVED endpoints near their targets -- the Q-angle / aLDFA
% build line of the loaded design.  Score term only, never a
% constraint.  Weight 10 with a 1 cm norm: 1 cm drift on both ends
% costs 20, while 1 N m of worst-angle torque shortfall costs ~4300,
% so torque always wins when they conflict.
if isfield(ctx, 'bpa2TargetP1') && ~isempty(ctx.bpa2TargetP1) && ...
        isfield(pred, 'p1B')
    anchorScale1 = norm(pred.p1B(:) - ctx.bpa2TargetP1(:))/ctx.bpa2AnchorNorm;
    anchorScale2 = norm(pred.p2B(:) - ctx.bpa2TargetP2(:))/ctx.bpa2AnchorNorm;
    Jbpa2Anchor = ctx.bpa2AnchorWeight*(anchorScale1^2 + anchorScale2^2);
else
    Jbpa2Anchor = 0;
end

J = 1e5*Jworst + 1e3*Jtorque + 1e-2*Jshape + ...
    JoffAxis + JrestLength + JstrainHi + JstrainLo + ...
    JwrapRelease + Jgeom + Jlen + Jbpa2Anchor;

if ~isfinite(J)
    J = 1e12;
end
end


function pred = predictKneeFlexorPulley(x, ctx)
%PREDICTKNEEFLEXORPULLEY Evaluate one flexor design over the full angle
% grid with MonoPam_pulley objects. x supplies the eight optimizer
% variables ([p1(1:3), pEnd(1:3), rest, tendon]) plus the pulley gain as
% the ninth variable when ctx.optimizePulleyGain is true. The two routes
% are built by buildKneeFlexorRoute20mm exactly as predictKneeFlexor20mm
% does; each route is wrapped in a MonoPam_pulley with the run's
% transmission config (nPulleyBPA, tackleLineParts, routing mode, exit
% row, Bowden envelope). Route torques add across the pair; the pulley
% outputs are reported per route and summed as scalars.

    p1       = x(1:3);       % parent-frame / femur-side attachment
    p2       = x(4:6);       % theta1-frame end/insertion design variable
    rest     = x(7);         % active BPA rest length
    tendon   = x(8);         % physical tendon length

    % The travel gain: the ninth variable when it is optimized, else the
    % run's tackle line parts (tackleLineParts = G, 1 = pulley disabled).
    if isfield(ctx, 'optimizePulleyGain') && ctx.optimizePulleyGain && numel(x) > 8
        G = x(9);
    else
        G = ctx.tackleLineParts;
    end

    % Fixed identified stiffness parameters (Xi3/BendMeasure are NOT part
    % of the MonoPam_pulley contract).
    Xi0 = ctx.Xi0;
    Xi1 = ctx.Xi1;
    Xi2 = ctx.Xi2;

    BPAcount = ctx.BPAcount;

    KMAX = ctx.KMAX;
    kmax = rest*(1-KMAX);      % measured free-contracted BPA length at 620 kPa

    pred.ok = true;
    pred.failReason = "";

    if rest <= 0 || kmax <= 0 || kmax >= rest || tendon < 0
        pred.ok = false;
        pred.failReason = "Invalid length parameters";
        return
    end
    if ~(isscalar(G) && isfinite(G) && G >= 1)
        pred.ok = false;
        pred.failReason = "Pulley gain must be finite and >= 1";
        return
    end

    % Run-level transmission config for the MonoPam_pulley constructor.
    % Shared-tackle aggregation (Ben's rig: ONE block-and-tackle fed by
    % the parallel BPAs): with the two-route pipeline (BPAcount = 2) the
    % physical parallel BPAs ARE the two routes, so each route object
    % carries ONE BPA line (nPulleyBPA = 1) and the route SUM reproduces
    % the rig's total pull F_t = (F1 + F2)/G = nPulleyBPA*F_single/G.
    % Folding ctx.pulleyBPACount into EVERY route object would double-
    % count (2 routes x 2 BPAs = 4). The class's nPulleyBPA > 1 path is
    % used only when the design is modeled as ONE corridor (BPAcount =
    % 1), where the identical parallel BPAs fold into the transmission
    % (biPulleySpecsFromOpenSim symmetric-bundle convention).
    if BPAcount == 1
        perRouteN = ctx.pulleyBPACount;
    else
        perRouteN = 1;
        if ctx.pulleyBPACount ~= BPAcount
            warning(['predictKneeFlexorPulley:ConfigMismatch', ...
                'ctx.pulleyBPACount = %d but the pipeline models %d ' ...
                'route corridors; the route SUM carries the parallel-BPA ' ...
                'pull (total = %d BPAs feeding one tackle).'], ...
                ctx.pulleyBPACount, BPAcount, BPAcount)
        end
    end
    pulleyConfig = struct( ...
        'nPulleyBPA', perRouteN, ...
        'tackleLineParts', ctx.tackleLineParts, ...
        'gain', G, ...
        'routingMode', ctx.pulleyRoutingMode, ...
        'pulleyExitIndex', ctx.pulleyExitIndex, ...
        'bowdenBossDia', ctx.bowdenBossDia, ...
        'bowdenRunClearance', ctx.bowdenRunClearance);
    pred.pulleyConfig = pulleyConfig;
    pred.pulleyGain = G;

    try
        ctxUsed = ctx;
        % routeInfo.pWrapT1 is the moving t1-frame wrap-point array.
        [Location1, ~, routeInfo1] = ...
            buildKneeFlexorRoute20mm(p1, p2, tendon, ctxUsed);

        if BPAcount == 2
            [p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
            [Location2, ~, routeInfo2] = ...
                buildKneeFlexorRoute20mm(p1B, p2B, tendon, ctxUsed);
            bpa1 = MonoPam_pulley( ...
                ctx.Name, Location1, ctx.CrossPoint, ctx.Dia, ctx.T_Pam, ...
                rest, kmax, tendon, ctx.fitting, ctx.targetPressure, ...
                Xi0, Xi1, Xi2, ctx.wraps, pulleyConfig);
            bpa2 = MonoPam_pulley( ...
                ctx.Name, Location2, ctx.CrossPoint, ctx.Dia, ctx.T_Pam, ...
                rest, kmax, tendon, ctx.fitting, ctx.targetPressure, ...
                Xi0, Xi1, Xi2, ctx.wraps, pulleyConfig);
        elseif BPAcount == 1
            bpa1 = MonoPam_pulley( ...
                ctx.Name, Location1, ctx.CrossPoint, ctx.Dia, ctx.T_Pam, ...
                rest, kmax, tendon, ctx.fitting, ctx.targetPressure, ...
                Xi0, Xi1, Xi2, ctx.wraps, pulleyConfig);
            bpa2 = [];
            routeInfo2 = [];
            p1B = p1;
            p2B = p2;
            Location2 = [];
        else
            error('predictKneeFlexorPulley:BPAcount', ...
                'BPAcount must be 1 or 2 for the current flexor model.')
        end
    catch ME
        pred.ok = false;
        pred.failReason = string(ME.message);
        return
    end

    pred.bpa1 = bpa1;
    if BPAcount == 2
        pred.bpa2 = bpa2;
    end

    % Torques add across the pair (NaN rows propagate: infeasible frames
    % stay nonfinite, which the objective prices).
    pred.Torque   = bpa1.Torque_p + bpa2.Torque_p;
    pred.TorqueX  = bpa1.Torque_p(:,1) + bpa2.Torque_p(:,1);
    pred.TorqueY  = bpa1.Torque_p(:,2) + bpa2.Torque_p(:,2);
    pred.TorqueZ  = bpa1.Torque_p(:,3) + bpa2.Torque_p(:,3);
    pred.offAxisTorque = hypot(pred.TorqueX, pred.TorqueY);

    % Pulley-honest outputs: the objective's torque series is the
    % INSERTION torque (after the tackle's 1/G force loss).
    pred.TorqueIns  = bpa1.Torque_ins + bpa2.Torque_ins;
    pred.TorqueInsZ = pred.TorqueIns(:,3);
    pred.FinsMag = vecnorm(bpa1.F_ins, 2, 2) + vecnorm(bpa2.F_ins, 2, 2);
    pred.ReactionMag = bpa1.ReactionFmag + bpa2.ReactionFmag;
    pred.infeasibleCount = nnz(bpa1.PulleyInfeasible | bpa2.PulleyInfeasible);
    pred.slackCount = nnz(bpa1.PulleySlack | bpa2.PulleySlack);

    % Route-based bookkeeping of BPA 1 (constraint/display continuity with
    % predictKneeFlexor20mm; the route strains exclude Xi3 by design).
    pred.Location = Location1;   % BPA 1 route; nonlcon checks both routes
    pred.LocationAll = {Location1};
    if BPAcount == 2
        pred.LocationAll = {Location1; Location2};
        pred.Location2 = Location2;
        pred.routeInfo2 = routeInfo2;
        pred.p1B = p1B;
        pred.p2B = p2B;
    end
    pred.routeInfo = routeInfo1;
    pred.geo = ctx.geo;
    pred.bpaRadiusMode = string(ctx.geo.bpaRadiusMode);
    pred.bpaRadius = routeInfo1.bpaRs;
    pred.bpaRadiusIteration = 0;
    pred.bpaRadiusConverged = true;
    pred.bpaRadiusChange = 0;

    pred.Contraction = bpa1.strain_p(:);
    pred.strain_p = bpa1.strain_p(:);
    pred.strain_f = bpa1.strain_p(:);   % no Xi3 term in the pulley class
    pred.strain = pred.strain_p;
    pred.forceMismatch = zeros(size(pred.Contraction));
    pred.momentArmVectorAll = bpa1.mA_p;
    pred.momentArmVector = bpa1.mA_p;
    pred.pathLength0 = bpa1.MuscleLength(:);
    pred.pathLength = bpa1.Lmt_p(:) + Xi0;
    pred.delta_L = bpa1.deltaL(:);
    pred.gama = bpa1.gama(:);
    pred.sContraction = bpa1.sContraction(:);
    pred.PulleyTravel = bpa1.PulleyTravel(:);
    pred.Ftendon = bpa1.Ftendon(:);
    pred.relativeContraction = pred.Contraction ./ KMAX;
    pred.relativeStrainF = pred.strain_f ./ KMAX;
    pred.relativeStrainP = pred.strain_p ./ KMAX;
    pred.relativeStrain = pred.relativeStrainP;
    pred.activeLength = rest .* (1 - pred.strain_p);
    pred.momentArm = hypot(pred.momentArmVector(:,1), pred.momentArmVector(:,2));

    % Store design variables for optimizer output.
    pred.p1 = p1;
    pred.p2 = p2;
    pred.pEnd = p2;
    pred.pWrap = routeInfo1.pWrapT1(ctx.idxExtension,:);
    pred.rest = rest;
    pred.tendon = tendon;
    pred.KMAX = KMAX;
    pred.kmax = kmax;
    pred.BPAcount = BPAcount;

    % Extension-frame geometry check.
    idx = ctx.idxExtension;

    % v2 is the final insertion point in the knee/ICR frame.
    v2 = Location1(3,:,idx);

    % w2 is that same point transformed into the femur frame.
    w2 = RowVecTrans(ctx.T_Pam(:,:,idx), v2);

    pred.v2 = v2;
    pred.w2 = w2;

    pred.extensionDistance = pred.pathLength0(idx);

    % Diagnostic only: this is the modeled zero-strain musculotendon length.
    % Positive Xi0 makes the model behave as if the required Lmt is longer.
    pred.restLmt = rest + tendon + 2*ctx.fitting + Xi0;

    % Constraint value <= 0 is feasible:
    % distance(p1, v2 at extension) <= rest + tendon + 2*fitting
    pred.cRestLength = pred.extensionDistance - pred.restLmt;

end


function out = objconstrExclusionPulley(x, obj, geo, ctx, idxP2)
% Objective/constraint wrapper for surrogateopt (pulley spin-off).
%
% surrogateopt wants:
%   out.Fval
%   out.Ineq <= 0

out.Fval = obj(x);

[c, ~] = nonlconExclusionPulley(x, geo, ctx, idxP2);

out.Ineq = c;

end


function [c, ceq, info] = nonlconExclusionPulley(x, geo, ctx, idxP2)
%NONLCONEXCLUSIONPULLEY The 7 base exclusion constraints (nonlconExclusion,
% called unmodified) PLUS the spin-off's route-to-route clearance
% constraint as c(8): the moving via sweeps an envelope as the joint
% rotates and the nPulleyBPA parallel BPAs each occupy a body diameter, so
% the two routes' segment geometry (including the swept via segments) must
% stay at least the sum of the two BPA collision radii plus
% ctx.routeClearanceMargin apart, worst over the orientation sweep --
% same worst-over-sweep style as the existing exclusion constraints.
% ceq is returned empty: patternsearch's nonlcon contract is [c, ceq].

ceq = [];

[cBase, ~, info] = nonlconExclusion(x, geo, ctx, idxP2);

% "" is a 1x1 string (not isempty); test the text length instead.
if strlength(string(info.failReason)) > 0 || any(~isfinite(cBase))
    c = ones(8,1);
    info.failReason = string(info.failReason);
    info.routeClearanceRequired = NaN;
    info.minRouteSeparation = NaN;
    return
end

[cRoute, infoRoute] = routeClearanceConstraint(x, geo, ctx, idxP2);

info.routeClearanceRequired = infoRoute.requiredWorst;
info.minRouteSeparation = infoRoute.minSeparation;

c = [cBase(:); cRoute];

end


function [cRoute, infoRoute] = routeClearanceConstraint(x, geo, ctx, idxP2)
%ROUTECLEARANCECONSTRAINT Worst-over-sweep route-to-route clearance
% violation (c <= 0 feasible): min over orientations of the distance
% between the two routes' segment geometry (both polylines evaluated at
% the SAME orientation -- the two BPAs physically exist simultaneously;
% the min over orientations covers the swept via envelope) must stay at
% least the sum of the two BPA collision radii + margin.
%
% Note: routes 1 and 2 are kept DISTINCT even when nPulleyBPA > 1 -- the
% nPulleyBPA multiplier folds the identical parallel BPAs into the
% transmission math (biPulleySpecsFromOpenSim places them symmetrically
% about the OpenSim line), while the two modeled routes are the two
% physical corridors that must not collide.

p1 = x(1:3);
p2 = x(idxP2);
tendon = x(8);

try
    ctxUsed = ctx;
    ctxUsed.geo = geo;
    [Location1, ~, ~] = buildKneeFlexorRoute20mm(p1, p2, tendon, ctxUsed);

    if ctx.BPAcount ~= 2
        cRoute = -1;    % no second corridor: constraint trivially satisfied
        infoRoute = struct('requiredWorst', NaN, 'minSeparation', inf);
        return
    end

    [p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
    [Location2, ~, ~] = buildKneeFlexorRoute20mm(p1B, p2B, tendon, ctxUsed);

    % Per-frame BPA collision radius (scalar mode: the hand-set geo value;
    % may still be an N-array).
    bpaRs = geo.bpaRs;
    if isscalar(bpaRs)
        bpaRs = repmat(bpaRs, ctx.N, 1);
    end
    bpaRs = bpaRs(:);
    required = 2*bpaRs + ctx.routeClearanceMargin;

    worst = -inf;
    dMinWorst = inf;
    for ii = 1:ctx.N
        PA = routeInFemurFrame(Location1(:,:,ii), ctx.T_Pam(:,:,ii), ctx.CrossPoint);
        PB = routeInFemurFrame(Location2(:,:,ii), ctx.T_Pam(:,:,ii), ctx.CrossPoint);
        dMin = inf;
        for ia = 1:size(PA,1)-1
            for ib = 1:size(PB,1)-1
                dMin = min(dMin, segSegDistance( ...
                    PA(ia,:), PA(ia+1,:), PB(ib,:), PB(ib+1,:)));
            end
        end
        % Degenerate polylines (all segments zero-length) collapse to the
        % point-to-point distance.
        if ~isfinite(dMin)
            dMin = min(vecnorm(PA - PB, 2, 2));
        end
        worst = max(worst, required(ii) - dMin);
        dMinWorst = min(dMinWorst, dMin);
    end
    cRoute = worst;
    infoRoute = struct('requiredWorst', max(required), ...
        'minSeparation', dMinWorst);
catch
    cRoute = 1;    % infeasible: the sweep could not be evaluated
    infoRoute = struct('requiredWorst', NaN, 'minSeparation', NaN);
end

end


function P = routeInFemurFrame(Li, Tii, crossPoint)
% One route page (Mx3) with the distal-body rows transformed into the
% femur frame so both routes live in ONE frame for the distance test.
P = Li;
for j = crossPoint:size(P,1)
    P(j,:) = RowVecTrans(Tii, P(j,:));
end
end


function dist = segSegDistance(p1q, p2q, q1, q2)
% Closest distance between 3-D segments [p1q,p2q] and [q1,q2]
% (standard clamped closest-point construction). Zero-length segments
% return inf so the caller can fall back to point distances.
d1 = p2q - p1q;
d2 = q2 - q1;
r  = p1q - q1;
a = dot(d1, d1);
e = dot(d2, d2);
if a <= 1e-14 || e <= 1e-14
    dist = inf;
    return
end
f = dot(d2, r);
c = dot(d1, r);
b = dot(d1, d2);
denom = a*e - b*b;
if denom > 1e-14
    s = min(1, max(0, (b*f - c*e)/denom));
else
    s = 0;
end
t = (b*s + f)/e;
if t < 0
    t = 0;
    s = min(1, max(0, -c/a));
elseif t > 1
    t = 1;
    s = min(1, max(0, (b - c)/a));
end
dist = norm((p1q + s*d1) - (q1 + t*d2));
end


function smokeSanityCheck()
%SMOKESANITYCHECK The Opt_sanity_pulley.m construction (92-orientation
% knee geometry, straight-tendon route whose crossing segment is the
% tendon line) run through MonoPam_pulley: (a) the G = 1 regression
% identity vs MonoPamDataExplicit_balance (rel tol 1e-8 on F_p, mA_p,
% Torque_p) and (b) a G = 2 closure spot check (tol 1e-6). Xi uses the
% deterministic fallback: the identity is Xi-independent.

fprintf('SMOKE: building the Opt_sanity_pulley geometry (92 orientations).\n')

positions = 92;
R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
T_ICR_t1 = zeros(4, 4, positions);
c = pi/180;

knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

kneeMin = -45*c;
kneeMax = 10*c;
phi = linspace(kneeMin, kneeMax, positions);
[~, posHome] = min(abs(phi));
phi(posHome) = 0;

for i = 1:positions
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');
    t1toICR = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR');
end

% Straight-tendon route (3 rows: origin, tackle exit, insertion),
% CrossPoint = 3, PulleyExitIndex = 2, rows 1-2 constant.
p1 = [-0.050, 0.390, 0.050];
pExit = [-0.040, 0.360, 0.048];
p2 = [0.0574, 0.0355, 0.005];
Location = zeros(3, 3, positions);
routeLength = zeros(positions, 1);
for i = 1:positions
    p2ICR = RowVecTrans(T_ICR_t1(:,:,i), p2);
    Location(:,:,i) = [p1; pExit; p2ICR];
    routeLength(i) = norm(p1 - pExit) + ...
        norm(pExit - RowVecTrans(T_Pam(:,:,i), p2ICR));
end

Name = 'Bicep Femoris (Short Head)';
Dia = 20;
tendon0 = 0.015;
fitting = 0.021;
pres = 620;
wraps = 6;
KMAX = 0.255;
Xi0 = 0.00394; Xi1 = 3.998e4; Xi2 = 1.473e4;
try
    S = load('minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat', ...
        'filtered_results', 'xCols');
    g = S.filtered_results(77, S.xCols);
    Xi0 = g(1); Xi1 = g(2); Xi2 = g(3);
    clear S g
catch
    fprintf('SMOKE: stiffness mat not found; fallback Xi (identity is Xi-independent).\n')
end

rest0 = max(routeLength) - Xi0 - tendon0 - 2*fitting;
kmax0 = (1 - KMAX)*rest0;
cfgG1 = struct('nPulleyBPA', 1, 'tackleLineParts', 1, 'pulleyExitIndex', 2);
cfgG2 = struct('nPulleyBPA', 1, 'tackleLineParts', 2, 'pulleyExitIndex', 2);

base_chk = MonoPamDataExplicit_balance(Name, Location, 3, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps);
pul_chk = MonoPam_pulley(Name, Location, 3, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfgG1);

for fld = ["F_p", "mA_p", "Torque_p"]
    relErr = norm(pul_chk.(fld)(:) - base_chk.(fld)(:))/ ...
        max(norm(base_chk.(fld)(:)), eps);
    fprintf('SMOKE identity check: %-8s relative difference = %.3e\n', ...
        fld, relErr)
    if ~(relErr < 1e-8)
        error(['SMOKE identity check FAILED: G = 1 %s relative ' ...
            'difference %.3e exceeds 1e-8.'], fld, relErr)
    end
end

pul2_chk = MonoPam_pulley(Name, Location, 3, Dia, ...
    T_Pam, rest0, kmax0, tendon0, fitting, pres, Xi0, Xi1, Xi2, wraps, cfgG2);
if any(pul2_chk.PulleyInfeasible) || any(pul2_chk.PulleySlack)
    error('SMOKE closure check FAILED: the G = 2 sweep must be feasible and taut.')
end
Fbpa_chk = festo4(Dia, pul2_chk.sContraction(:)/rest0/KMAX, pres) ...
    .* pul2_chk.Fmax;
resFt = max(abs(pul2_chk.kSpr .* pul2_chk.gama(:) - Fbpa_chk/2));
resClosure = max(abs((2 .* pul2_chk.PulleyTravel(:) - pul2_chk.deltaL(:)) ...
    - pul2_chk.gama(:)));
fprintf(['SMOKE closure check (G = 2): max |kSpr*gama - F_BPA/G| = %.3e N, ' ...
    'max |(2*PulleyTravel - DeltaL) - gama| = %.3e m\n'], resFt, resClosure)
if resFt >= 1e-6 || resClosure >= 1e-6
    error('SMOKE closure check FAILED: residuals %.3e N / %.3e m exceed 1e-6.', ...
        resFt, resClosure)
end
fprintf('SMOKE sanity-geometry check PASSED.\n')
end
