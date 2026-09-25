% Opt Run Bi Pulley
% Author: Ben Bolen
% Date: 2026-09-24
% Description: BATCH driver for the biarticular reverse-pulley actuator
% class (BiPam_pulley). Loops a manifest of muscle specs -- built by
% biPulleySpecsFromMasters.m (master's attachment data,
% "D:\Bipedal humanoid\Matlab code") or biPulleySpecsFromOpenSim.m
% (gait2392_robotbody.osim) -- builds each two-hinge transform grid, sizes
% the BPA from the route, and optimizes [BPA resting length; tendon length]
% to maximize the worst-grid torque margin about the ACTIVE crossing, with
% the strain window and pulley feasibility as nonlinear constraints.
%
% The transform grid uses the BiPamData storage convention: dim 3 = angle
% index, dim 4 = JOINT index -- T(:,:,i1,1) is joint 1's transform at angle
% i1 and T(:,:,i2,2) is joint 2's at angle i2 (two stacked slabs; the
% class chains them per crossing in its frame-mapping helpers).
%
% MANIFEST FIELDS (struct array; the builders fill these):
%   name / source  'semimem_r' + 'masters' | 'med_gas_r' + 'opensim' ...
%   points         Npts x 3 route rows grouped by frame (BiPamData rule)
%   cross          1x2 crossing row indices
%   diameter       default 20 mm
%   pulley         transmission config struct (nPulleyBPA, tackleLineParts,
%                  routingMode; defaults 1 / 1 / 'moving_via'; scalar
%                  fields apply to both crossings) -- passed through as
%                  BiPam_pulley's 15th argument
%   bundleOffsets  nPulleyBPA x 3 symmetric bundle placement about the
%                  original line (empty when nPulleyBPA = 1; OpenSim
%                  builder only -- the master's routes are single-line)
%   bundleMargin   m between adjacent BPA skins (spacing verification)
%   footprintWidth m total bundle footprint tangent to the anatomy
%                  ((n-1)*spacing + diameter; ECHOED IN THE SUMMARY CSV)
%   fmaxSource     char documenting the Fmax origin; fmax = N
%   frames         1x3 cell of frame labels
%   kin            hinge data: pivots (2x3), axes (2x3), types (1x2 cell),
%                  thetaRangesDeg (2x2), kneeSlot (0/1/2), optional kneeFit
%                  (masters rolling knee) / transTheta-transX-transY (osim
%                  rolling knee), pointFolds (osim subtalar fold; applied
%                  to slot-3 rows here in the driver)
%   notes          char provenance string
%
% INTERFERENCE (Ben Q, 2026-09-24): the batch's routes come from the
% manifest (fixed), so route-to-route clearance cannot be SEARCHED here
% (that belongs to the Opt_run exclusion family); it is CHECKED and
% reported instead: min over orientation pairs of the distance between any
% two routes' segment geometry (each route swept over its whole
% orientation grid, which includes the moving-via envelope) must stay >=
% the sum of the two BPA radii + margin; each muscle is checked against
% the routes of the muscles already executed in the same batch, and the
% bundle's own parallel members are spacing-verified. Violations are
% reported in the console and the summary CSV; they do not change the
% [rest, tendon] objective (which cannot fix them).
%
% RUN MODES:
%   RUN_BATCH = false (default): DRY-RUN -- resolves the manifest, prints
%     the plan, runs nothing, writes nothing. (Console is NOT cleared so
%     batch logs survive; deviation from Opt_run.m noted on purpose.)
%   RUN_BATCH = true + FULL_RUN = true: real budgets (surrogateopt 1000 +
%     patternsearch 5000, parallel) -- EASTEREGG2 ONLY (parpool(10)).
%   RUN_BATCH = true + FULL_RUN = false: check-scale budgets, serial.
%   env OPT_BIPULLEY_SMOKE=1: one muscle, tiny budget, UseParallel false,
%     completes in minutes.
%
% RESULT ARTIFACTS (per executed muscle, Results\):
%   BiPulley_<source>_<name>_<stamp>.mat  dated FULL-WORKSPACE bare save
%     behind a per-iteration liveRun flag (Ben directive: a variable-list
%     save or an exist()-guard is a BUG; liveRun is cleared immediately
%     before the save so rerunning a section from a loaded mat cannot mint
%     a new dated mat).
%   BiPulley_batch_summary.csv  one appended row per muscle: name, source,
%     exit flags, objective, worst constraint, torque margin, gain,
%     footprint width, min route clearance, stamp.
%   BiPulley_checkpoint.mat     done-muscle checkpoint so an interrupted
%     batch resumes; delete it to re-run everything.

% DRY-RUN by default. Set true to execute the batch.
RUN_BATCH = false;

% Real optimizer budgets belong on easteregg2 (10 cores, 128 GB RAM).
FULL_RUN = false;

smoke = strcmp(getenv('OPT_BIPULLEY_SMOKE'), '1');
% OPT_BIPULLEY_SMOKE=1 FORCES the one-muscle tiny-budget run regardless of
% RUN_BATCH (the env flag is the batch-contract's smoke switch).
if smoke
    RUN_BATCH = true;
end

%% Path setup (repo pattern: Mesh_Optimization wins shadowing, Festo appended)
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root) || isempty(parent)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

resDir = fullfile(scriptDir, 'Results');
if ~exist(resDir, 'dir')
    mkdir(resDir)
end

%% Resolve the manifest (builders extract + document the sources)
manifest = defaultManifest();

N1 = 5;     % joint-1 angles (deg range spec.kin.thetaRangesDeg(1,:))
N2 = 5;     % joint-2 angles

% Shared stiffness defaults (Xi in the hand-tune neighborhood; per AGENTS
% notes Xi1 sits on a flat likelihood valley, so keep it explicit here).
baseXi0 = 0.005;        % m
baseXi1 = 5e5;          % N/m
baseXi2 = 1e4;          % N/m
baseWraps = 6;          % 20 mm cable-wrap multiplier (Spr)
basePressure = 620;     % kPa
tauArmNominal = 0.05;   % m, nominal moment arm for the torque target

%% Print the resolved plan (DRY-RUN output; also the batch preamble)
fprintf('\n================ BI PULLEY BATCH PLAN ================\n')
fprintf('Mode                = %s\n', modeString(RUN_BATCH, FULL_RUN, smoke))
fprintf('Orientation grid    = %d x %d angles per muscle\n', N1, N2)
fprintf('Stiffness defaults  = Xi0 %g m, Xi1 %g N/m, Xi2 %g N/m, wraps %d\n', ...
    baseXi0, baseXi1, baseXi2, baseWraps)
fprintf('Torque target       = fmax * %.0f cm nominal arm\n', 100 * tauArmNominal)
fprintf('Artifacts           = Results\\BiPulley_batch_summary.csv + ')
fprintf('dated mats + checkpoint\n')
fprintf('Muscles in manifest = %d\n', numel(manifest))
for k = 1:numel(manifest)
    printPlanItem(manifest(k), N1, N2)
end
fprintf('======================================================\n')

if ~RUN_BATCH
    fprintf(['DRY-RUN complete: the plan above is fully resolved. ' ...
        'Set RUN_BATCH = true to execute (full budgets belong on ' ...
        'easteregg2), or OPT_BIPULLEY_SMOKE=1 for a minutes-long smoke.\n'])
    return
end

%% Batch execution
if smoke
    fprintf(['SMOKE MODE: executing muscle 1 only with a tiny budget and ' ...
        'UseParallel = false.\n'])
    manifest = manifest(1);
    useParallel = false;
    budgetS = 0;       % surrogateopt stage skipped in smoke
    budgetP = 60;      % patternsearch evaluations
else
    useParallel = FULL_RUN;   % parallel only on easteregg2 (FULL_RUN)
    if FULL_RUN
        budgetS = 1000;    % surrogateopt (Opt_run.m scale)
        budgetP = 5000;    % patternsearch refinement
    else
        budgetS = 30;      % check scale, not a real run
        budgetP = 100;
        fprintf(['NOTE: FULL_RUN = false -- these are check-scale budgets. ' ...
            'Real runs belong on easteregg2 with FULL_RUN = true.\n'])
    end
end

% Checkpoint (resume support)
ckptFile = fullfile(resDir, 'BiPulley_checkpoint.mat');
doneKeys = {};
if isfile(ckptFile)
    Sck = load(ckptFile, 'doneKeys');
    doneKeys = Sck.doneKeys;
    fprintf('Checkpoint found: %d muscle(s) already done.\n', numel(doneKeys))
end

summaryFile = fullfile(resDir, 'BiPulley_batch_summary.csv');
if ~isfile(summaryFile)
    summaryFid = fopen(summaryFile, 'a');
    fprintf(summaryFid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'name', ...
        'source', 'exitflagS', 'exitflagP', 'objective', 'worstConstraint', ...
        'torqueMargin', 'pulleyGain', 'footprintWidth_mm', ...
        'minRouteClearance_mm', 'timestamp');
    fclose(summaryFid);
end

% Routes already executed in THIS batch (for the route-to-route clearance
% check; see the INTERFERENCE note in the header).
doneRoutes = {};

for k = 1:numel(manifest)
    spec = manifest(k);
    key = sprintf('%s_%s', spec.source, spec.name);
    if any(strcmp(doneKeys, key))
        fprintf('[%d/%d] %s (%s): already done (checkpoint), skipping.\n', ...
            k, numel(manifest), spec.name, spec.source)
        continue
    end
    fprintf('\n[%d/%d] ===== %s (%s) =====\n', k, numel(manifest), ...
        spec.name, spec.source)
    try
        runMuscle(spec, N1, N2, baseXi0, baseXi1, baseXi2, baseWraps, ...
            basePressure, tauArmNominal, useParallel, budgetS, budgetP, ...
            resDir, summaryFile, doneRoutes, smoke);
        doneRoutes{end + 1} = struct('name', spec.name, ...
            'points', pointsWithFold(spec, N1, N2), ...
            'cross', spec.cross);   %#ok<AGROW>
        doneKeys{end + 1} = key;   %#ok<AGROW>
        save(ckptFile, 'doneKeys');
        fprintf('Checkpoint updated (%d done).\n', numel(doneKeys))
    catch err
        fprintf('MUSCLE FAILED (%s): %s\n', spec.name, err.message)
        fprintf('%s\n', getReport(err, 'basic'))
        summaryFid = fopen(summaryFile, 'a');
        fprintf(summaryFid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
            spec.name, spec.source, 'error', 'error', 'nan', 'nan', ...
            'nan', gainString(spec.pulley.nPulleyBPA), ...
            sprintf('%.3f', 1000 * spec.footprintWidth), 'nan', ...
            char(string(datetime('now'), 'yyyyMMdd_HHmm:ss')));
        fclose(summaryFid);
    end
end

fprintf('\nBatch complete: %d/%d muscles done. See %s\n', ...
    numel(doneKeys), numel(manifest), summaryFile)

%% =====================================================================
%% Local functions
%% =====================================================================
function mode = modeString(runBatch, fullRun, smoke) %#ok<INUSD>
if ~runBatch
    mode = 'DRY-RUN (plan only)';
elseif smoke
    mode = 'SMOKE (1 muscle, tiny budget, serial)';
elseif fullRun
    mode = 'FULL RUN (easteregg2 budgets, parallel)';
else
    mode = 'CHECK RUN (small budgets, serial)';
end
end

function out = clean(txt)
out = strtrim(regexprep(txt, '\(.*$', ''));
end

function gs = gainString(g)
if nargin < 1 || isempty(g)
    gs = '1';
else
    gs = strjoin(arrayfun(@(v) sprintf('%g', v), g, 'UniformOutput', false), '/');
end
end

function manifest = defaultManifest()
% Masters routes first, then the OpenSim set. med_gas_r exists in BOTH
% sources (different provenance, both kept); the checkpoint keys on
% source_name so they cannot shadow each other.
manifest = [];
try
    manifest = [manifest, biPulleySpecsFromMasters()];   %#ok<AGROW>
catch err
    fprintf('MASTER''S BUILDER UNAVAILABLE: %s\n', err.message)
end
try
    manifest = [manifest, biPulleySpecsFromOpenSim()];   %#ok<AGROW>
catch err
    fprintf('OSIM BUILDER UNAVAILABLE: %s\n', err.message)
end
if isempty(manifest)
    error('Opt_run_BiPulley:noManifest', ...
        'Neither builder produced specs; check the data sources.')
end
end

function printPlanItem(spec, N1, N2)
fprintf(['  %-11s %-7s pts=%d cross=[%d %d] frames=[%s %s %s] Fmax=%5g N ' ...
    'gain=%s mode=%s\n'], spec.name, spec.source, size(spec.points, 1), ...
    spec.cross(1), spec.cross(2), clean(spec.frames{1}), ...
    clean(spec.frames{2}), clean(spec.frames{3}), spec.fmax, ...
    gainString(spec.pulley.nPulleyBPA), char(spec.pulley.routingMode));
fprintf('    tackleLineParts=%s  footprint=%.1f mm  bundle offsets=%d\n', ...
    gainString(spec.pulley.tackleLineParts), 1000 * spec.footprintWidth, ...
    size(spec.bundleOffsets, 1));
fprintf('    joints: [%s; %s]  ranges(deg): [%g %g; %g %g]\n', ...
    spec.kin.types{1}, spec.kin.types{2}, ...
    spec.kin.thetaRangesDeg(1, :), spec.kin.thetaRangesDeg(2, :));
fprintf('    grid %dx%d angles; fmaxSource: %s\n', N1, N2, spec.fmaxSource);
fprintf('    notes: %s\n', spec.notes);
end

function points = pointsWithFold(spec, N1, N2) %#ok<INUSD>
% The route points the driver actually feeds the class (with the osim
% subtalar fold applied), for the route-to-route clearance archive.
points = spec.points;
if isfield(spec.kin, 'pointFolds') && ...
        any(any(abs(spec.kin.pointFolds{2}) > 0))
    points(spec.cross(2):end, :) = ...
        points(spec.cross(2):end, :) + spec.kin.pointFolds{2}';
end
end

function runMuscle(spec, N1, N2, baseXi0, baseXi1, baseXi2, baseWraps, ...
    basePressure, tauArmNominal, useParallel, budgetS, budgetP, ...
    resDir, summaryFile, doneRoutes, isSmoke)
% One muscle end to end: transforms, BPA sizing, optimization, artifacts.

% --- Two-hinge transform grid (BiPamData slabs: dim3 = angle, dim4 = joint)
[T, theta1, theta2] = buildBiTransforms(spec, N1, N2);

% Apply the osim subtalar fold to slot-3 rows (documented simplification).
points = pointsWithFold(spec, N1, N2);

% --- Characteristic route length for sizing/bounds --------------------
Lrig = rigidRouteLengths(points, spec.cross, T, N1, N2);
Lchar = max(Lrig(:));

% --- Interference checks (fixed routes: check + report, see header) ----
fprintf(['Clearance checks (BPA diameter %d mm, margin %.1f mm):\n'], ...
    spec.diameter, 1000 * spec.bundleMargin)
if size(spec.bundleOffsets, 1) > 1
    nOff = size(spec.bundleOffsets, 1);
    minSelf = inf;
    for a = 1:nOff - 1
        for b = a + 1:nOff
            minSelf = min(minSelf, norm(spec.bundleOffsets(b, :) - ...
                spec.bundleOffsets(a, :)));
        end
    end
    needSelf = spec.diameter * 1e-3 + spec.bundleMargin;
    fprintf(['  bundle self-clearance: min member spacing %.1f mm ' ...
        '(need >= %.1f mm): %s\n'], 1000 * minSelf, 1000 * needSelf, ...
        verdictString(minSelf >= needSelf))
end
minPair = inf;
pairName = 'none (first muscle)';
for d = 1:numel(doneRoutes)
    c = routePairClearance(points, spec.cross, doneRoutes{d}.points, ...
        doneRoutes{d}.cross, T, N1, N2);
    if c < minPair
        minPair = c;
        pairName = doneRoutes{d}.name;
    end
end
if ~isempty(doneRoutes)
    needPair = 2 * (spec.diameter * 1e-3 / 2) + spec.bundleMargin;
    fprintf(['  route-to-route clearance vs executed muscles: min %.1f mm ' ...
        '(worst pair vs %s; need >= %.1f mm): %s\n'], 1000 * minPair, ...
        pairName, 1000 * needPair, verdictString(minPair >= needPair))
end

% --- Objective/constraint over x = [rest; tendon] ---------------------
% Torque target: fmax times a nominal 5 cm arm (documented heuristic;
% per-joint human targets are the Opt_run family's job, not the batch's).
tauTarget = spec.fmax * tauArmNominal;

lb = [0.30 * Lchar; 0.02];
ub = [0.90 * Lchar; max(0.30 * Lchar, 0.06)];

    function [state, ok] = evaluate(xq)
    rest = xq(1);
    tendon = xq(2);
    kmax = 0.75 * rest;      % 25% maximum contraction fraction (KMAX)
    obj = BiPam_pulley(spec.name, points, spec.cross, 20, T, ...
        rest, kmax, tendon, 0.025, basePressure, ...
        baseXi0, baseXi1, baseXi2, baseWraps, spec.pulley);
    aK = obj.ActiveCrossing;
    tauAbs = abs(obj.Torque_ins(:, :, :, aK));
    validTau = tauAbs(~isnan(tauAbs));
    if isempty(validTau)
        state = struct('margin', -1e3, 'worstC', 1e3, 'obj', [], ...
            'tauAbs', tauAbs);
        ok = false;
        return
    end
    margin = min(validTau) / tauTarget - 1;
    % Constraints (c <= 0 feasible): strain floor per grid cell, and zero
    % pulley-infeasible cells.
    strainGrid = reshape(obj.strain_p, size(Lrig));
    c1 = obj.strainFloor() - min(strainGrid(:));
    c2 = double(sum(obj.PulleyInfeasible(:)));
    state = struct('margin', margin, 'worstC', max([c1, c2, 0]), ...
        'obj', obj, 'tauAbs', tauAbs);
    ok = true;
    end

    function f = batchObjective(xq)
    [st, ok] = evaluate(xq);
    if ~ok
        f = 1e3;
    else
        f = -st.margin + 10 * max(0, st.worstC);
    end
    end

    function [c, ceq] = batchNonlcon(xq)
    % patternsearch calls nonlcon with two outputs ([c, ceq]).
    [st, ok] = evaluate(xq);
    if ~ok
        c = 1e3;
    else
        c = st.worstC;
    end
    ceq = [];
    end

    function [fVal, c, ceq] = objconstrWrap(xq)
    % surrogateopt objconstr signature: [f, c, ceq].
    fVal = batchObjective(xq);
    c = batchNonlcon(xq);
    ceq = [];
    end

x0 = [0.60 * Lchar; 0.05];
f0 = batchObjective(x0);
fprintf('Initial: x = [%.4f %.4f] m, penalized margin = %.4f\n', x0, f0)

% --- Optimization ------------------------------------------------------
if useParallel && isempty(gcp('nocreate'))
    parpool;     % easteregg2 only (FULL_RUN)
end

x = x0;
f = f0;
exitflagS = NaN;
exitflagP = NaN;

% Stage 1: surrogateopt global pass (skipped in smoke so a smoke completes
% in minutes with patternsearch alone).
if budgetS > 0
    optsS = optimoptions('surrogateopt', ...
        'Display', 'iter', ...
        'UseParallel', useParallel, ...
        'MaxFunctionEvaluations', budgetS, ...
        'ConstraintTolerance', 1e-6);
    [xS, fS, exitflagS] = surrogateopt(@objconstrWrap, lb, ub, optsS);
    if fS < f
        x = xS;
        f = fS;
    end
end

% Stage 2: patternsearch refinement.
optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', useParallel, ...
    'MaxFunctionEvaluations', budgetP, ...
    'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4, ...
    'ConstraintTolerance', 1e-6);
[xBest, fBest, exitflagP] = patternsearch( ...
    @batchObjective, x, [], [], [], [], lb, ub, @batchNonlcon, optsP);

% --- Final evaluation and reporting ------------------------------------
[st, ok] = evaluate(xBest);
if ~ok
    error('Opt_run_BiPulley:finalEval', ...
        'Final evaluation failed for %s', spec.name)
end
finalObj = st.obj;
fprintf('\nBest: rest = %.4f m, tendon = %.4f m, margin = %+.4f, ', ...
    xBest, st.margin)
fprintf('worst constraint = %+.6g\n', st.worstC)
fprintf('Exit flags: surrogateopt %g, patternsearch %g\n', exitflagS, exitflagP)
fprintf('Active crossing = %d, infeasible cells = %d, slack cells = %d\n', ...
    finalObj.ActiveCrossing, sum(finalObj.PulleyInfeasible(:)), ...
    sum(finalObj.PulleySlack(:)))

% Reverse-pulley transmission + packaging echoes (Ben's 2026-09-24 notes).
fprintf('Transmission: nPulleyBPA = %s, tackleLineParts = %s, gain = %s, mode = %s\n', ...
    gainString(finalObj.NPulleyBPA), gainString(finalObj.TackleLineParts), ...
    gainString(finalObj.PulleyGain), strjoin(finalObj.RoutingMode, '/'))
if any(strcmp(finalObj.RoutingMode, 'bowden'))
    fprintf(['Bowden envelope (Shimano-type): boss dia %.1f mm, run ' ...
        'clearance %.1f mm -- verify the housing path against both\n'], ...
        1000 * finalObj.BowdenBossDia, 1000 * finalObj.BowdenRunClearance)
end
fprintf('Bundle footprint width = %.1f mm (diameter %d mm)\n', ...
    1000 * spec.footprintWidth, spec.diameter)

% Dated FULL-WORKSPACE result capture (Ben directive: bare save behind
% liveRun; cleared just before saving so section reruns from a loaded mat
% cannot mint a new dated mat). SMOKE never saves (liveRun stays
% undefined), matching the Opt_run_pulley convention.
liveRun = ~isSmoke;
if exist('liveRun', 'var') && liveRun
    stamp = char(string(datetime('now'), 'yyyyMMdd_HHmm'));
    resultFile = fullfile(resDir, ...
        sprintf('BiPulley_%s_%s_%s.mat', spec.source, spec.name, stamp));
    clear liveRun
    save(resultFile)
    fprintf('Saved %s\n', resultFile)
end

% Summary CSV row (open/append/close: fprintf needs a file ID, not a path).
summaryFid = fopen(summaryFile, 'a');
fprintf(summaryFid, '%s,%s,%g,%g,%.6g,%.6g,%.6g,%s,%.3f,%.3f,%s\n', ...
    spec.name, spec.source, exitflagS, exitflagP, fBest, st.worstC, ...
    st.margin, gainString(finalObj.PulleyGain), ...
    1000 * spec.footprintWidth, 1000 * minPair, ...
    char(string(datetime('now'), 'yyyyMMdd_HHmm:ss')));
fclose(summaryFid);
end

function s = verdictString(ok)
if ok
    s = 'OK';
else
    s = 'VIOLATION';
end
end

function c = routePairClearance(pointsA, crossA, pointsB, crossB, T, N1, N2)
% Min over orientation pairs (and segment pairs) of the distance between
% two routes' segment geometry, each swept over its whole grid -- the
% moving-via envelope is included by taking the min across cells.
segsA = routeSegments(pointsA, crossA, T, N1, N2);
segsB = routeSegments(pointsB, crossB, T, N1, N2);
c = inf;
for a = 1:size(segsA, 1)
    for b = 1:size(segsB, 1)
        c = min(c, segSegDist(segsA(a, 1:3), segsA(a, 4:6), ...
            segsB(b, 1:3), segsB(b, 4:6)));
    end
end
end

function S = routeSegments(points, cross, T, N1, N2)
% All route segments at every grid cell, expressed in the proximal frame:
% rows = cell*segment (column-major over the grid), [pA pB] flattened.
nPts = size(points, 1);
S = zeros(N1 * N2 * (nPts - 1), 6);
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for iii = 1:N2
        T2 = T(:, :, iii, 2);
        for r = 1:nPts - 1
            row = (iii - 1) * N1 * (nPts - 1) + (ii - 1) * (nPts - 1) + r;
            S(row, 1:3) = mapRow(points, cross, r, T1, T2, 1);
            S(row, 4:6) = mapRow(points, cross, r + 1, T1, T2, 1);
        end
    end
end
end

function d = segSegDist(p1, q1, p2, q2)
% Minimum distance between two 3-D segments (clamped closest points,
% standard Ericson "Real-Time Collision Detection" formulation).
d1 = q1 - p1;
d2 = q2 - p2;
r = p1 - p2;
a = dot(d1, d1);
e = dot(d2, d2);
f = dot(d2, r);
c = dot(d1, r);
b = dot(d1, d2);
if a <= 1e-18 && e <= 1e-18
    d = norm(r);
    return
elseif a <= 1e-18
    s = 0;
    t = min(1, max(0, f / e));
elseif e <= 1e-18
    t = 0;
    s = min(1, max(0, -c / a));
else
    denom = a * e - b * b;
    if denom > 1e-18
        s = min(1, max(0, (b * f - c * e) / denom));
    else
        s = 0;
    end
    t = (b * s + f) / e;
    if t < 0
        t = 0;
        s = min(1, max(0, -c / a));
    elseif t > 1
        t = 1;
        s = min(1, max(0, (b - c) / a));
    end
end
d = norm((p1 + s * d1) - (p2 + t * d2));
end

function [T, theta1, theta2] = buildBiTransforms(spec, N1, N2)
% Two-hinge transform grid in the BiPamData storage convention:
%   T(:,:,i1,1) = joint 1 at angle i1 (maps middle -> proximal frame)
%   T(:,:,i2,2) = joint 2 at angle i2 (maps distal -> middle frame)
% Translations are the joint pivots (child-frame origin = joint), so the
% class's moment arms are about the joints. The rolling knee (masters
% kneeFit / osim transTheta-transX-transY) moves the KNEE pivot with its
% angle; kneeSlot says which joint the knee is (0 = none).
kin = spec.kin;
r1 = rangeOrDefault(kin, 1, [-30, 90]);
r2 = rangeOrDefault(kin, 2, [-110, 20]);
theta1 = linspace(r1(1), r1(2), N1) * pi / 180;
theta2 = linspace(r2(1), r2(2), N2) * pi / 180;

% Rolling-knee data (either source's convention flattens to theta/x/y).
thK = [];
xK = [];
yK = [];
if isfield(kin, 'kneeSlot') && ~isempty(kin.kneeSlot) && kin.kneeSlot > 0
    if isfield(kin, 'kneeFit') && ~isempty(kin.kneeFit) && ...
            ~isempty(kin.kneeFit.thetaY)
        thK = kin.kneeFit.thetaY;
        xK = interp1(kin.kneeFit.thetaX, kin.kneeFit.x, thK, ...
            'linear', 'extrap');
        yK = kin.kneeFit.y;
    elseif isfield(kin, 'transTheta') && ~isempty(kin.transTheta)
        thK = kin.transTheta;
        xK = interp1(thK, kin.transX, thK, 'linear', 'extrap');
        yK = kin.transY;
    end
end

T = zeros(4, 4, N1, N2);
for i1 = 1:N1
    p1 = kin.pivots(1, :);
    if kin.kneeSlot == 1 && ~isempty(thK)
        p1 = [interp1(thK, xK, theta1(i1), 'linear', 'extrap'), ...
            interp1(thK, yK, theta1(i1), 'linear', 'extrap'), p1(3)];
    end
    T(:, :, i1, 1) = hingeTrans(kin.types{1}, kin.axes(1, :), p1, theta1(i1));
end
for i2 = 1:N2
    p2 = kin.pivots(2, :);
    if kin.kneeSlot == 2 && ~isempty(thK)
        p2 = [interp1(thK, xK, theta2(i2), 'linear', 'extrap'), ...
            interp1(thK, yK, theta2(i2), 'linear', 'extrap'), p2(3)];
    end
    T(:, :, i2, 2) = hingeTrans(kin.types{2}, kin.axes(2, :), p2, theta2(i2));
end
end

function T = hingeTrans(jtype, axis, pivot, theta)
% RpToTrans-built hinge; 'fixed' returns identity (child origin still at
% the pivot point so moment arms stay about the joint).
switch jtype
    case 'fixed'
        R = eye(3);
    otherwise
        R = axisRot(axis, theta);
end
T = RpToTrans(R, pivot(:));
end

function R = axisRot(axis, theta)
% Rodrigues rotation about a unit axis.
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end

function r = rangeOrDefault(kin, slot, def)
if isfield(kin, 'thetaRangesDeg') && size(kin.thetaRangesDeg, 1) >= slot
    r = kin.thetaRangesDeg(slot, :);
else
    r = def;
end
end

function L = rigidRouteLengths(points, cross, T, N1, N2)
% Rigid route length per grid cell (driver-side sizing helper; mirrors the
% class's segment chaining without constructing the object).
nPts = size(points, 1);
L = zeros(N1, N2);
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for iii = 1:N2
        T2 = T(:, :, iii, 2);
        tot = 0;
        for r = 1:nPts - 1
            pA = mapRow(points, cross, r, T1, T2, 1);
            pB = mapRow(points, cross, r + 1, T1, T2, 1);
            tot = tot + norm(pA - pB);
        end
        L(ii, iii) = tot;
    end
end
end

function p = mapRow(points, cross, r, T1, T2, frame)
% Driver-side twin of the class's frame chaining for one row.
if r < cross(1)
    f = 1;
elseif cross(1) == cross(2)
    f = 3;
elseif r < cross(2)
    f = 2;
else
    f = 3;
end
p = points(r, :);
while f > frame
    if f == 3
        p = RowVecTrans(T2, p);
        f = 2;
    else
        p = RowVecTrans(T1, p);
        f = 1;
    end
end
end
