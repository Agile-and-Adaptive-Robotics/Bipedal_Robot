% Opt Run Bi Pulley survey 20261022
% Author: Ben Bolen (machinery of record); survey driver assembled
% 2026-10-22 by the BiPulley survey engineer per the small-budget
% per-muscle capability-survey ask (workflow, 2026-10-22).
%
% Description: NEW DATED DRIVER (copy-and-mark-edits convention; no Ben
% file edited). Adapts the 2026-10-22 BiPulley machinery --
%   tmp/paper_bridge/opt_improve_20261022/Opt_run_BiPulley2_20261022.m
%   (config injection, escalation-ladder loop, runMuscle skeleton),
%   Code/Matlab/Mesh_Optimization/BiPam_pulley2_20261022.m (variant
%   class; called with its DEFAULT config so the executed arithmetic is
%   bit-identical to BiPam_pulley -- no mount/bowden compliance fields),
%   Code/Matlab/Mesh_Optimization/biPulleySpecsFromOpenSim.m (spec
%   builder, one name per call so a failing seed never kills the batch),
%   tmp/paper_bridge/opt_improve_20261022/escalation_advisor_20261022.m
%   (fuzzy ladder: torque-short -> +1 BPA, RoM/slack-short -> +1 line
%   part; taut-cell-only torque test),
%   hipAdductionEnvelope_20261022.m / softCornerWeights_20261022.m /
%   taperKernel_20261022.m (envelope + guard machinery, reused as the
%   cellWeights hook so a ladder step can never be justified by an
%   out-of-envelope corner cell),
%   Code/Matlab/Mesh_Optimization/PULLEY_SWEEP_DIRECTIVE_20261021.md
%   (sizes re-solved per configuration; force gain = nBPA/G, travel
%   gain = G; ladder instead of blind sweep),
%   Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json
%   (the 27 actuators, robot_seed_preference order, rom_table_deg).
%
% SURVEY CONTRACT (the ask, 2026-10-22):
%   1. SITES: attempt ALL 27 actuators; feed each actuator's
%      robot_seed_preference names to biPulleySpecsFromOpenSim IN ORDER
%      until one resolves; muscles that cannot build are recorded and
%      skipped, never fatal. A resolved spec whose PRIMARY dof is not
%      modeled by the builder's hinge chain (e.g. subtalar-primary
%      routes that the builder reduces to the ankle alias, hip-adduction
%      primary vs the builder's flexion-first hip hinge, or the
%      torso-start psoas route whose crossing-1 key 'torso_to_femur_r'
%      has no hinge) is recorded grade-unable and skipped the same way.
%   2. CONFIG: free per muscle, nBPA {1,2,3} x tackle line parts
%      G {1,2,3}, driven by escalation_advisor_20261022; sizes
%      [rest; tendon] RE-SOLVED per rung per the directive.
%   3. OBJECTIVE (soft, nothing rejected): minimize
%         -(1000*meetFrac + activeFrac)
%      meetFrac  = graded-RoM cells where the BPA torque about the
%                  primary-dof axis has the same sign as the human
%                  torque AND |tau_BPA| >= |tau_human| (the human curve
%                  from the CSV; cells with |tau_human| <= 1e-9 count
%                  as met when |tau_BPA| >= |tau_human|, sign vacuous);
%      activeFrac= graded-RoM cells with positive strain AND taut
%                  torque > 0. POSITIVE STRAIN = the class's strain_p
%                  > 0 (contraction-positive convention: the effective
%                  BPA length sits below RestingL, i.e. the route
%                  demand exceeds the laid-out BPA -- the ask's "length
%                  > resting length" phrasing mapped onto the class's
%                  sign convention; the literal elongation reading is
%                  unsatisfiable because the closure clamps s >= 0).
%                  Taut torque > 0 = ~PulleySlack & ~PulleyInfeasible
%                  & |tau_BPA| > 0.
%      NO hard performance constraints, NO nonlcon: pulley-infeasible
%      cells (out-of-envelope moving-via point travel) are priced
%      softly as not-active/not-meeting cells and the best design is
%      kept -- exactly the ask's "soft and nothing rejected".
%   4. BUDGET CAPS: <= 300 optimizer evaluations per muscle TOTAL
%      across all rungs (default rung split 60 x at most 5 rungs:
%      (1,1) -> at most 4 monotone +1 steps -> (3,3)), ~10 minutes wall
%      per muscle (checked between rungs), parpool <= 8. This is a
%      capability survey, not a design campaign.
%   5. BOUNDS: rest length in [0.30*Lchar, 0.90*Lchar] (the listed
%      resting-length ceiling is an UPPER limit -- the solver may
%      shorten, never lengthen; 0.30*Lchar is the class's existing
%      lower convention), tendon in [0.02, max(0.30*Lchar, 0.06)]
%      (Opt_run_BiPulley.m bounds, unchanged). Route points are NOT
%      design variables: the builder's robotbody route is used verbatim
%      (origin/insertion envelope respected trivially); the 10-22
%      envelope/guard machinery enters as the advisor's cellWeights.
%   6. PRESSURE: 620 kPa everywhere.
%   7. GRADING: human_torque_curves_20261022.csv (written concurrently
%      by the OpenSim stage) is the single source of truth. Every
%      candidate is graded on exactly that 25-point primary-dof grid;
%      the secondary joint is pinned at its default (0 deg) per the
%      builder convention. The driver refuses to run without the CSV.
%   8. OUTPUTS: dated per-muscle mats in Code/Matlab/Mesh_Optimization/
%      Results/ named BiPulley_survey_<site>_nB<n>_G<g>_<stamp>.mat
%      (new names, never overwritten -- stamp yyyyMMdd_HHmmss; each
%      (site, nB, G) pair is visited at most once per ladder, and a new
%      run gets a new stamp), plus driver-level survey + ladder CSVs in
%      tmp/paper_bridge/capability_table_20261022/.
%
% The deterministic replay gate is NOT re-run here: this workflow's
% verify phase already passed it for this machinery (ask, 2026-10-22);
% Opt_run_BiPulley2_20261022.m mode 'gate' remains the gate of record.
%
% FULL-LAUNCH EDIT (2026-10-22, launch ask): (a) resume checkpoint --
% the full-run survey CSV carries the FIXED name
% BiPulley_survey_full_results_20261022.csv, one row appended per
% actuator at its terminal status; a restarted driver reads it back and
% skips finished actuators (per-muscle caps are per-launch: an
% interrupted muscle restarts from rung 1 with a fresh budget); smoke
% runs keep stamped throwaway names. (b) site-level convergence
% bookkeeping: stopReason in {converged, maxrungs, eval_budget,
% wall_budget} + convergedSite recorded per row (the LAST rung's
% advisor verdict, not the best rung's) so the harvester can flag
% cap-hit-unconverged muscles. (c) slot-2 graded-grid layout fixed per
% the class's real T contract (joint-2 transforms indexed along dim 3;
% see buildGradedTransforms) -- the first slot-2 probe (med_gas) had
% crashed on it. (d) tie-break-scale dead-cell price in the OPTIMIZER
% objective only (see surveyObj); reported scores stay pure.
% (e) 5-point box multi-start per rung (see runRung) so the flat
% all-dead plateau cannot strand the search where a taut window exists.
% Everything else unchanged from the smoke-verified driver (smoke
% 2026-10-22 13:34: tib_ant_r 96.0% meet / 80.0% active at (2,1);
% deterministic mat re-grade reproduced both).
%
% RUN MODES (env):
%   default          DRY-RUN: resolve all 27 sites through the builder
%                    (the buildable / cannot-build lists), print the
%                    plan, exit. No CSV needed, no optimization.
%   OPT_SURVEY_RUN=1   execute the survey (requires the human CSV).
%   OPT_SURVEY_SMOKE=1 one site (tib_ant_r), serial, tiny budget
%                    (15 evals x 2 rungs), end to end incl. grading;
%                    forces RUN.
%   OPT_SURVEY_PARALLEL=1  parpool(min(8, ...)) for the full run.
%   OPT_SURVEY_MAXEVAL   per-muscle total eval cap (default 300).
%   OPT_SURVEY_RUNGEVAL  per-rung eval budget (default 60; smoke 15).
%   OPT_SURVEY_TIMEMIN   per-muscle wall budget, minutes (default 10).
%   OPT_SURVEY_SITES     comma-separated master names to restrict.
%   OPT_SURVEY_HUMANCSV  override CSV path.

smoke = strcmp(getenv('OPT_SURVEY_SMOKE'), '1');
RUN_BATCH = smoke || strcmp(getenv('OPT_SURVEY_RUN'), '1');
useParallel = ~smoke && strcmp(getenv('OPT_SURVEY_PARALLEL'), '1');

maxEvalsSite = iEnvNum('OPT_SURVEY_MAXEVAL', 300);
rungEval = iEnvNum('OPT_SURVEY_RUNGEVAL', 60);
timeBudgetMin = iEnvNum('OPT_SURVEY_TIMEMIN', 10);
maxRungs = 5;
if smoke
    rungEval = iEnvNum('OPT_SURVEY_RUNGEVAL', 15);
    maxRungs = 2;
end
rungEval = min(rungEval, maxEvalsSite);

%% Path setup -----------------------------------------------------------
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
addpath(fullfile(root, 'tmp', 'paper_bridge', 'opt_improve_20261022')); ...
    % escalation_advisor / softCornerWeights / hipAdductionEnvelope /
    % taperKernel live there
addpath(scriptDir);

resultsRoot = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
if ~exist(resultsRoot, 'dir')
    mkdir(resultsRoot)
end
outDir = scriptDir;   % driver-level CSVs live next to this driver

%% Shared physical constants (Opt_run_BiPulley.m defaults, verbatim) ----
baseXi0 = 0.005;          % m
baseXi1 = 5e5;            % N/m
baseXi2 = 1e4;            % N/m
baseWraps = 6;            % 20 mm cable-wrap multiplier (Spr)
basePressure = 620;       % kPa (ask item 6: 620 kPa everywhere)

%% Actuator map (the 27) ------------------------------------------------
mapFile = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', ...
    'Routing27_v2', 'actuator_map_27.json');
if ~isfile(mapFile)
    error('Opt_run_BiPulley_survey_20261022:noMap', ...
        'actuator_map_27.json not found at %s', mapFile)
end
txt = fileread(mapFile);
amap = jsondecode(txt);
actuators = amap.actuators;
nAct = numel(actuators);

siteFilter = strsplit(strtrim(getenv('OPT_SURVEY_SITES')), ',');
hasFilter = ~isempty(strtrim(getenv('OPT_SURVEY_SITES')));

%% Human torque curves (single source of truth for grading) -------------
humanCsv = getenv('OPT_SURVEY_HUMANCSV');
if isempty(humanCsv)
    humanCsv = fullfile(scriptDir, 'human_torque_curves_20261022.csv');
end
if RUN_BATCH
    % Poll briefly here (the driver must also work when launched after
    % the CSV landed; the long poll belongs to the operator).
    waited = 0;
    while ~isfile(humanCsv) && waited < 60
        pause(10);
        waited = waited + 10;
    end
    if ~isfile(humanCsv)
        error('Opt_run_BiPulley_survey_20261022:noHumanCsv', ...
            ['Human torque curves not found at %s -- the OpenSim stage ' ...
             'writes them concurrently; poll and relaunch.'], humanCsv)
    end
    curves = loadHumanCurves(humanCsv);
    fprintf('Human curves loaded: %d actuator curve set(s) from %s\n', ...
        numel(curves), humanCsv)
else
    curves = [];
    if isfile(humanCsv)
        fprintf('Human CSV present (dry-run does not grade): %s\n', humanCsv)
    else
        fprintf('Human CSV not yet present (dry-run does not need it): %s\n', ...
            humanCsv)
    end
end

%% Resolve every site through the builder (ask item 1) ------------------
fprintf('\n================ BI PULLEY SURVEY (20261022) PLAN ================\n')
fprintf('Mode          = %s\n', iModeString(RUN_BATCH, smoke, useParallel))
fprintf('Budgets       = %d evals/muscle TOTAL, %d evals/rung, %d rungs max, %.0f min wall\n', ...
    maxEvalsSite, rungEval, maxRungs, timeBudgetMin)
fprintf('Config ladder = seed (1,1), advisor steps +1 BPA / +1 line part, cap {3,3}\n')
fprintf('Pressure      = %d kPa; sizes [rest; tendon] re-solved per rung\n', ...
    basePressure)
fprintf('Sites         = %d actuators (robot_seed_preference order)\n', nAct)

sites = struct([]);
buildable = {};
cannotBuild = {};
gradeUnable = {};
for a = 1:nAct
    A = actuators(a);
    if hasFilter && ~any(strcmp(siteFilter, A.master_name))
        continue
    end
    if smoke && ~any(strcmp(A.robot_seed_preference, 'tib_ant_r'))
        % Smoke site (ask, 2026-10-22): tib_ant_r = actuator 16.
        continue
    end
    fprintf('[%2d] %-34s primary %-18s seeds: %s\n', A.id, A.master_name, ...
        A.primary_dof, strjoin(A.robot_seed_preference, ' > '))
    seedUsed = '';
    seedErrors = cell(numel(A.robot_seed_preference), 1);
    spec0 = [];
    for s = 1:numel(A.robot_seed_preference)
        nm = A.robot_seed_preference{s};
        try
            cand = biPulleySpecsFromOpenSim({nm});
            spec0 = cand(1);
            seedUsed = nm;
            break
        catch err
            seedErrors{s} = sprintf('%s: %s', nm, err.message);
        end
    end
    if isempty(seedUsed)
        cannotBuild{end + 1} = A.master_name;   %#ok<AGROW>
        fprintf('     CANNOT BUILD (%s)\n', strjoin(seedErrors, ' | '))
        sites = iAppendSite(sites, iSiteRow(A, '', spec0, 'cannot_build', ...
            strjoin(seedErrors, ' | ')));
        continue
    end
    fprintf('     seed resolved: %s (%d points, cross [%d %d])\n', ...
        seedUsed, size(spec0.points, 1), spec0.cross(1), spec0.cross(2))

    % Primary-dof slot: which crossing's hinge models this dof?
    [slot, slotReason, slotCoords] = primarySlotForDof(spec0, ...
        A.primary_dof, fullfile(root, 'Solid_Models', 'OpenSim', ...
        'Gait2392_Robotbody', 'gait2392_robotbody.osim'));
    if slot == 0
        gradeUnable{end + 1} = A.master_name;   %#ok<AGROW>
        fprintf('     GRADE-UNABLE: %s\n', slotReason)
        sites = iAppendSite(sites, iSiteRow(A, seedUsed, spec0, ...
            'grade_unable', slotReason));
        continue
    end
    fprintf('     graded at crossing %d (modeled coords: {%s})\n', ...
        slot, strjoin(slotCoords, ', '))
    row = iSiteRow(A, seedUsed, spec0, 'ready', sprintf('slot %d', slot));
    row.slot = slot;
    sites = iAppendSite(sites, row);
    buildable{end + 1} = A.master_name;         %#ok<AGROW>
end
fprintf('====================================================================\n')
fprintf('Build summary: %d buildable, %d cannot-build, %d grade-unable\n', ...
    numel(buildable), numel(cannotBuild), numel(gradeUnable))

if ~RUN_BATCH
    fprintf(['DRY-RUN complete (no optimization; the human CSV is not ' ...
        'needed here). Set OPT_SURVEY_RUN=1 (or OPT_SURVEY_SMOKE=1) to ' ...
        'execute the survey.\n'])
    return
end

%% Parpool (ask item 4: at most 8) --------------------------------------
if useParallel
    want = min(8, feature('numcores'));
    if isempty(gcp('nocreate'))
        parpool(want);
    else
        fprintf('parpool already running (%d workers).\n', ...
            numWorkers(gcp('nocreate')))
    end
end

%% Driver-level CSVs -----------------------------------------------------
% FULL-LAUNCH EDIT (2026-10-22): the full-run CSVs carry FIXED names and
% act as the RESUME CHECKPOINT. One row is appended per actuator the
% moment its site finishes (any terminal status: surveyed / cannot_build
% / grade_unable / no_human_curve); a restarted driver reads the file
% back and skips every actuator already recorded, so a crashed campaign
% resumes without re-running finished muscles. Per-muscle caps are
% per-launch: a muscle interrupted mid-ladder (no row yet) restarts from
% rung 1 with a fresh budget -- the cap bounds one muscle's survey, not
% the campaign. Smoke runs keep their stamped throwaway names.
if smoke
    stamp0 = char(string(datetime('now'), 'yyyyMMdd_HHmmss'));
    surveyCsv = fullfile(outDir, sprintf('BiPulley_survey_results_%s.csv', stamp0));
    ladderCsv = fullfile(outDir, sprintf('BiPulley_survey_ladder_%s.csv', stamp0));
else
    surveyCsv = fullfile(outDir, 'BiPulley_survey_full_results_20261022.csv');
    ladderCsv = fullfile(outDir, 'BiPulley_survey_full_ladder_20261022.csv');
end
if ~isfile(surveyCsv)
    sFid = fopen(surveyCsv, 'w');
    fprintf(sFid, ['actuator_id,master_name,primary_dof,seed,status,reason,' ...
        'rungs,chosen_nB,chosen_G,rest_m,tendon_m,meetPct,activePct,score,' ...
        'evalsUsed,wallMin,converged,stopReason,ladderChain,matFile\n']);
    fclose(sFid);
end
if ~isfile(ladderCsv)
    lFid = fopen(ladderCsv, 'w');
    fprintf(lFid, ['actuator_id,master_name,seed,rung,nBPA,G,evals,rest_m,' ...
        'tendon_m,meetPct,activePct,score,slackCells,infeasibleCells,' ...
        'tautCells,minTautTau_Nm,targetTau_Nm,verdict,stepNBPA,stepG,' ...
        'nextNBPA,nextG,converged,smoke,timestamp\n']);
    fclose(lFid);
end
doneIds = [];
if ~smoke && isfile(surveyCsv)
    Tprev = readtable(surveyCsv, 'VariableNamingRule', 'preserve', ...
        'Delimiter', ',', 'TextType', 'char');
    doneIds = double(Tprev.actuator_id(:));
    fprintf('RESUME: %d actuator(s) already recorded in %s\n', ...
        numel(doneIds), surveyCsv);
end

%% =================== THE SURVEY (per muscle) ==========================
for a = 1:numel(sites)
    S = sites(a);
    if any(doneIds == S.id)
        fprintf('\n[%d] %s already recorded (checkpoint) -- skipping.\n', ...
            S.id, S.master_name)
        continue
    end
    if ~strcmp(S.status, 'ready')
        % recorded and skipped, never fatal (ask item 1)
        iAppendSurveyRow(surveyCsv, S, []);
        continue
    end
    fprintf('\n########## SURVEY %s [%d] (%s, seed %s) ##########\n', ...
        S.master_name, S.id, S.primary_dof, S.seed)
    tSite = tic;
    evalsUsed = 0;

    % Human graded curve for THIS actuator (25-point primary-dof grid)
    HC = iCurvesFor(curves, S.id, S.master_name, S.primary_dof);
    if isempty(HC)
        fprintf('  no human curve for actuator %d / %s -- recorded, skipped.\n', ...
            S.id, S.primary_dof)
        iAppendSurveyRow(surveyCsv, iSetStatus(S, 'no_human_curve', ...
            sprintf('no curve for dof %s in the CSV', S.primary_dof)), []);
        continue
    end
    anglesDeg = HC.anglesDeg(:)';    % 1 x N (the graded grid, exact)
    tauHuman = reshape(HC.tauNm, size(anglesDeg));
    fprintf('  graded grid: %d points, %.1f..%.1f deg; |tau_human| peak %.2f N*m\n', ...
        numel(anglesDeg), anglesDeg(1), anglesDeg(end), max(abs(tauHuman)))

    % Transforms on the graded grid: primary slot swept over the CSV
    % angles, the other slot pinned at its default (0 deg) -- ask item 7.
    [T, N1, N2, slotSwept] = buildGradedTransforms(S.spec, S.slot, anglesDeg);
    % Human curve onto the class grid: the primary pose varies along
    % dim 3 for slot-1 sites and along dim 4 (columns) for slot-2
    % sites; the pinned axis is duplicated (see the layout contract in
    % buildGradedTransforms), so repmat along the pinned axis.
    mGrid = numel(anglesDeg);
    if slotSwept == 1
        tauHuman = repmat(reshape(tauHuman, mGrid, 1), 1, N2);
    else
        tauHuman = repmat(reshape(tauHuman, 1, mGrid), N1, 1);
    end
    points = pointsWithFold(S.spec);
    Lrig = rigidRouteLengths(points, S.spec.cross, T, N1, N2);
    Lchar = max(Lrig(:));

    % Envelope / guard machinery (ask item 5): cellWeights for the
    % advisor; on the pinned-secondary graded grid these evaluate to 1
    % inside the graded RoM (inclusive kernel) and taper only through
    % the hip-coupling mask at hip sites -- a documented, benign reuse.
    W = gradedCellWeights(S, N1, N2, slotSwept, anglesDeg);

    % Advisor scalar target: the peak |human torque| over the graded
    % grid (the ladder is fuzzy; the peak is the conservative reading).
    targetTau = max(abs(tauHuman(:)));

    axisU = S.spec.kin.axes(S.slot, :);
    axisU = axisU / norm(axisU);

    % ----- Ladder (ask items 2-4) --------------------------------------
    nB = 1; G = 1;
    chain = '';
    best = [];
    stopReason = '';   % FULL-LAUNCH EDIT: site-level convergence bookkeeping
    for rung = 1:maxRungs
        remaining = maxEvalsSite - evalsUsed;
        if remaining <= 2
            fprintf('  eval budget exhausted (%d used) -- stopping ladder.\n', ...
                evalsUsed)
            stopReason = 'eval_budget';
            break
        end
        elapsedMin = toc(tSite) / 60;
        if elapsedMin > timeBudgetMin
            fprintf('  wall budget %.1f min exceeded -- stopping ladder.\n', ...
                elapsedMin)
            stopReason = 'wall_budget';
            break
        end
        budget = min(rungEval, remaining - 5);   % -5: room for the box
        % seed evaluations so the per-muscle total never exceeds the cap
        fprintf('\n  ===== RUNG %d | nBPA = %d, G = %d (sizes re-solved, %d evals) =====\n', ...
            rung, nB, G, budget)
        spec = injectConfig(S.spec, nB, G);
        R = runRung(spec, points, T, N1, N2, Lchar, axisU, tauHuman, ...
            S.slot, basePressure, baseXi0, baseXi1, baseXi2, baseWraps, ...
            budget, useParallel, S, nB, G, rung, resultsRoot, smoke, ...
            ladderCsv, targetTau, W);
        evalsUsed = evalsUsed + R.evals;
        chain = sprintf('%s(%d,%d)[%.0f%%/%.0f%%]', chain, nB, G, ...
            R.meetPct, R.activePct);   %#ok<AGROW>
        if isempty(best) || R.score < best.score || ...
                (R.score == best.score && (nB + G) < (best.nB + best.G))
            best = R;   % keep the best design (ask item 3); tie -> simpler
        end
        if R.advice.converged
            fprintf('  ladder converged at rung %d (config (%d,%d)).\n', ...
                rung, nB, G)
            chain = [chain, ' CONV'];   %#ok<AGROW>
            stopReason = 'converged';
            break
        end
        nB = nB + R.advice.stepNBPA;
        G = G + R.advice.stepG;
        chain = [chain, ' -> '];        %#ok<AGROW>
    end
    if isempty(stopReason)
        stopReason = 'maxrungs';   % ladder length cap reached, advisor
        % still wanted a step -- recorded as an unconverged cap hit
    end
    convergedSite = strcmp(stopReason, 'converged');

    wallMin = toc(tSite) / 60;
    fprintf(['  SITE RESULT %s: chosen (%d,%d), rest %.4f m, tendon %.4f m, ' ...
        'meet %.1f%%, active %.1f%%, %d evals, %.1f min (%s)\n'], ...
        S.master_name, best.nB, best.G, best.xBest(1), best.xBest(2), ...
        best.meetPct, best.activePct, evalsUsed, wallMin, stopReason)
    S.best = best;
    S.evalsUsed = evalsUsed;
    S.wallMin = wallMin;
    S.chain = chain;
    S.stopReason = stopReason;
    S.convergedSite = convergedSite;
    iAppendSurveyRow(surveyCsv, S, best);
end

fprintf('\nSurvey complete. Driver CSVs:\n  %s\n  %s\n', surveyCsv, ladderCsv)

%% =====================================================================
%% Local functions
%% =====================================================================

function R = runRung(spec, points, T, N1, N2, Lchar, axisU, tauHuman, ...
    slot, pres, xi0, xi1, xi2, wraps, budget, useParallel, S, nB, G, ...
    rung, resDir, smoke, ladderCsv, targetTau, W)
% One ladder rung end to end: sizes [rest; tendon] re-solved at THIS
% configuration on the graded grid, soft survey objective, advisor
% diagnosis, dated mat. No hard constraints anywhere (ask item 3).
% Eval accounting uses patternsearch's own output.funccount (exact in
% serial AND parallel -- a worker-side nested counter never reaches the
% client), plus the one client-side initial evaluation below.

lb = [0.30 * Lchar; 0.02];
ub = [0.90 * Lchar; max(0.30 * Lchar, 0.06)];
x0 = [0.60 * Lchar; 0.05];

    function [g, ok] = grade(xq)
    % Grading math (ask items 3 and 7) -- used by BOTH the objective and
    % the final evaluation so the score and the reported percentages can
    % never diverge.
        rest = xq(1);
        tendon = xq(2);
        kmax = 0.75 * rest;
        try
            obj = BiPam_pulley2_20261022(S.seed, points, spec.cross, 20, ...
                T, rest, kmax, tendon, 0.025, pres, xi0, xi1, xi2, ...
                wraps, spec.pulley);
        catch err
            g = iBadGrade(size(tauHuman), err.message);
            ok = false;
            return
        end
        T3 = obj.Torque_ins(:, :, :, slot);          % N1 x 3 x N2 x 1
        tauB = reshape(sum(T3 .* reshape(axisU, 1, 3, 1), 2), N1, N2);
        strain = reshape(obj.strain_p, N1, N2);
        slack = logical(obj.PulleySlack);
        infeas = logical(obj.PulleyInfeasible);
        taut = ~slack & ~infeas & isfinite(tauB);
        meet = taut & isfinite(tauHuman) & ...
            (sign(tauB) == sign(tauHuman) | abs(tauHuman) <= 1e-9) & ...
            abs(tauB) >= abs(tauHuman) - 1e-12;
        active = taut & isfinite(strain) & strain > 0 & abs(tauB) > 0;
        g = struct('tauB', tauB, 'strain', strain, 'slack', slack, ...
            'infeas', infeas, 'taut', taut, 'meet', meet, ...
            'active', active, 'obj', obj);
        ok = true;
    end

    function f = surveyObj(xq)
        [g, ok] = grade(xq);
        if ~ok
            f = 1e6;
            return
        end
        meetFrac = sum(g.meet(:)) / numel(g.meet);
        activeFrac = sum(g.active(:)) / numel(g.active);
        % Ask's suggested score -(1000*meetFrac + activeFrac), PLUS a
        % tie-break-scale dead-cell price 1e-3*(slackFrac+infeasFrac)
        % (FULL-LAUNCH EDIT): without it an all-dead design scores
        % exactly 0 -- BETTER than any partially-active design under
        % minimization -- and patternsearch parks on that plateau (the
        % med_gas probe at x0 was 100% infeasible). The price is <=
        % 0.001 total, ~1000x finer than one active cell, so it never
        % trades capability for dead-cell avoidance; the REPORTED score
        % in every mat/CSV row stays the pure suggested score,
        % recomputed at xBest by grade().
        deadFrac = (sum(g.slack(:)) + sum(g.infeas(:))) / numel(g.meet);
        f = -(1000 * meetFrac + activeFrac) + 1e-3 * deadFrac;
    end

opts = optimoptions('patternsearch', ...
    'Display', 'none', ...
    'UseParallel', useParallel, ...
    'MaxFunctionEvaluations', budget, ...
    'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4);
% FULL-LAUNCH EDIT (e): 5-point box multi-start. The all-dead plateau
% (every cell slack or infeasible) is FLAT even under the dead-cell
% tie-break, so a single blind x0 can strand patternsearch there when a
% taut window does exist in bounds (med_gas, by contrast, has provably
% NO taut window at knee = 0: laidOut - mL_min <= 0.25*rest reduces to
% 0.75*rest + tendon <= 0.0358, unsatisfiable with rest >= 0.0345,
% tendon >= 0.02 -- recorded as the honest 0% capability). Seeds span
% the box corners; the best seed starts the search. Counted in the eval
% accounting below.
ubT = ub(2);
seeds = [0.60 * Lchar, 0.05; ...
         0.30 * Lchar, 0.02; ...
         0.90 * Lchar, 0.02; ...
         0.30 * Lchar, ubT; ...
         0.90 * Lchar, ubT];
fSeed = zeros(size(seeds, 1), 1);
for sIdx = 1:size(seeds, 1)
    fSeed(sIdx) = surveyObj(seeds(sIdx, :));
end
[f0, ixSeed] = min(fSeed);
x0 = seeds(ixSeed, :);
fprintf('  rung initial (best of %d box seeds): x = [%.4f %.4f] m, score %.4f\n', ...
    size(seeds, 1), x0, f0)
% Sign-convention diagnostic (builder axis projection vs the OpenSim
% native-sign human curve): printed every rung; a systematic 0% would
% flag an axis-sign mismatch, not a capability result.
[g0, ok0] = grade(x0);
if ok0
    mSign = g0.taut & abs(tauHuman) > 1;
    if any(mSign(:))
        agree = mean(sign(g0.tauB(mSign)) == sign(tauHuman(mSign)));
        fprintf(['  sign check (initial x): %d taut cells with ' ...
            '|tau_h| > 1 N*m, %.0f%% agree with the human sign\n'], ...
            nnz(mSign), 100 * agree)
    end
end
[xBest, fBest, ~, pout] = patternsearch(@surveyObj, x0, [], [], [], [], ...
    lb, ub, [], opts);
evals = pout.funccount + size(seeds, 1);   % + the box seed evaluations

[g, ok] = grade(xBest);
if ~ok
    error('Opt_run_BiPulley_survey_20261022:finalGrade', ...
        'Final grading failed for %s: %s', S.master_name, g.reason)
end
meetFrac = sum(g.meet(:)) / numel(g.meet);
activeFrac = sum(g.active(:)) / numel(g.active);
score = -(1000 * meetFrac + activeFrac);

% Fuzzy ladder diagnosis (ask item 2): the advisor consumes the graded
% grids (taut-cell-only torque test, slack routed to the RoM side) with
% the 10-22 guard cellWeights so corner cells cannot drive a step.
advStruct = struct('TorqueIns', reshape(abs(g.tauB), N1, N2), ...
    'PulleySlack', g.slack, 'PulleyInfeasible', g.infeas);
[advice, rep] = escalation_advisor_20261022(advStruct, targetTau, ...
    struct('curNBPA', nB, 'curG', G, 'verbose', true, ...
    'cellWeights', W));

fprintf(['  rung %d best: rest %.4f m, tendon %.4f m -> meet %.1f%%, ' ...
    'active %.1f%% (score %.2f), %d evals\n'], rung, xBest, ...
    100 * meetFrac, 100 * activeFrac, score, evals)
fprintf('  rung %d advisor: %s -> steps (+%d,+%d) | %s\n', rung, ...
    advice.verdict, advice.stepNBPA, advice.stepG, advice.reason)

% Dated mat (ask item 8): new name, never overwritten.
stamp = char(string(datetime('now'), 'yyyyMMdd_HHmmss'));
matFile = fullfile(resDir, sprintf('BiPulley_survey_%s_nB%d_G%d_%s.mat', ...
    S.seed, nB, G, stamp));
surveyMeta = struct('actuator_id', S.id, 'master_name', S.master_name, ...
    'primary_dof', S.primary_dof, 'seed', S.seed, 'slot', slot, ...
    'nBPA', nB, 'G', G, 'rung', rung, 'smoke', smoke, ...
    'evals', evals, 'meetPct', 100 * meetFrac, ...
    'activePct', 100 * activeFrac, 'score', score, ...
    'targetTauPeak', targetTau, 'created', datetime('now'));
save(matFile, 'surveyMeta', 'spec', 'points', 'T', 'xBest', 'fBest', ...
    'score', 'tauHuman', 'g', 'advice', 'rep', 'W');
fprintf('  saved %s\n', matFile)

% Ladder CSV row
lFid = fopen(ladderCsv, 'a');
fprintf(lFid, ['%d,%s,%s,%d,%d,%d,%d,%.6g,%.6g,%.4f,%.4f,%.4f,%d,%d,%d,' ...
    '%.6g,%.6g,%s,%d,%d,%s,%s,%d,%d,%s\n'], S.id, S.master_name, S.seed, ...
    rung, nB, G, evals, xBest(1), xBest(2), 100 * meetFrac, ...
    100 * activeFrac, score, sum(g.slack(:)), sum(g.infeas(:)), ...
    sum(g.taut(:)), rep.minTautTau, targetTau, advice.verdict, ...
    advice.stepNBPA, advice.stepG, iNumOrDash(advice.nextNBPA), ...
    iNumOrDash(advice.nextG), advice.converged, smoke, ...
    char(string(datetime('now'), 'yyyyMMdd_HHmmss')));
fclose(lFid);

R = struct('xBest', xBest, 'fBest', fBest, 'score', score, ...
    'meetPct', 100 * meetFrac, 'activePct', 100 * activeFrac, ...
    'evals', evals, 'nB', nB, 'G', G, 'advice', advice, 'rep', rep, ...
    'grade', g, 'matFile', matFile);
end

%% ------------ Graded transform grid (ask item 7) ---------------------
function [T, N1, N2, slotSwept] = buildGradedTransforms(spec, slot, anglesDeg)
% Two-hinge transforms on the EXACT human-CSV grid for the primary
% slot; the other slot pinned at its default (0 deg). Rolling-knee
% translation fits carried verbatim from Opt_run_BiPulley.m:585-634.
%
% CLASS LAYOUT CONTRACT (established 2026-10-22 smoke + med_gas probe,
% debug_grade_{tib_ant,med_gas}_20261022.log): joint-1 transforms live
% at T(:,:,ii,1) and joint-2 transforms at T(:,:,iii,2) -- BOTH indexed
% along dim 3, dim 4 selects the joint (slot 2 is the only dim-4 index
% ever read). Hence dim3 (the class's N1) must cover the JOINT-2 angle
% count, and dim4 (N2) IS the joint-2 count. A pinned joint contributes
% a duplicated transform at every index it occupies:
%   slot-1 swept: N1 = 25 primary angles, N2 = 2 (pinned duplicated);
%   slot-2 swept: N1 = N2 = 25 -- the pinned joint-1 transform is
%     duplicated across ALL 25 dim-3 slices, so the 625-cell grid has
%     25 identical rows per primary pose and the graded percentages
%     equal the 25-point curve exactly (each pose counted 25 times).
% (The 5x5 base driver never revealed this because N1 == N2 == 5.)
kin = spec.kin;
thetaP = anglesDeg(:)' * pi / 180;      % primary sweep (exact grid)
nA = numel(thetaP);
slotSwept = slot;
if slot == 1
    N1 = nA;
    N2 = 2;
    th1v = thetaP;
    th2v = [0, 0];
else
    N1 = nA;
    N2 = nA;
    th1v = zeros(1, N1);   % pinned joint 1, duplicated across dim 3
    th2v = thetaP;
end

thK = []; xK = []; yK = [];
if isfield(kin, 'kneeSlot') && ~isempty(kin.kneeSlot) && kin.kneeSlot > 0
    if isfield(kin, 'kneeFit') && ~isempty(kin.kneeFit) && ...
            ~isempty(kin.kneeFit.thetaY)
        thK = kin.kneeFit.thetaY;
        xK = interp1(kin.kneeFit.thetaX, kin.kneeFit.x, thK, 'linear', 'extrap');
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
        p1 = [interp1(thK, xK, th1v(i1), 'linear', 'extrap'), ...
            interp1(thK, yK, th1v(i1), 'linear', 'extrap'), p1(3)];
    end
    T(:, :, i1, 1) = hingeTrans(kin.types{1}, kin.axes(1, :), p1, th1v(i1));
end
for i2 = 1:N2
    p2 = kin.pivots(2, :);
    if kin.kneeSlot == 2 && ~isempty(thK)
        p2 = [interp1(thK, xK, th2v(i2), 'linear', 'extrap'), ...
            interp1(thK, yK, th2v(i2), 'linear', 'extrap'), p2(3)];
    end
    T(:, :, i2, 2) = hingeTrans(kin.types{2}, kin.axes(2, :), p2, th2v(i2));
end
end

%% ------------ Config injection (Opt_run_BiPulley2 injectConfig) ------
function spec = injectConfig(spec0, nB, G)
spec = spec0;
spec.pulley = struct('nPulleyBPA', nB, 'tackleLineParts', [G, 1], ...
    'routingMode', 'moving_via');
spacing = spec.diameter * 1e-3 + 0.005;   % builder's bundleMargin = 5 mm
if nB > 1
    u = [];
    for r = size(spec.points, 1) - 1:-1:1
        cand = spec.points(r + 1, :) - spec.points(r, :);
        if norm(cand) > 1e-12
            u = cand / norm(cand);
            break
        end
    end
    if abs(u(3)) < 0.9
        perp = cross(u, [0, 0, 1]);
    else
        perp = cross(u, [0, 1, 0]);
    end
    perp = perp / norm(perp);
    offs = ((1:nB)' - (nB + 1) / 2) * spacing;
    spec.bundleOffsets = offs .* perp;
    spec.footprintWidth = (nB - 1) * spacing + spec.diameter * 1e-3;
else
    spec.bundleOffsets = zeros(0, 3);
    spec.footprintWidth = spec.diameter * 1e-3;
end
end

%% ------------ Route points with the osim subtalar fold ---------------
function points = pointsWithFold(spec)
points = spec.points;
if isfield(spec.kin, 'pointFolds') && ...
        any(any(abs(spec.kin.pointFolds{2}) > 0))
    points(spec.cross(2):end, :) = ...
        points(spec.cross(2):end, :) + spec.kin.pointFolds{2}';
end
end

%% ------------ Guard cellWeights (ask item 5) -------------------------
function W = gradedCellWeights(S, N1, N2, slotSwept, anglesDeg)
% Reuse of the 10-22 envelope/guard machinery as the advisor's
% cellWeights: soft corner taper over the graded RoM (inclusive kernel:
% the graded grid IS the envelope, so W = 1 on it) plus the hip
% flexion-x-adduction coupling mask at hip sites (min AND semantics).
% With the secondary pinned at its default the coupling mask is 1 across
% the band -- a documented, benign reuse that still fires the moment a
% future 2-D grading grid pairs flexion with adduction.
if slotSwept == 1
    t1 = repmat(anglesDeg(:), 1, N2);       % N1 x N2 (N2 = 1)
    t2 = zeros(N1, N2);
else
    t1 = zeros(N1, N2);
    t2 = repmat(anglesDeg, N1, 1);          % N1 x N2 (N1 = 1)
end
env1 = [min(anglesDeg), max(anglesDeg)];
env2 = [-1, 1];   % pinned joint: single interior point, kernel = 1
opts = struct('band', 10, 'wMin', 0.1);   % PENDING-BEN defaults (10-22)
W = softCornerWeights_20261022(t1, t2, env1, env2, opts);
isHip = strncmp(S.primary_dof, 'hip_', 4);
if isHip
    if slotSwept == 1
        flexDeg = t1;
    else
        flexDeg = t2;
    end
    adducDeg = zeros(size(flexDeg));        % adduction pinned at default
    [hipMask, hinfo] = hipAdductionEnvelope_20261022(flexDeg, adducDeg);
    if ~isempty(hinfo.banner)
        fprintf('  %s\n', hinfo.banner);    % envelope contract: print it
    end
    W = min(W, hipMask);   % AND semantics: the tighter guard wins
end
end

%% ------------ Primary-dof slot matching ------------------------------
function [slot, reason, slotCoords] = primarySlotForDof(spec, dofName, osimFile)
% Which builder hinge (if any) models dofName. Replicates the builder's
% own JointSet reduction (biPulleySpecsFromOpenSim.osimKinematics: each
% CustomJoint -> hinge on its FIRST coordinate) and its ankle alias
% (tibia_r>calcn_r reads the ankle joint; header note 6), then matches
% the dof against the coordinate modeled at each crossing's body pair.
jmap = jointFirstCoords(osimFile);
b = cellfun(@(f) regexprep(strtrim(f), '\s*\(.*$', ''), spec.frames, ...
    'UniformOutput', false);
pairs = {b{1}, b{2}; b{2}, b{3}};
coords = cell(1, 2);
for s = 1:2
    if strcmpi(strtrim(spec.kin.types{s}), 'fixed')
        coords{s} = '';
        continue
    end
    coords{s} = iCoordForPair(jmap, pairs{s, 1}, pairs{s, 2});
end
slotCoords = coords(:);
for s = 1:2
    if strcmp(coords{s}, dofName)
        slot = s;
        reason = sprintf('crossing %d models %s', s, dofName);
        return
    end
end
slot = 0;
reason = sprintf(['primary dof %s not modeled by this spec''s hinge ' ...
    'chain (modeled coords: slot1 "%s", slot2 "%s")'], dofName, ...
    coords{1}, coords{2});
end

function [jmap] = jointFirstCoords(osimFile)
doc = xmlread(osimFile);
model = doc.getElementsByTagName('Model').item(0);
js = model.getElementsByTagName('JointSet').item(0);
jNodes = js.getElementsByTagName('CustomJoint');
jmap = containers.Map('KeyType', 'char', 'ValueType', 'any');
for i = 0:jNodes.getLength() - 1
    j = jNodes.item(i);
    jname = char(j.getAttribute('name'));
    parentOff = iOffsetFrame(j, 'socket_parent_frame');
    childOff = iOffsetFrame(j, 'socket_child_frame');
    if isempty(parentOff.frame) || isempty(childOff.frame)
        continue
    end
    coordName = '';
    taNodes = j.getElementsByTagName('TransformAxis');
    for t = 0:taNodes.getLength() - 1
        cn = iChildText(taNodes.item(t), 'coordinates');
        if isempty(cn)
            continue
        end
        coordName = regexprep(cn, ',', ' ');
        break
    end
    parentBody = regexprep(parentOff.frame, '_offset$', '');
    childBody = regexprep(childOff.frame, '_offset$', '');
    E = struct('name', jname, 'coord', coordName);
    jmap(sprintf('%s_to_%s', parentBody, childBody)) = E;
    if strcmp(jname, 'ankle_r')
        jmap('tibia_r_to_calcn_r') = E;   % builder note 6 alias
    end
end
end

function c = iCoordForPair(jmap, b1, b2)
k1 = sprintf('%s_to_%s', b1, b2);
k2 = sprintf('%s_to_%s', b2, b1);
if isKey(jmap, k1)
    c = jmap(k1).coord;
elseif isKey(jmap, k2)
    c = jmap(k2).coord;
else
    c = '';
end
end

function off = iOffsetFrame(jointNode, socketTag)
off = struct('frame', '');
socks = jointNode.getElementsByTagName(socketTag);
if socks.getLength() == 0
    return
end
target = regexprep(strtrim(char(socks.item(0).getFirstChild().getData())), ...
    '.*/', '');
frames = jointNode.getElementsByTagName('PhysicalOffsetFrame');
for i = 0:frames.getLength() - 1
    f = frames.item(i);
    if strcmp(char(f.getAttribute('name')), target)
        off.frame = target;
        return
    end
end
end

function childText = iChildText(node, tagName)
childText = '';
nodes = node.getElementsByTagName(tagName);
if nodes.getLength() == 0
    return
end
n0 = nodes.item(0);
if ~n0.hasChildNodes()
    return
end
childText = strtrim(char(n0.getFirstChild().getData()));
end

%% ------------ Human curve loading (single source of truth) -----------
function curves = loadHumanCurves(csvPath)
% Loads human_torque_curves_20261022.csv (written by the concurrent
% OpenSim stage, header inspected 2026-10-22). Layout: WIDE, one row
% per actuator -- actuator_id, master_name, primary_dof,
% theta_deg_01..theta_deg_25, tau_h_Nm_01..tau_h_Nm_25,
% artifact_flag_caveat -- with '#' comment lines and quoted names
% (master names contain commas). The tau_h block is the summed
% human_group curve (arm x Fmax, OpenSim native signs) on the 25-point
% rom_table_deg grid. Returns a containers.Map keyed by actuator id
% (double) -> struct(.dof, .anglesDeg 1xN, .tauNm 1xN).
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
cols = raw.Properties.VariableNames;
angCols = iOrderedCols(cols, 'theta_deg_');
tauCols = iOrderedCols(cols, 'tau_h_Nm_');
if isempty(angCols) || isempty(tauCols) || ...
        ~any(strcmpi(cols, 'actuator_id')) || ...
        ~any(strcmpi(cols, 'primary_dof'))
    error('Opt_run_BiPulley_survey_20261022:csvFormat', ...
        ['Unrecognized human CSV layout. Columns found: %s. Need ' ...
         'actuator_id/primary_dof/theta_deg_*/tau_h_Nm_* columns.'], ...
        strjoin(cols, ', '))
end
if numel(angCols) ~= numel(tauCols)
    error('Opt_run_BiPulley_survey_20261022:csvGrid', ...
        'Angle (%d) and torque (%d) column counts differ.', ...
        numel(angCols), numel(tauCols))
end
curves = containers.Map('KeyType', 'double', 'ValueType', 'any');
for r = 1:height(raw)
    idn = double(raw.actuator_id(r));
    if isnan(idn)
        continue
    end
    angles = zeros(1, numel(angCols));
    taus = zeros(1, numel(tauCols));
    for c = 1:numel(angCols)
        angles(c) = raw.(angCols{c})(r);
        taus(c) = raw.(tauCols{c})(r);
    end
    e = struct('dof', strtrim(raw.primary_dof{r}), ...
        'anglesDeg', angles, 'tauNm', taus);
    if isKey(curves, idn)
        curves(idn) = [curves(idn), e];
    else
        curves(idn) = e;
    end
end
end

function picked = iOrderedCols(cols, prefix)
% Columns named <prefix><2-digit index>, returned in index order.
idx = [];
for c = 1:numel(cols)
    if strncmpi(cols{c}, prefix, numel(prefix))
        tail = strtrim(cols{c}(numel(prefix) + 1:end));
        n = str2double(tail);
        if ~isnan(n)
            idx(end + 1, :) = [n, c];   %#ok<AGROW>
        end
    end
end
idx = sortrows(idx, 1);
picked = cols(idx(:, 2));
end

function HC = iCurvesFor(curves, actId, masterName, dof) %#ok<INUSL>
HC = [];
if isempty(curves)
    return
end
if isKey(curves, double(actId))
    es = curves(double(actId));
else
    return
end
for e = 1:numel(es)
    if strcmp(es(e).dof, dof)
        HC = es(e);
        return
    end
end
end

%% ------------ Small helpers ------------------------------------------
function s = iSiteRow(A, seed, spec0, status, reason)
s = struct('id', A.id, 'master_name', A.master_name, ...
    'primary_dof', A.primary_dof, 'seed', seed, 'spec', spec0, ...
    'status', status, 'reason', reason, 'slot', 0, 'best', [], ...
    'evalsUsed', 0, 'wallMin', 0, 'chain', '');
end

function sites = iAppendSite(sites, row)
% Struct-array growth from an EMPTY struct([]) errors in MATLAB; branch.
if isempty(sites)
    sites = row;
else
    sites(end + 1) = row;   %#ok<AGROW>
end
end

function s = iSetStatus(s, status, reason)
s.status = status;
s.reason = reason;
end

function iAppendSurveyRow(surveyCsv, S, best)
% One checkpoint row per actuator, appended the moment the site reaches
% a terminal status (FULL-LAUNCH EDIT: converged/stopReason are
% site-level -- the LAST rung's advisor verdict and why the ladder
% stopped, not the best rung's).
fid = fopen(surveyCsv, 'a');
if isempty(best)
    fprintf(fid, ['%d,%s,%s,%s,%s,"%s",%d,%d,%d,%s,%s,%s,%s,%s,%d,%.2f,' ...
        '%d,%s,%s,%s\n'], S.id, S.master_name, S.primary_dof, S.seed, ...
        S.status, iCsvSafe(S.reason), 0, 0, 0, '', '', '', '', '', 0, 0, ...
        0, S.status, '', '');
    fclose(fid);
    return
end
fprintf(fid, ['%d,%s,%s,%s,%s,"%s",%d,%d,%d,%.6g,%.6g,%.4f,%.4f,%.4f,%d,' ...
    '%.2f,%d,%s,"%s",%s\n'], S.id, S.master_name, S.primary_dof, S.seed, ...
    'surveyed', iCsvSafe(S.reason), numel(strsplit(S.chain, ' -> ')), ...
    best.nB, best.G, best.xBest(1), best.xBest(2), best.meetPct, ...
    best.activePct, best.score, S.evalsUsed, S.wallMin, ...
    S.convergedSite, S.stopReason, iCsvSafe(S.chain), best.matFile);
fclose(fid);
end

function s = iCsvSafe(txt)
s = strrep(txt, '"', '''');
end

function s = iNumOrDash(v)
if isnan(v)
    s = '-';
else
    s = sprintf('%d', v);
end
end

function g = iBadGrade(szd, msg)
g = struct('tauB', nan(szd), 'strain', nan(szd), ...
    'slack', true(szd), 'infeas', true(szd), 'taut', false(szd), ...
    'meet', false(szd), 'active', false(szd), 'obj', [], ...
    'reason', msg);
end

function v = iEnvNum(name, def)
raw = getenv(name);
if isempty(raw)
    v = def;
else
    v = str2double(raw);
    if isnan(v)
        v = def;
    end
end
end

function mode = iModeString(runBatch, isSmoke, useParallel)
if isSmoke
    mode = 'SMOKE (1 site tib_ant_r, serial, tiny budget, end to end)';
elseif ~runBatch
    mode = 'DRY-RUN (buildability + plan; no optimization)';
elseif useParallel
    mode = 'SURVEY (all sites, parpool <= 8)';
else
    mode = 'SURVEY (all sites, serial)';
end
end

%% ------------ Rigid route lengths (Opt_run_BiPulley.m:663-681) -------
function L = rigidRouteLengths(points, cross, T, N1, N2)
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

function T = hingeTrans(jtype, axis, pivot, theta)
switch jtype
    case 'fixed'
        R = eye(3);
    otherwise
        R = axisRot(axis, theta);
end
T = RpToTrans(R, pivot(:));
end

function R = axisRot(axis, theta)
a = axis(:) / norm(axis);
Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end
