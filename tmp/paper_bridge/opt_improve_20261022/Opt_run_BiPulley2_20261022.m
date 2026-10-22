% Opt Run Bi Pulley2 20261022
% Author: Ben Bolen (base Opt_run_BiPulley.m, 2026-09-24); configuration-
% aware driver family added 2026-10-22 per OPTIMIZER_IMPROVEMENT_PLAN.md
% core section + PULLEY_SWEEP_DIRECTIVE_20261021.md.
% Description: CONFIGURATION-AWARE copy of Opt_run_BiPulley.m (doctrine
% copy-and-mark-edits convention, tmp/paper_bridge/Opt_run_BiPulley_
% doctrine_*_20261021.m). No Ben file is edited. Marked edits:
%
%   EDIT-1 (plan Item 1 / directive section 1): outer 3x3 sweep
%     nBPA {1,2,3} x tackleLineParts G {1,2,3}, config injected into the
%     single-builder manifest exactly like the doctrine copies
%     (struct('nPulleyBPA', nB, 'tackleLineParts', [G 1], 'routingMode',
%     'moving_via')); [rest; tendon] (the SIZES) re-solved per
%     configuration by the same surrogateopt+patternsearch pass -- the
%     directive's size re-solve requirement (fixed sizes go slack at
%     G >= 2, proven 2026-10-21). Dated mats
%     BiPulley_<source>_<name>_nB%d_G%d_<stamp>.mat; summary CSV with
%     nBPA/G columns; sweep33 checkpoint.
%   EDIT-2 (review blocker 2, 2026-10-22): DETERMINISTIC REPLAY
%     regression gate at 0.000 tolerance: load a recorded anchor mat
%     (2026-09-26 G = 1 baseline + 2026-10-21 doctrine), re-evaluate its
%     xBest through the UNCHANGED BiPam_pulley class exactly as
%     Opt_run_BiPulley.m:360-385 evaluate() does (deterministic: class
%     construction + arithmetic, no RNG), assert the recomputed
%     margin/worstC equal the recorded st.margin/st.worstC. This
%     REPLACES the plan's "run-to-run tolerance"/"sign and scale"
%     phrasing (a fresh stochastic run can never reproduce a recorded
%     optimum at 0.000 tolerance; the deterministic re-evaluation can).
%   EDIT-3 (plan Item 2 / directive section 5, review blocker 1): the
%     fuzzy ESCALATION LADDER as an outer loop: per site, diagnose the
%     solved object with escalation_advisor_20261022 (torque-short test
%     over TAUT FEASIBLE cells only; slack cells -- where Torque_ins is
%     identically 0 by construction, F_ins zeroed at slack,
%     BiPam_pulley.m:671 -- map to the RoM/slack side), log the choice,
%     step the configuration (torque-short -> nBPA+1, RoM/slack-short ->
%     G+1, both -> both, capped at the directive grid {3,3}), re-solve
%     the sizes at the new rung.
%   EDIT-4: integrality asserts on both config axes (the class accepts
%     any TackleLineParts >= 1, integer not enforced: BiPam_pulley.m:
%     252-260).
%   EDIT-5: results pinned to Code/Matlab/Mesh_Optimization/Results
%     (this copy lives in tmp/paper_bridge/opt_improve_20261022).
%
% UNCHANGED BY CONSTRUCTION: the objective/constraint pair is
% byte-identical to Opt_run_BiPulley.m:387-412 -- anti-gaming terms stay
% in: batchObjective = -margin + 10*max(0, worstC) with the strain-floor
% and zero-pulley-infeasible constraints; margin is the worst-cell
% insertion-torque margin over ALL cells (slack cells price in as 0
% torque, so a slack design can never read as feasible). The flexor-side
% companion copy Opt_run_pulley_sweep33_20261022.m carries the full
% J-stack (JoffAxis, JwrapRelease, JrestLength, JstrainHi/Lo,
% Jbpa2Anchor) verbatim.
%
% RUN MODES:
%   default: DRY-RUN -- resolve + print the plan (all 9 configs), run
%     the deterministic replay gate on the smoke anchor pair, nothing
%     else.
%   env OPT_BPULLEY2_MODE = 'sweep' (default) | 'escalate' | 'gate':
%     sweep    -- manifest x 9 configs, one dated mat + CSV row + advisor
%                 verdict per cell (the completeness check).
%     escalate -- per-site ladder: seed config, solve, diagnose, step,
%                 re-solve; every choice logged to console + ladder CSV.
%     gate     -- the deterministic replay gate over ALL anchor mats,
%                 no optimization.
%   env OPT_BPULLEY2_SMOKE=1: one site (med_gas_r), replay gate on the
%     anchor pair, sweep sub-grid {(1,1),(2,2)}, ladder seeded at the
%     doctrine-known hard cell (2,2) with maxRungs = 2, tiny budget,
%     serial. Proves every piece end to end in minutes.
%   RUN_BATCH = true + FULL_RUN = true: real budgets (surrogateopt 1000 +
%     patternsearch 5000, parallel) -- EASTEREGG2 ONLY (parpool <= 10).

% DRY-RUN by default. Set true to execute optimization modes.
% EDIT-6 (campaign launch 2026-10-22): env overrides so the real campaign
% launches without editing this file per run (defaults unchanged):
%   OPT_BPULLEY2_RUN=1  -> RUN_BATCH = true (execute optimization modes)
%   OPT_BPULLEY2_FULL=1 -> FULL_RUN  = true (easteregg2 real budgets,
%                          surrogateopt 1000 + patternsearch 5000, parallel)
RUN_BATCH = strcmp(getenv('OPT_BPULLEY2_RUN'), '1');

% Real optimizer budgets belong on easteregg2 (10 cores, 128 GB RAM).
FULL_RUN = strcmp(getenv('OPT_BPULLEY2_FULL'), '1');

smoke = strcmp(getenv('OPT_BPULLEY2_SMOKE'), '1');
modeSel = lower(strtrim(getenv('OPT_BPULLEY2_MODE')));
if isempty(modeSel)
    modeSel = 'sweep';
end
if ~ismember(modeSel, {'sweep', 'escalate', 'gate'})
    error('Opt_run_BiPulley2_20261022:mode', ...
        'OPT_BPULLEY2_MODE must be sweep|escalate|gate, got %s.', modeSel)
end
if smoke
    RUN_BATCH = true;
end

%% Path setup (EDIT-5: this copy lives in tmp/paper_bridge/opt_improve_20261022)
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
addpath(scriptDir);   % escalation_advisor_20261022 lives next to this copy

resultsRoot = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
if ~exist(resultsRoot, 'dir')
    mkdir(resultsRoot)
end

%% Configuration grid (EDIT-1 + EDIT-4)
gridNB = [1, 2, 3];    % nBPA in parallel feeding the tackle
gridG = [1, 2, 3];     % tackleLineParts G (travel gain)
validateattributes(gridNB, {'numeric'}, {'row', 'positive', 'integer'}, ...
    '', 'gridNB', 20261022);
validateattributes(gridG, {'numeric'}, {'row', 'positive', 'integer'}, ...
    '', 'gridG', 20261022);

% Smoke sub-grid: two cells prove the config injection changes the
% problem ((1,1) direct drive vs (2,2) parallel pair with tackle).
if smoke
    smokeCells = {[1, 1], [2, 2]};
else
    smokeCells = {};
end

N1 = 5;     % joint-1 angles (deg; spec.kin.thetaRangesDeg(1,:))
N2 = 5;     % joint-2 angles

% Shared stiffness defaults (identical to Opt_run_BiPulley.m).
baseXi0 = 0.005;        % m
baseXi1 = 5e5;          % N/m
baseXi2 = 1e4;          % N/m
baseWraps = 6;          % 20 mm cable-wrap multiplier (Spr)
basePressure = 620;     % kPa
tauArmNominal = 0.05;   % m, nominal moment arm for the torque target

%% Deterministic replay gate (EDIT-2) -- runs in EVERY mode ----------
% Anchor mats of record: the 2026-09-26 G = 1 baseline (gains 1/1) and
% the 2026-10-21 doctrine set. Read-only; asserts the UNCHANGED class
% reproduces each recorded optimum at 0.000 tolerance.
anchors = collectAnchorMats(resultsRoot);
if smoke || strcmp(modeSel, 'gate')
    anchorsSubset = anchors(contains(anchors, 'med_gas_r'));
    if isempty(anchorsSubset)          % non-easteregg2 fallback: any pair
        anchorsSubset = anchors(1:min(2, numel(anchors)));
    end
    replayGate(anchorsSubset, resultsRoot)
else
    replayGate(anchors, resultsRoot)
end

%% Resolve the manifest (single-builder; masters drive absent on easteregg2)
manifest = defaultManifest();
if smoke
    fprintf('SMOKE MODE: one site only (%s, manifest slot 1).\n', ...
        manifest(1).name)
    manifest = manifest(1);
end

fprintf('\n================ BI PULLEY2 (20261022) PLAN ================\n')
fprintf('Mode                = %s\n', modeString2(RUN_BATCH, FULL_RUN, smoke, modeSel))
fprintf('Config grid         = nBPA %s x G %s (%d cells)%s\n', ...
    mat2str(gridNB), mat2str(gridG), numel(gridNB) * numel(gridG), ...
    iSmokeCellsNote(smoke, smokeCells))
fprintf('Orientation grid    = %d x %d angles per muscle\n', N1, N2)
fprintf('Stiffness defaults  = Xi0 %g m, Xi1 %g N/m, Xi2 %g N/m, wraps %d\n', ...
    baseXi0, baseXi1, baseXi2, baseWraps)
fprintf('Torque target       = fmax * %.0f cm nominal arm (documented heuristic)\n', ...
    100 * tauArmNominal)
fprintf('Muscles in manifest = %d\n', numel(manifest))
for k = 1:numel(manifest)
    printPlanItem(manifest(k), N1, N2)
end
fprintf('=============================================================\n')

if ~RUN_BATCH
    fprintf(['DRY-RUN complete: plan resolved, replay gate executed. Set ' ...
        'RUN_BATCH = true (or OPT_BPULLEY2_SMOKE=1) to execute.\n'])
    return
end

%% Budgets
if smoke
    useParallel = false;
    budgetS = 0;       % surrogateopt stage skipped in smoke
    budgetP = 60;      % patternsearch evaluations
else
    useParallel = FULL_RUN;   % parallel only on easteregg2 (FULL_RUN)
    if FULL_RUN
        budgetS = 1000;
        budgetP = 5000;
    else
        budgetS = 30;
        budgetP = 100;
        fprintf(['NOTE: FULL_RUN = false -- check-scale budgets. Real runs ' ...
            'belong on easteregg2 with FULL_RUN = true.\n'])
    end
end
if useParallel && isempty(gcp('nocreate'))
    parpool(10);       % easteregg2 cap
end

%% ================= MODE: sweep (EDIT-1) =============================
if strcmp(modeSel, 'sweep')
    ckptFile = fullfile(resultsRoot, 'BiPulley2_checkpoint_sweep33_20261022.mat');
    doneKeys = {};
    if isfile(ckptFile)
        Sck = load(ckptFile, 'doneKeys');
        doneKeys = Sck.doneKeys;
        fprintf('Checkpoint found: %d cell(s) already done.\n', numel(doneKeys))
    end
    summaryFile = fullfile(resultsRoot, 'BiPulley_batch_summary_sweep33_20261022.csv');
    if ~isfile(summaryFile)
        fid = fopen(summaryFile, 'a');
        fprintf(fid, ['%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n'], ...
            'name', 'source', 'nBPA', 'G', 'exitflagS', 'exitflagP', ...
            'objective', 'worstConstraint', 'torqueMargin', 'slackCells', ...
            'infeasibleCells', 'verdict', 'stepNBPA', 'stepG', ...
            'pulleyGain', 'footprintWidth_mm', 'timestamp');
        fclose(fid);
    end
    adviceCsv = fullfile(resultsRoot, 'BiPulley_escalation_advice_sweep33_20261022.csv');

    for k = 1:numel(manifest)
        spec0 = manifest(k);
        cellsNow = iCellsForRun(smoke, smokeCells, gridNB, gridG);
        for c = 1:size(cellsNow, 1)
            nB = cellsNow(c, 1);
            G = cellsNow(c, 2);
            key = sprintf('%s_%s_nB%d_G%d', spec0.source, spec0.name, nB, G);
            if any(strcmp(doneKeys, key))
                fprintf('[%s]: already done (checkpoint), skipping.\n', key)
                continue
            end
            fprintf('\n===== SWEEP %s | nBPA = %d, G = %d =====\n', key, nB, G)
            try
                spec = injectConfig(spec0, nB, G);
                R = runMuscle(spec, N1, N2, baseXi0, baseXi1, baseXi2, ...
                    baseWraps, basePressure, tauArmNominal, useParallel, ...
                    budgetS, budgetP, resultsRoot, smoke, nB, G);
                % Advisor verdict per cell (the sweep stays the
                % completeness check; the advisor orders/prunes).
                [advice, ~] = escalation_advisor_20261022(R.st.obj, ...
                    R.tauTarget, struct('curNBPA', nB, 'curG', G, ...
                    'verbose', true, 'csvFile', adviceCsv));
                sfid = fopen(summaryFile, 'a');
                fprintf(sfid, ['%s,%s,%d,%d,%g,%g,%.6g,%.6g,%.6g,%d,%d,%s,' ...
                    '%d,%d,%s,%.3f,%s\n'], spec.name, spec.source, nB, G, ...
                    R.exitflagS, R.exitflagP, R.fBest, R.st.worstC, ...
                    R.st.margin, sum(R.st.obj.PulleySlack(:)), ...
                    sum(R.st.obj.PulleyInfeasible(:)), advice.verdict, ...
                    advice.stepNBPA, advice.stepG, ...
                    gainString(R.st.obj.PulleyGain), ...
                    1000 * spec.footprintWidth, ...
                    char(string(datetime('now'), 'yyyyMMdd_HHmmss')));
                fclose(sfid);
                doneKeys{end + 1} = key;   %#ok<AGROW>
                save(ckptFile, 'doneKeys');
                fprintf('Checkpoint updated (%d done).\n', numel(doneKeys))
            catch err
                fprintf('CELL FAILED (%s): %s\n%s\n', key, err.message, ...
                    getReport(err, 'basic'))
            end
        end
    end
    fprintf('\nSweep complete. See %s\n', summaryFile)
end

%% ================= MODE: escalate (EDIT-3) ==========================
if strcmp(modeSel, 'escalate')
    ladderCsv = fullfile(resultsRoot, 'BiPulley_escalation_ladder_20261022.csv');
    if ~isfile(ladderCsv)
        fid = fopen(ladderCsv, 'a');
        fprintf(fid, ['%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,' ...
            '%s,%s,%s,%s,%s,%s\n'], 'name', 'source', 'rung', 'nBPA', 'G', ...
            'sizesSolved', 'marginRaw', 'worstC', 'slackCells', ...
            'infeasibleCells', 'minTautTau', 'tautMargin', 'verdict', ...
            'stepNBPA', 'stepG', 'nextNBPA', 'nextG', 'converged', ...
            'timestamp');
        fclose(fid);
    end
    maxRungs = 6;
    if smoke
        maxRungs = 2;
    end
    for k = 1:numel(manifest)
        spec0 = manifest(k);
        if smoke
            seedNB = 2; seedG = 2;   % doctrine-known hard cell: a step is
            % guaranteed demonstrable (10 slack cells at (2,[2 1]),
            % doctrine_rerun_20261021.log:1160-1164).
        else
            seedNB = 1; seedG = 1;   % the 2026-09-26 baseline config
        end
        nB = seedNB; G = seedG;
        fprintf('\n########## LADDER %s (%s): seed (%d,%d), maxRungs %d ##########\n', ...
            spec0.name, spec0.source, nB, G, maxRungs)
        for rung = 1:maxRungs
            spec = injectConfig(spec0, nB, G);
            fprintf('\n===== RUNG %d | %s | nBPA = %d, G = %d (sizes re-solved) =====\n', ...
                rung, spec.name, nB, G)
            R = runMuscle(spec, N1, N2, baseXi0, baseXi1, baseXi2, ...
                baseWraps, basePressure, tauArmNominal, useParallel, ...
                budgetS, budgetP, resultsRoot, smoke, nB, G, ...
                'ladder_rung%d', rung);
            [advice, rep] = escalation_advisor_20261022(R.st.obj, ...
                R.tauTarget, struct('curNBPA', nB, 'curG', G, 'verbose', true));
            lfid = fopen(ladderCsv, 'a');
            fprintf(lfid, ['%s,%s,%d,%d,%d,%d,%.6g,%.6g,%d,%d,%.6g,%.6g,%s,' ...
                '%d,%d,%d,%d,%d,%s\n'], spec.name, spec.source, rung, nB, G, ...
                1, R.st.margin, R.st.worstC, ...
                sum(R.st.obj.PulleySlack(:)), ...
                sum(R.st.obj.PulleyInfeasible(:)), rep.minTautTau, ...
                rep.tautMargin, advice.verdict, advice.stepNBPA, ...
                advice.stepG, advice.nextNBPA, advice.nextG, ...
                advice.converged, char(string(datetime('now'), ...
                'yyyyMMdd_HHmmss')));
            fclose(lfid);
            fprintf(['LADDER CHOICE (rung %d, %s): verdict = %s -> steps ' ...
                '(nBPA,G) = (+%d,+%d)%s | %s\n'], rung, spec.name, ...
                advice.verdict, advice.stepNBPA, advice.stepG, ...
                iNextStr(advice), advice.reason);
            if advice.converged
                fprintf('LADDER %s: converged at rung %d (config (%d,%d)).\n', ...
                    spec.name, rung, nB, G)
                break
            end
            nB = nB + advice.stepNBPA;
            G = G + advice.stepG;
        end
    end
    fprintf('\nEscalation ladders complete. See %s\n', ladderCsv)
end

fprintf('\nOpt_run_BiPulley2_20261022 done (mode %s, smoke %d).\n', ...
    modeSel, smoke)

%% =====================================================================
%% Local functions
%% =====================================================================
function mode = modeString2(runBatch, fullRun, isSmoke, modeSel) %#ok<INUSD>
if ~runBatch
    mode = 'DRY-RUN (plan + replay gate only)';
elseif isSmoke
    mode = sprintf('SMOKE (%s path, 1 site, tiny budget, serial)', modeSel);
elseif fullRun
    mode = sprintf('FULL RUN (%s, easteregg2 budgets, parallel)', modeSel);
else
    mode = sprintf('CHECK RUN (%s, small budgets, serial)', modeSel);
end
end

function note = iSmokeCellsNote(isSmoke, smokeCells)
if ~isSmoke
    note = '';
    return
end
parts = cellfun(@(c) sprintf('(%d,%d)', c(1), c(2)), smokeCells, ...
    'UniformOutput', false);
note = [' -- SMOKE sub-grid: ', strjoin(parts, ' ')];
end

function cellsNow = iCellsForRun(isSmoke, smokeCells, gridNB, gridG)
if isSmoke
    % vertcat, NOT cell2mat: cell2mat({[1 1],[2 2]}) flattens to a 1x4
    % row (the exact column-vector trap Opt_run_pulley.m:196-199 warns
    % about -- one iteration binding the whole vector).
    cellsNow = vertcat(smokeCells{:});
    return
end
cellsNow = zeros(numel(gridNB) * numel(gridG), 2);
c = 0;
for nB = gridNB
    for G = gridG
        c = c + 1;
        cellsNow(c, :) = [nB, G];
    end
end
end

function s = iNextStr(advice)
if isnan(advice.nextNBPA) && isnan(advice.nextG)
    s = ' (at directive grid cap)';
else
    s = sprintf(' -> next (%s,%s)', iNumStr(advice.nextNBPA), ...
        iNumStr(advice.nextG));
end
end

function s = iNumStr(v)
if isnan(v)
    s = '-';
else
    s = sprintf('%d', v);
end
end

function manifest = defaultManifest()
% Single-builder manifest (the doctrine-copy pattern; the master's data
% drive is absent on easteregg2). The CONFIG is injected per cell/rung by
% injectConfig; the builder default here only resolves names/geometry.
manifest = [];
try
    manifest = biPulleySpecsFromOpenSim();   %#ok<AGROW>
catch err
    fprintf('OSIM BUILDER UNAVAILABLE: %s\n', err.message)
end
if isempty(manifest)
    error('Opt_run_BiPulley2_20261022:noManifest', ...
        'The builder produced no specs; check the data sources.')
end
end

function spec = injectConfig(spec0, nB, G)
% Doctrine-copy injection, generalized: per-crossing gains [G 1]
% (BiPam_pulley TwoTackles rejects both > 1; doctrine convention has the
% tackle at crossing 1). nBPA = 1 has NO bundle (single line); nBPA > 1
% re-derives the symmetric bundle exactly as the builder does.
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

function anchors = collectAnchorMats(resDir)
% Anchor mats of record (READ-ONLY): 2026-09-26 G = 1 baseline + the
% 2026-10-21 doctrine set (plan Item 1 / blocker 2).
d = [dir(fullfile(resDir, 'BiPulley_opensim_*_20260926_*.mat')); ...
    dir(fullfile(resDir, 'BiPulley_opensim_*_2BPA_G2_20261021_*.mat')); ...
    dir(fullfile(resDir, 'BiPulley_opensim_*_2BPA_G4_20261021_*.mat')); ...
    dir(fullfile(resDir, 'BiPulley_opensim_*_4BPA_G4_20261021_*.mat'))];
anchors = cellfun(@(f) fullfile(resDir, f), {d.name}, ...
    'UniformOutput', false);
anchors = anchors(:);
if isempty(anchors)
    error('Opt_run_BiPulley2_20261022:noAnchors', ...
        'No anchor mats found in %s', resDir)
end
end

function replayGate(anchorFiles, resDir) %#ok<INUS>
% EDIT-2: deterministic 0.000-tolerance replay. Load each anchor mat,
% re-evaluate its recorded xBest through the UNCHANGED BiPam_pulley class
% exactly as Opt_run_BiPulley.m:360-385 evaluate() does, and assert the
% recomputed margin/worstC EQUAL the recorded st.margin/st.worstC.
% Deterministic by construction (festo4 + fzero inside the class, no
% RNG), which is why 0.000 tolerance is meaningful here -- unlike a
% fresh stochastic rerun.
nFail = 0;
fprintf('\n---- DETERMINISTIC REPLAY GATE (0.000 tolerance) ----\n')
for a = 1:numel(anchorFiles)
    f = anchorFiles{a};
    S = load(f);
    rest = S.xBest(1);
    tendon = S.xBest(2);
    kmax = 0.75 * rest;      % Opt_run_BiPulley.m:363
    obj = BiPam_pulley(S.spec.name, S.points, S.spec.cross, 20, S.T, ...
        rest, kmax, tendon, 0.025, S.basePressure, ...
        S.baseXi0, S.baseXi1, S.baseXi2, S.baseWraps, S.spec.pulley);
    aK = obj.ActiveCrossing;
    tauAbs = abs(obj.Torque_ins(:, :, :, aK));
    validTau = tauAbs(~isnan(tauAbs));
    margin = min(validTau) / S.tauTarget - 1;
    strainGrid = reshape(obj.strain_p, size(S.Lrig));
    c1 = obj.strainFloor() - min(strainGrid(:));
    c2 = double(sum(obj.PulleyInfeasible(:)));
    worstC = max([c1, c2, 0]);
    dMargin = margin - S.st.margin;
    dWorstC = worstC - S.st.worstC;
    [pt, ~, ext] = fileparts(f);
    okGate = (dMargin == 0) && (dWorstC == 0);   % 0.000 tolerance
    if okGate
        fprintf(['  PASS %-58s margin %+.17g (d = %+g), worstC %+.17g ' ...
            '(d = %+g)\n'], [pt, ext], margin, dMargin, worstC, dWorstC)
    else
        nFail = nFail + 1;
        fprintf(['  FAIL %-58s margin %+.17g vs recorded %+.17g (d = ' ...
            '%+g), worstC %+.17g vs recorded %+.17g (d = %+g)\n'], ...
            [pt, ext], margin, S.st.margin, dMargin, worstC, ...
            S.st.worstC, dWorstC)
    end
end
if nFail > 0
    error('Opt_run_BiPulley2_20261022:replayGate', ...
        'Deterministic replay gate FAILED on %d of %d anchor mat(s).', ...
        nFail, numel(anchorFiles))
end
fprintf('  REPLAY GATE PASS: %d of %d anchor mat(s) reproduced exactly.\n', ...
    numel(anchorFiles), numel(anchorFiles))
end

function printPlanItem(spec, N1, N2)
fprintf(['  %-11s %-7s pts=%d cross=[%d %d] frames=[%s %s %s] Fmax=%5g N ' ...
    'mode=%s\n'], spec.name, spec.source, size(spec.points, 1), ...
    spec.cross(1), spec.cross(2), clean(spec.frames{1}), ...
    clean(spec.frames{2}), clean(spec.frames{3}), spec.fmax, ...
    char(spec.pulley.routingMode));
fprintf('    joints: [%s; %s]  ranges(deg): [%g %g; %g %g]\n', ...
    spec.kin.types{1}, spec.kin.types{2}, ...
    spec.kin.thetaRangesDeg(1, :), spec.kin.thetaRangesDeg(2, :));
fprintf('    grid %dx%d angles; fmaxSource: %s\n', N1, N2, spec.fmaxSource);
fprintf('    notes: %s\n', spec.notes);
end

function points = pointsWithFold(spec) %#ok<INUSD>
% The route points the driver actually feeds the class (osim subtalar
% fold applied), as in Opt_run_BiPulley.m:292-301.
points = spec.points;
if isfield(spec.kin, 'pointFolds') && ...
        any(any(abs(spec.kin.pointFolds{2}) > 0))
    points(spec.cross(2):end, :) = ...
        points(spec.cross(2):end, :) + spec.kin.pointFolds{2}';
end
end

function R = runMuscle(spec, N1, N2, baseXi0, baseXi1, baseXi2, ...
    baseWraps, basePressure, tauArmNominal, useParallel, budgetS, ...
    budgetP, resDir, isSmoke, nB, G, varargin)
% One muscle end to end at ONE configuration: transforms, sizing-bound
% construction, optimization (sizes = [rest; tendon] RE-SOLVED here),
% artifacts. Byte-identical objective/constraints to
% Opt_run_BiPulley.m:303-503 except the marked artifact-name edits.
% Returns R with st, xBest, fBest, exitflags, tauTarget so the ladder
% outer loop can diagnose without touching the optimizer.

if ~isempty(varargin)
    tagFmt = varargin{1};      % e.g. 'ladder_rung%d'
    tagArg = varargin{2};
    cellTag = sprintf('_%s_nB%d_G%d', sprintf(tagFmt, tagArg), nB, G);
else
    cellTag = sprintf('_nB%d_G%d', nB, G);
end

% --- Two-hinge transform grid (BiPamData slabs) ------------------------
[T, ~, ~] = buildBiTransforms(spec, N1, N2);

% --- Route points with the osim subtalar fold --------------------------
points = pointsWithFold(spec);

% --- Characteristic route length for sizing/bounds --------------------
Lrig = rigidRouteLengths(points, spec.cross, T, N1, N2);
Lchar = max(Lrig(:));

% --- Interference checks: bundle self-clearance ------------------------
fprintf(['Clearance checks (BPA diameter %d mm, margin %.1f mm):\n'], ...
    spec.diameter, 1000 * 0.005)
if size(spec.bundleOffsets, 1) > 1
    nOff = size(spec.bundleOffsets, 1);
    minSelf = inf;
    for a = 1:nOff - 1
        for b = a + 1:nOff
            minSelf = min(minSelf, norm(spec.bundleOffsets(b, :) - ...
                spec.bundleOffsets(a, :)));
        end
    end
    needSelf = spec.diameter * 1e-3 + 0.005;
    fprintf(['  bundle self-clearance: min member spacing %.1f mm ' ...
        '(need >= %.1f mm): %s\n'], 1000 * minSelf, 1000 * needSelf, ...
        verdictString(minSelf >= needSelf))
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
    [st, ok] = evaluate(xq);
    if ~ok
        c = 1e3;
    else
        c = st.worstC;
    end
    ceq = [];
    end

    function [fVal, c, ceq] = objconstrWrap(xq)
    fVal = batchObjective(xq);
    c = batchNonlcon(xq);
    ceq = [];
    end

x0 = [0.60 * Lchar; 0.05];
f0 = batchObjective(x0);
fprintf('Initial: x = [%.4f %.4f] m, penalized margin = %.4f\n', x0, f0)

% --- Optimization (sizes re-solved per configuration) ------------------
x = x0;
f = f0;
exitflagS = NaN;
exitflagP = NaN;

% Stage 1: surrogateopt global pass (skipped in smoke).
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
    error('Opt_run_BiPulley2_20261022:finalEval', ...
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
fprintf('Transmission: nPulleyBPA = %s, tackleLineParts = %s, gain = %s, mode = %s\n', ...
    gainString(finalObj.NPulleyBPA), gainString(finalObj.TackleLineParts), ...
    gainString(finalObj.PulleyGain), strjoin(finalObj.RoutingMode, '/'))
fprintf('Bundle footprint width = %.1f mm (diameter %d mm)\n', ...
    1000 * spec.footprintWidth, spec.diameter)

% Dated FULL-WORKSPACE result capture (Ben directive: bare save behind
% liveRun; cleared just before saving). SMOKE never saves.
liveRun = ~isSmoke;
if exist('liveRun', 'var') && liveRun
    stamp = char(string(datetime('now'), 'yyyyMMdd_HHmm'));
    resultFile = fullfile(resDir, ...
        sprintf('BiPulley_%s_%s%s_%s.mat', spec.source, spec.name, ...
        cellTag, stamp));
    clear liveRun
    save(resultFile)
    fprintf('Saved %s\n', resultFile)
end

R = struct('st', st, 'xBest', xBest, 'fBest', fBest, ...
    'exitflagS', exitflagS, 'exitflagP', exitflagP, ...
    'tauTarget', tauTarget, 'nB', nB, 'G', G, 'Lrig', Lrig, 'T', T, ...
    'points', points, 'spec', spec);
end

function s = verdictString(ok)
if ok
    s = 'OK';
else
    s = 'VIOLATION';
end
end

function [T, theta1, theta2] = buildBiTransforms(spec, N1, N2)
% Identical to Opt_run_BiPulley.m:585-634.
kin = spec.kin;
r1 = rangeOrDefault(kin, 1, [-30, 90]);
r2 = rangeOrDefault(kin, 2, [-110, 20]);
theta1 = linspace(r1(1), r1(2), N1) * pi / 180;
theta2 = linspace(r2(1), r2(2), N2) * pi / 180;

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

function r = rangeOrDefault(kin, slot, def)
if isfield(kin, 'thetaRangesDeg') && size(kin.thetaRangesDeg, 1) >= slot
    r = kin.thetaRangesDeg(slot, :);
else
    r = def;
end
end

function L = rigidRouteLengths(points, cross, T, N1, N2)
% Identical to Opt_run_BiPulley.m:663-681.
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
