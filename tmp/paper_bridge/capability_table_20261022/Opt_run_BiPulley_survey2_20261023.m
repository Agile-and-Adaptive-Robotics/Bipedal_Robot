% Opt Run Bi Pulley survey2 20261023
% Author: Ben Bolen (machinery of record); survey2 driver assembled
% 2026-10-23 by the BiPulley survey engineer per the capability-table v2
% optimized-row ask (FIX_SPEC_20261023.md section 4, workflow 2026-10-23).
%
% Description: NEW DATED DRIVER (copy-and-mark-edits convention; no file
% of Ben's edited, no v1 artifact touched). Dated copy of
%   tmp/paper_bridge/capability_table_20261022/Opt_run_BiPulley_survey_20261022.m
% (the 1-D survey driver) RE-GRAFTED onto the 2-D machinery:
%   tmp/paper_bridge/capability_table_20261022/biPulleySpecsFromMasters_20261023.m
%   (extended masters builder, FIX_SPEC section 3; ALL 27 routes; the
%   patella rule; rowFolds applied by the builder, driver pointFolds zero
%   by contract),
%   Code/Matlab/Mesh_Optimization/BiPam_pulley2_20261022.m (variant class,
%   default config = bit-identical base arithmetic),
%   tmp/paper_bridge/opt_improve_20261022/escalation_advisor_20261022.m
%   (fuzzy ladder: taut-cell torque-short -> +1 BPA, slack/infeasible ->
%   +1 line part),
%   tmp/old_code_copy_20261022/PAM_Prop.mat (the masters data of record).
%
% SURVEY2 CONTRACT (the ask + FIX_SPEC section 4, 2026-10-23) -- every
% numbered edit below is marked [S2-n] at its site:
%   [S2-1] GRADING = the 2-D SINGLE-MUSCLE SURFACES,
%      human_surfaces_5x5_grading_20261023.csv (5x5 two-DOF grid per
%      actuator, OpenSim native signs, in_graded_domain envelope mask) --
%      the SINGLE source of truth for BOTH capability-table rows. Every
%      candidate is graded cell-by-cell against tau_h at the SAME pose:
%      primary nodes on the slot that models the primary dof, secondary
%      nodes on the other slot when that dof is modeled; sites whose
%      secondary dof is NOT one of the class's two hinges (bifemsh_r /
%      vas_int_r rigid clone; soleus_r knee; psoas_r hip_adduction) pin
%      the remaining joint at 0 and the robot surface is CONSTANT along
%      that axis -- exactly the hand-tuned row's own old-code behavior
%      (ht04_grade_hand_tuned_2d_20261023.m header; soleus documented
%      as-is there) -- recorded in the row note, never silently.
%      in_graded_domain = 0 cells (hip flexion-adduction coupled envelope,
%      FIX_SPEC 1b; acts 24-27 only, 10 cells each) are EXCLUDED from the
%      area denominator for BOTH metrics and from the advisor -- never
%      counted as failures (envelope guard).
%   [S2-2] OBJECTIVE (soft, nothing rejected): maximize MEET AREA share
%      FIRST, ACTIVE AREA share second:
%          f = -(1000*meetFrac + activeFrac) + 1e-3*deadFrac
%      over DOMAIN cells (25 minus excluded); meetFrac/activeFrac = the
%      10-22 survey's exact booleans (grade() in
%      Opt_run_BiPulley_survey_20261022.m:524-527) masked by the domain;
%      the 1e-3 dead-cell price is optimizer-only tie-break scale (the
%      all-dead plateau guard); REPORTED scores stay pure
%      -(1000*meetFrac + activeFrac). NO hard constraints, NO nonlcon:
%      pulley-infeasible cells are priced softly as not-active/meeting.
%   [S2-3] WINNER = the BEST-MEETING design (FIX_SPEC section 4: "an
%      active-percentage tie-break must never record a lower-meeting
%      design"). Rung results compare LEXICOGRAPHICALLY on INTEGER domain
%      cell counts -- meetCells desc, then activeCells desc, then simpler
%      config (nB+G) asc, then earlier rung -- never on the float score.
%   [S2-4] BPA SIZE per muscle = its OLD-CODE size, BOTH rows same
%      hardware (Ben, 2026-10-23): inside diameter from the Size.m rule
%      verbatim (Fc = 0.9 on the Festo theoretical maxima F1 630 / F2
%      1500 / F3 6000 N) applied to the muscle's PAM_Prop.mat MIF (the
%      builder's own spec.fmax load). Resting length stays a DESIGN
%      VARIABLE re-solved per configuration.
%   [S2-5] CONFIG free: nBPA {1,2,3} x tackle G {1,2,3} via the
%      escalation ladder (seed (1,1), advisor steps +1 BPA / +1 line
%      part, caps {3,3}, at most 5 rungs); sizes [rest; tendon]
%      re-solved per rung (PULLEY_SWEEP_DIRECTIVE_20261021); RESTING
%      LENGTH CEILING IS AN UPPER LIMIT: rest in [0.30*Lchar, 0.90*Lchar]
%      (solver may shorten, never lengthen), tendon in
%      [0.02, max(0.30*Lchar, 0.06)] (10-22 bounds unchanged); Lchar =
%      max rigid route length over the 2-D graded grid.
%   [S2-6] CAPS UNCHANGED: <= 300 optimizer evaluations per muscle TOTAL,
%      ~10 min wall per muscle (checked between rungs), parpool <= 8.
%   [S2-7] PRESSURE 620 kPa everywhere.
%   [S2-8] OUTPUTS: dated per-muscle mats in Code/Matlab/Mesh_Optimization/
%      Results/ named BiPulley_survey2_<site>_nB<n>_G<g>_<stamp>.mat
%      (new names, never overwritten; <site> = the masters spec name);
%      per-muscle results CSV
%      BiPulley_survey2_perMuscle_20261023.csv (this folder) with ONE ROW
%      PER ACTUATOR (27, built:false rows included) and the ladder CSV
%      BiPulley_survey2_ladder_20261023.csv. RESUME-SAFE: the perMuscle
%      CSV is the checkpoint (a restarted driver skips finished
%      actuators; per-muscle caps are per-launch). Smoke runs keep
%      stamped throwaway CSV names.
%   Route points are NOT design variables (masters route verbatim);
%   nothing is rejected -- every actuator gets a row.
%
% The deterministic replay gate is NOT re-run here: this workflow's
% verify phase already passed it for this machinery (ask, 2026-10-23);
% the surfaces CSV is the graded input, read-only.
%
% RUN MODES (env):
%   default             DRY-RUN: build all 27 sites, parse + verify the
%                       surfaces, print the plan, exit. No optimization.
%   OPT2_SURVEY_RUN=1   execute the survey (all sites; requires CSV).
%   OPT2_SURVEY_SMOKE=1 one site (tib_ant_r), serial, tiny budget
%                       (15 evals x 2 rungs), end to end incl. grading;
%                       forces RUN.
%   OPT2_SURVEY_PARALLEL=1  parpool(min(8, cores)) for the full run.
%   OPT2_SURVEY_MAXEVAL   per-muscle total eval cap (default 300).
%   OPT2_SURVEY_RUNGEVAL  per-rung eval budget (default 60; smoke 15).
%   OPT2_SURVEY_TIMEMIN   per-muscle wall budget, minutes (default 10).
%   OPT2_SURVEY_SITES     comma-separated masters spec names to restrict.

smoke = strcmp(getenv('OPT2_SURVEY_SMOKE'), '1');
RUN_BATCH = smoke || strcmp(getenv('OPT2_SURVEY_RUN'), '1');
useParallel = ~smoke && strcmp(getenv('OPT2_SURVEY_PARALLEL'), '1');

maxEvalsSite = iEnvNum('OPT2_SURVEY_MAXEVAL', 300);
rungEval = iEnvNum('OPT2_SURVEY_RUNGEVAL', 60);
timeBudgetMin = iEnvNum('OPT2_SURVEY_TIMEMIN', 10);
maxRungs = 5;
if smoke
    rungEval = iEnvNum('OPT2_SURVEY_RUNGEVAL', 15);
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
addpath(genpath(fullfile(root, 'Code', 'Matlab')));   % class + festo4 + helpers
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'tmp', 'paper_bridge', 'opt_improve_20261022')); ...
    % escalation_advisor_20261022 lives there
addpath(scriptDir);                                   % masters builder 20261023

resultsRoot = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
if ~exist(resultsRoot, 'dir')
    mkdir(resultsRoot)
end
outDir = scriptDir;   % driver-level CSVs live next to this driver

%% Shared physical constants (Opt_run_BiPulley.m defaults, verbatim) ----
baseXi0 = 0.005;          % m
baseXi1 = 5e5;            % N/m
baseXi2 = 1e4;            % N/m
baseWraps = 6;            % cable-wrap multiplier (Spr)
basePressure = 620;       % kPa (ask item: 620 kPa everywhere)
baseFitn = 0.025;         % one fitting length, m

%% [S2-1] Surfaces CSV (single source of truth for BOTH rows) -----------
surfCsv = fullfile(scriptDir, 'human_surfaces_5x5_grading_20261023.csv');
if ~isfile(surfCsv)
    error('Opt_run_BiPulley_survey2_20261023:noSurfaces', ...
        'Grading surfaces not found at %s', surfCsv)
end
surfaces = loadSurfaces(surfCsv);
fprintf('Surfaces loaded: %d actuator surface(s) from %s\n', ...
    surfaces.Count, surfCsv)

%% Resolve every site through the masters builder ------------------------
fprintf('\n================ BI PULLEY SURVEY2 (20261023) PLAN ================\n')
fprintf('Mode          = %s\n', iModeString(RUN_BATCH, smoke, useParallel))
fprintf('Budgets       = %d evals/muscle TOTAL, %d evals/rung, %d rungs max, %.0f min wall\n', ...
    maxEvalsSite, rungEval, maxRungs, timeBudgetMin)
fprintf('Config ladder = seed (1,1), advisor steps +1 BPA / +1 line part, cap {3,3}\n')
fprintf('Grading       = 2-D 5x5 surfaces, domain-masked (hip envelope excluded)\n')
fprintf('BPA size      = OLD-CODE per-muscle diameter (Size rule on PAM_Prop MIF)\n')
fprintf('Pressure      = %d kPa; sizes [rest; tendon] re-solved per rung\n', ...
    basePressure)

siteFilter = strsplit(strtrim(getenv('OPT2_SURVEY_SITES')), ',');
hasFilter = ~isempty(strtrim(getenv('OPT2_SURVEY_SITES')));

% All 27 spec names in actuator-map order (builder's own roster).
roster = biPulleySpecsFromMasters_20261023();
allNames = {roster.name};
allIds = [roster.actuatorId];
clear roster

sites = struct([]);
nBuild = 0; nCannot = 0;
for a = 1:numel(allNames)
    nm = allNames{a};
    if hasFilter && ~any(strcmp(siteFilter, nm))
        continue
    end
    if smoke && ~strcmp(nm, 'tib_ant_r')   % smoke site (ask, 2026-10-23)
        continue
    end
    row = struct('id', allIds(a), 'spec_name', nm, 'status', '', ...
        'reason', '', 'spec', [], 'surf', [], 'dia', NaN, 'note', '', ...
        'best', [], 'evalsUsed', 0, 'wallMin', 0, 'chain', '', ...
        'stopReason', '', 'convergedSite', false);
    try
        spec0 = biPulleySpecsFromMasters_20261023({nm});
        row.spec = spec0(1);
    catch err
        row.status = 'cannot_build';
        row.reason = err.message;
        nCannot = nCannot + 1;
        fprintf('[%2d] %-12s CANNOT BUILD: %s\n', row.id, nm, err.message)
        sites = iAppendSite(sites, row);
        continue
    end
    sp = row.spec;

    % [S2-4] old-code diameter (Size.m rule verbatim on the builder's own
    % PAM_Prop MIF load). 'Use additional muscles' -> cannot build (and
    % per the ask, a built:false row must be escalated by the operator).
    [dia, diaWhy] = oldCodeDiameter(sp.fmax);
    if isnan(dia)
        row.status = 'cannot_build';
        row.reason = sprintf('old-code size rule: %s', diaWhy);
        nCannot = nCannot + 1;
        fprintf('[%2d] %-12s CANNOT BUILD (size rule): %s\n', row.id, nm, diaWhy)
        sites = iAppendSite(sites, row);
        continue
    end
    row.dia = dia;

    % [S2-1] surface + dof mapping: primary must be a modeled hinge at the
    % builder's primarySlot; secondary modeled or constant-axis caveat.
    [row.surf, mapNote] = surfaceForSite(surfaces, sp, dia);
    if strcmp(mapNote, 'no surface')
        row.status = 'no_surface';
        row.reason = sprintf('no surfaces rows for actuator %d', sp.actuatorId);
        nCannot = nCannot + 1;
        fprintf('[%2d] %-12s NO SURFACE ROWS\n', row.id, nm)
        sites = iAppendSite(sites, row);
        continue
    end
    if row.surf.priSlot == 0
        row.status = 'grade_unable';
        row.reason = mapNote;
        nCannot = nCannot + 1;
        fprintf('[%2d] %-12s GRADE-UNABLE: %s\n', row.id, nm, mapNote)
        sites = iAppendSite(sites, row);
        continue
    end
    row.note = mapNote;
    row.status = 'ready';
    nBuild = nBuild + 1;
    fprintf(['[%2d] %-12s %-18s x %-18s pri slot %d; dia %d mm; ' ...
        'domain %d/25; %s\n'], row.id, nm, row.surf.priDof, ...
        row.surf.secDof, row.surf.priSlot, dia, row.surf.nDom, mapNote)
    sites = iAppendSite(sites, row);
end
fprintf('====================================================================\n')
fprintf('Build summary: %d ready, %d cannot-build/grade-unable/no-surface\n', ...
    nBuild, nCannot)

if ~RUN_BATCH
    fprintf(['DRY-RUN complete (no optimization). Set OPT2_SURVEY_RUN=1 ' ...
        '(or OPT2_SURVEY_SMOKE=1) to execute the survey.\n'])
    return
end

%% Parpool ([S2-6]: at most 8) -------------------------------------------
if useParallel
    want = min(8, feature('numcores'));
    if isempty(gcp('nocreate'))
        parpool(want);
    else
        fprintf('parpool already running (%d workers).\n', ...
            numWorkers(gcp('nocreate')))
    end
end

%% Driver-level CSVs ([S2-8]: resume checkpoint) --------------------------
if smoke
    stamp0 = char(string(datetime('now'), 'yyyyMMdd_HHmmss'));
    surveyCsv = fullfile(outDir, sprintf('BiPulley_survey2_smoke_results_%s.csv', stamp0));
    ladderCsv = fullfile(outDir, sprintf('BiPulley_survey2_smoke_ladder_%s.csv', stamp0));
else
    surveyCsv = fullfile(outDir, 'BiPulley_survey2_perMuscle_20261023.csv');
    ladderCsv = fullfile(outDir, 'BiPulley_survey2_ladder_20261023.csv');
end
if ~isfile(surveyCsv)
    sFid = fopen(surveyCsv, 'w');
    fprintf(sFid, ['actuator_id,master_name,site,built,primary_dof,secondary_dof,' ...
        'bpa_size_mm,domain_cells,excluded_cells,rungs,chosen_nB,chosen_G,' ...
        'rest_m,tendon_m,meetPct,activePct,meetPct_unflagged,activePct_unflagged,' ...
        'sign_agree_pct,score,evalsUsed,wallMin,converged,stopReason,' ...
        'ladderChain,matFile,note\n']);
    fclose(sFid);
end
if ~isfile(ladderCsv)
    lFid = fopen(ladderCsv, 'w');
    fprintf(lFid, ['actuator_id,master_name,site,rung,nBPA,G,evals,rest_m,' ...
        'tendon_m,meetCells,activeCells,domainCells,meetPct,activePct,score,' ...
        'slackCells,infeasibleCells,tautCells,minTautTau_Nm,targetTau_Nm,' ...
        'verdict,stepNBPA,stepG,converged,smoke,timestamp\n']);
    fclose(lFid);
end
doneIds = [];
if ~smoke && isfile(surveyCsv)
    Tprev = readtable(surveyCsv, 'VariableNamingRule', 'preserve', ...
        'Delimiter', ',', 'TextType', 'char');
    doneIds = double(Tprev.actuator_id(:));
    fprintf('RESUME: %d actuator(s) already recorded in %s\n', ...
        numel(doneIds), surveyCsv)
end

%% =================== THE SURVEY (per muscle) ==========================
for a = 1:numel(sites)
    S = sites(a);
    if any(doneIds == S.id)
        fprintf('\n[%d] %s already recorded (checkpoint) -- skipping.\n', ...
            S.id, S.spec_name)
        continue
    end
    if ~strcmp(S.status, 'ready')
        % recorded and skipped, never fatal; built:false rows carry the
        % reason and are escalated by the operator per the ask.
        iAppendSurveyRow(surveyCsv, S);
        continue
    end
    fprintf('\n########## SURVEY2 %s [%d] ##########\n', S.spec_name, S.id)
    tSite = tic;
    evalsUsed = 0;

    sp = S.spec;
    SU = S.surf;                          % [S2-1] the 2-D surface
    th1v = SU.th1; th2v = SU.th2;         % per-slot angle vectors (deg)
    [T, N1, N2] = gridT(sp, th1v, th2v);  % graded 2-D transforms
    points = sp.points;                   % rowFolds applied by the builder
    assert(~any(any(abs(sp.kin.pointFolds{1}) > 0)) && ...
           ~any(any(abs(sp.kin.pointFolds{2}) > 0)), ...
           'masters builder contract: driver pointFolds must be zero')
    Lrig = rigidRouteLengths(points, sp.cross, T, N1, N2);
    Lchar = max(Lrig(:));                 % [S2-5] ceiling basis

    tauHuman = SU.H;                      % 5x5, slot-dim layout
    dom = SU.dom;                         % logical 5x5
    nDom = SU.nDom;
    flg = SU.flg;
    W = double(dom);                      % [S2-1] envelope guard for the
                                          % advisor: excluded cells count
                                          % nowhere in a ladder step
    targetTau = max(abs(tauHuman(dom)));  % fuzzy-ladder target: peak |tau_h|
                                          % over the graded domain
    slot = SU.priSlot;
    axisU = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));

    % ----- Ladder ([S2-5]) ---------------------------------------------
    nB = 1; G = 1;
    chain = '';
    best = [];
    stopReason = '';
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
        spec = injectConfig(sp, nB, G, S.dia);
        R = runRung(spec, points, T, N1, N2, Lchar, axisU, tauHuman, dom, ...
            slot, basePressure, baseXi0, baseXi1, baseXi2, baseWraps, ...
            baseFitn, budget, useParallel, S, nB, G, rung, resultsRoot, ...
            smoke, ladderCsv, targetTau, W);
        evalsUsed = evalsUsed + R.evals;
        chain = sprintf('%s(%d,%d)[%d/%d]', chain, nB, G, ...
            R.meetCells, R.activeCells);   %#ok<AGROW>
        % [S2-3] BEST-MEETING winner: lexicographic on INTEGER domain
        % cell counts; an active tie-break can NEVER record a
        % lower-meeting design.
        if isempty(best) || iBetterMeet(R, best)
            best = R;
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
        S.spec_name, best.nB, best.G, best.xBest(1), best.xBest(2), ...
        best.meetPct, best.activePct, evalsUsed, wallMin, stopReason)
    S.best = best;
    S.evalsUsed = evalsUsed;
    S.wallMin = wallMin;
    S.chain = chain;
    S.stopReason = stopReason;
    S.convergedSite = convergedSite;
    iAppendSurveyRow(surveyCsv, S);
end

fprintf('\nSurvey2 complete. Driver CSVs:\n  %s\n  %s\n', surveyCsv, ladderCsv)

%% =====================================================================
%% Local functions
%% =====================================================================

function R = runRung(spec, points, T, N1, N2, Lchar, axisU, tauHuman, dom, ...
    slot, pres, xi0, xi1, xi2, wraps, fitn, budget, useParallel, S, nB, G, ...
    rung, resDir, smoke, ladderCsv, targetTau, W)
% One ladder rung end to end: sizes [rest; tendon] re-solved at THIS
% configuration on the 2-D graded grid, soft survey2 objective, advisor
% diagnosis, dated mat. No hard constraints anywhere ([S2-2]). Eval
% accounting uses patternsearch's own output.funccount (exact in serial
% AND parallel) plus the client-side box-seed evaluations below.

lb = [0.30 * Lchar; 0.02];
ub = [0.90 * Lchar; max(0.30 * Lchar, 0.06)];
nDom = sum(dom(:));

    function [g, ok] = grade(xq)
    % Grading math ([S2-1]/[S2-2]) -- the 10-22 survey's exact booleans,
    % masked by the surfaces domain; used by BOTH the objective and the
    % final evaluation so the score and the reported percentages can
    % never diverge.
        rest = xq(1);
        tendon = xq(2);
        kmax = 0.75 * rest;
        try
            obj = BiPam_pulley2_20261022(S.spec_name, points, spec.cross, ...
                S.dia, T, rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, ...
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
        % [S2-2] MEET AREA first (1000x), ACTIVE AREA second, dead-cell
        % price at 1e-3 tie-break scale (optimizer only; reported scores
        % stay pure).
        meetFrac = sum(g.meet(dom)) / nDom;
        activeFrac = sum(g.active(dom)) / nDom;
        deadFrac = (sum(g.slack(dom)) + sum(g.infeas(dom))) / nDom;
        f = -(1000 * meetFrac + activeFrac) + 1e-3 * deadFrac;
    end

opts = optimoptions('patternsearch', ...
    'Display', 'none', ...
    'UseParallel', useParallel, ...
    'MaxFunctionEvaluations', budget, ...
    'MeshTolerance', 1e-4, ...
    'StepTolerance', 1e-4);
% 5-point box multi-start (10-22 convention): the all-dead plateau is
% flat even under the dead-cell price, so a single blind x0 can strand
% patternsearch where a taut window exists in bounds.
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
% Sign-convention diagnostic (printed every rung; a systematic 0% flags
% an axis-sign mismatch, not a capability result).
[g0, ok0] = grade(x0);
if ok0
    mSign = dom & g0.taut & abs(tauHuman) > 1;
    if any(mSign(:))
        agree = mean(sign(g0.tauB(mSign)) == sign(tauHuman(mSign)));
        fprintf(['  sign check (initial x): %d taut domain cells with ' ...
            '|tau_h| > 1 N*m, %.0f%% agree with the human sign\n'], ...
            nnz(mSign), 100 * agree)
    end
end
[xBest, fBest, ~, pout] = patternsearch(@surveyObj, x0, [], [], [], [], ...
    lb, ub, [], opts);
evals = pout.funccount + size(seeds, 1);   % + the box seed evaluations

[g, ok] = grade(xBest);
if ~ok
    error('Opt_run_BiPulley_survey2_20261023:finalGrade', ...
        'Final grading failed for %s: %s', S.spec_name, g.reason)
end
meetCells = sum(g.meet(dom));
activeCells = sum(g.active(dom));
meetPct = 100 * meetCells / nDom;
activePct = 100 * activeCells / nDom;
score = -(1000 * meetCells / nDom + activeCells / nDom);

% Unflagged shares ([S2-1]: flagged cells still count in the headline;
% both numbers reported per FIX_SPEC section 2).
du = dom & ~S.surf.flg;
if any(du(:))
    meetPctU = 100 * sum(g.meet(du)) / sum(du(:));
    activePctU = 100 * sum(g.active(du)) / sum(du(:));
else
    meetPctU = NaN; activePctU = NaN;
end
mSk = dom & g.taut & abs(tauHuman) > 1;
if any(mSk(:))
    signAgree = 100 * mean(sign(g.tauB(mSk)) == sign(tauHuman(mSk)));
else
    signAgree = NaN;
end

% Fuzzy ladder diagnosis: the advisor consumes the graded grids
% (taut-cell-only torque test, slack routed to the RoM side) with
% cellWeights = the domain mask ([S2-1] envelope guard).
advStruct = struct('TorqueIns', reshape(abs(g.tauB), N1, N2), ...
    'PulleySlack', g.slack, 'PulleyInfeasible', g.infeas);
[advice, rep] = escalation_advisor_20261022(advStruct, targetTau, ...
    struct('curNBPA', nB, 'curG', G, 'verbose', true, ...
    'cellWeights', W));

fprintf(['  rung %d best: rest %.4f m, tendon %.4f m -> meet %d/%d (%.1f%%), ' ...
    'active %d/%d (%.1f%%), %d evals\n'], rung, xBest, meetCells, nDom, ...
    meetPct, activeCells, nDom, activePct, evals)
fprintf('  rung %d advisor: %s -> steps (+%d,+%d)\n', rung, ...
    advice.verdict, advice.stepNBPA, advice.stepG)

% Dated mat ([S2-8]): new name, never overwritten.
stamp = char(string(datetime('now'), 'yyyyMMdd_HHmmss'));
matFile = fullfile(resDir, sprintf('BiPulley_survey2_%s_nB%d_G%d_%s.mat', ...
    S.spec_name, nB, G, stamp));
surveyMeta = struct('actuator_id', S.id, 'master_name', ...
    S.spec.masterName, 'spec_name', S.spec_name, ...
    'primary_dof', S.surf.priDof, 'secondary_dof', S.surf.secDof, ...
    'bpa_size_mm', S.dia, 'slot', slot, 'nBPA', nB, 'G', G, 'rung', rung, ...
    'smoke', smoke, 'evals', evals, 'meetCells', meetCells, ...
    'activeCells', activeCells, 'domainCells', nDom, ...
    'excludedCells', numel(dom) - nDom, 'meetPct', meetPct, ...
    'activePct', activePct, 'score', score, 'targetTauPeak', targetTau, ...
    'surfaces_csv', 'human_surfaces_5x5_grading_20261023.csv', ...
    'note', S.note, 'created', datetime('now'));
save(matFile, 'surveyMeta', 'spec', 'points', 'T', 'xBest', 'fBest', ...
    'score', 'tauHuman', 'dom', 'g', 'advice', 'rep', 'W');
fprintf('  saved %s\n', matFile)

% Ladder CSV row (master_name quoted: names contain commas)
lFid = fopen(ladderCsv, 'a');
fprintf(lFid, ['%d,"%s",%s,%d,%d,%d,%d,%.6g,%.6g,%d,%d,%d,%.4f,%.4f,' ...
    '%.4f,%d,%d,%d,%.6g,%.6g,%s,%d,%d,%d,%d,%s\n'], S.id, ...
    iCsvSafe(S.spec.masterName), S.spec_name, rung, nB, G, evals, ...
    xBest(1), xBest(2), meetCells, activeCells, nDom, meetPct, ...
    activePct, score, sum(g.slack(dom)), sum(g.infeas(dom)), ...
    sum(g.taut(dom)), rep.minTautTau, targetTau, advice.verdict, ...
    advice.stepNBPA, advice.stepG, advice.converged, smoke, ...
    char(string(datetime('now'), 'yyyyMMdd_HHmmss')));
fclose(lFid);

R = struct('xBest', xBest, 'fBest', fBest, 'score', score, ...
    'meetCells', meetCells, 'activeCells', activeCells, ...
    'meetPct', meetPct, 'activePct', activePct, 'meetPctU', meetPctU, ...
    'activePctU', activePctU, 'signAgree', signAgree, ...
    'evals', evals, 'nB', nB, 'G', G, 'advice', advice, 'rep', rep, ...
    'grade', g, 'matFile', matFile);
end

function better = iBetterMeet(R, best)
% [S2-3] BEST-MEETING comparison, strictly lexicographic on INTEGER
% domain cell counts: meet cells first, active cells second, simpler
% configuration third, earlier rung fourth. An active-percentage
% tie-break can never record a lower-meeting design (integer counts --
% no float-scale hazard at all).
better = R.meetCells > best.meetCells;
if ~better && R.meetCells == best.meetCells
    better = R.activeCells > best.activeCells;
    if ~better && R.activeCells == best.activeCells
        better = (R.nB + R.G) < (best.nB + best.G);
    end
end
end

%% ------------ 2-D graded transforms (both slots swept) ----------------
function [T, N1, N2] = gridT(sp, th1v, th2v)
% [S2-1] 4x4xN1xN2: T(:,:,ii,1) = hinge 1 at th1v(ii) deg; T(:,:,iii,2)
% = hinge 2 at th2v(iii) deg (both hinge angle grids indexed along dim 3,
% dim 4 selects the joint -- the class's T layout contract, established
% 2026-10-22 and reused by sanity_biPulleySpecs27_20261023 gridT_).
% A pinned slot (unmodeled secondary dof, or rigid clone) passes a
% zero-vector -- duplicated neutral transforms, constant robot surface
% along that dimension.
N1 = numel(th1v);
N2 = numel(th2v);
T = zeros(4, 4, N1, N2);
for ii = 1:N1
    T(:, :, ii, 1) = hingeT(sp, 1, th1v(ii));
end
for iii = 1:N2
    T(:, :, iii, 2) = hingeT(sp, 2, th2v(iii));
end
end

function T = hingeT(sp, slot, thetaDeg)
th = thetaDeg * pi / 180;
if strcmpi(strtrim(sp.kin.types{slot}), 'fixed')
    R = eye(3);
    p = [0, 0, 0];
else
    a = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));
    K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
    R = eye(3) + sin(th) * K + (1 - cos(th)) * (K * K);
    p = sp.kin.pivots(slot, :);
    if sp.kin.kneeSlot == slot && isfield(sp.kin, 'kneeFit') && ...
            ~isempty(sp.kin.kneeFit.thetaY)
        thK = sp.kin.kneeFit.thetaY;
        xK = interp1(sp.kin.kneeFit.thetaX, sp.kin.kneeFit.x, thK, ...
            'linear', 'extrap');
        p = [interp1(thK, xK, th, 'linear', 'extrap'), ...
            interp1(thK, sp.kin.kneeFit.y, th, 'linear', 'extrap'), p(3)];
    end
end
T = [R, p(:); 0, 0, 0, 1];
end

%% ------------ Config injection (10-22 injectConfig, size-aware) -------
function spec = injectConfig(spec0, nB, G, dia)
spec = spec0;
spec.diameter = dia;   % [S2-4] both rows same hardware (old-code size)
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

%% ------------ [S2-1] Surfaces loading + per-site dof mapping ----------
function surfaces = loadSurfaces(csvPath)
% Parses human_surfaces_5x5_grading_20261023.csv (header inspected
% 2026-10-23): one row per CELL, columns actuator_id, primary_dof,
% secondary_dof, p_node_deg, s_node_deg, in_graded_domain, tau_h_Nm,
% flag_via/flag_step/flag_slack/flag_short. Returns a containers.Map
% actuator id -> struct with sorted node axes and the 5x5 H / dom / flg
% grids in CANONICAL layout: rows = primary node index, cols = secondary
% node index. The per-site mapper then re-lays H onto the slot dims.
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
need = {'actuator_id', 'primary_dof', 'secondary_dof', 'p_node_deg', ...
    's_node_deg', 'in_graded_domain', 'tau_h_Nm'};
for k = 1:numel(need)
    if ~any(strcmp(raw.Properties.VariableNames, need{k}))
        error('Opt_run_BiPulley_survey2_20261023:surfacesFormat', ...
            'Surfaces CSV lacks column %s', need{k})
    end
end
surfaces = containers.Map('KeyType', 'double', 'ValueType', 'any');
for r = 1:height(raw)
    idn = double(raw.actuator_id(r));
    if isnan(idn)
        continue
    end
    e = struct('priDof', strtrim(raw.primary_dof{r}), ...
        'secDof', strtrim(raw.secondary_dof{r}), ...
        'p', iNum(raw.p_node_deg(r)), 's', iNum(raw.s_node_deg(r)), ...
        'dom', iNum(raw.in_graded_domain(r)) == 1, ...
        'tau', iNum(raw.tau_h_Nm(r)), ...
        'flg', false);
    for f = {'flag_via', 'flag_step', 'flag_slack', 'flag_short'}
        if any(strcmp(raw.Properties.VariableNames, f{1})) && ...
                iNum(raw.(f{1})(r)) == 1
            e.flg = true;
        end
    end
    if isKey(surfaces, idn)
        surfaces(idn) = [surfaces(idn), e];
    else
        surfaces(idn) = e;
    end
end
end

function v = iNum(x)
if isnumeric(x)
    v = double(x);
elseif iscell(x) || isstring(x)
    v = str2double(char(x));
else
    v = double(x);
end
end

function [SU, note] = surfaceForSite(surfaces, sp, dia) %#ok<INUSL>
% Maps one actuator's surface onto the builder's slot layout.
%  - priSlot: the slot whose kin.dofs models the CSV primary (asserted
%    equal to the builder's own kin.primarySlot);
%  - secSlot: the slot modeling the CSV secondary, 0 when none (rigid
%    clone / not in this chain) -- then the robot surface is constant
%    along the secondary axis (remaining joint pinned at 0), stated in
%    the note (the hand-tuned row's own old-code behavior).
% Returns SU with per-slot angle vectors th1/th2 (deg), the human torque
% grid H laid out on (slot-1 dim, slot-2 dim), domain/flg grids, nDom.
SU = struct('priDof', '', 'secDof', '', 'priSlot', 0, 'th1', [], ...
    'th2', [], 'H', [], 'dom', [], 'flg', [], 'nDom', 0);
note = '';
if ~isKey(surfaces, double(sp.actuatorId))
    note = 'no surface';
    return
end
es = surfaces(double(sp.actuatorId));
assert(numel(es) == 25, 'act %d: expected 25 surface rows, got %d', ...
    sp.actuatorId, numel(es))
SU.priDof = es(1).priDof;
SU.secDof = es(1).secDof;

priSlot = find(strcmp(sp.kin.dofs, SU.priDof), 1);
if isempty(priSlot)
    SU.priSlot = 0;
    note = sprintf(['primary dof %s not modeled by the spec hinges ' ...
        '(%s / %s)'], SU.priDof, sp.kin.dofs{1}, sp.kin.dofs{2});
    return
end
assert(priSlot == sp.kin.primarySlot, ...
    ['act %d: CSV primary %s maps to slot %d but builder primarySlot ' ...
     'is %d'], sp.actuatorId, SU.priDof, priSlot, sp.kin.primarySlot)
SU.priSlot = priSlot;
secSlot = find(strcmp(sp.kin.dofs, SU.secDof), 1);
otherSlot = 3 - priSlot;    % the grid dim the CSV secondary maps onto

% Sorted node axes (5 inclusive nodes per axis).
pNodes = sort(unique(round([es.p], 9)));
sNodes = sort(unique(round([es.s], 9)));
assert(numel(pNodes) == 5 && numel(sNodes) == 5, ...
    'act %d: not a 5x5 grid (p %d, s %d nodes)', sp.actuatorId, ...
    numel(pNodes), numel(sNodes))

% Canonical grids: rows = primary node, cols = secondary node.
Hc = nan(5, 5); domC = false(5, 5); flgC = false(5, 5);
for r = 1:25
    ip = find(abs(pNodes - es(r).p) < 1e-9, 1);
    is = find(abs(sNodes - es(r).s) < 1e-9, 1);
    assert(~isempty(ip) && ~isempty(is), 'act %d: node off-grid', ...
        sp.actuatorId)
    assert(isnan(Hc(ip, is)), 'act %d: duplicate cell (%d,%d)', ...
        sp.actuatorId, ip, is)
    Hc(ip, is) = es(r).tau;
    domC(ip, is) = es(r).dom;
    flgC(ip, is) = es(r).flg;
end
assert(all(isfinite(Hc(:))), 'act %d: non-finite tau_h', sp.actuatorId)
assert(any(domC(:)), 'act %d: empty graded domain', sp.actuatorId)

% Re-lay onto slot dims: primary axis on priSlot's dim, CSV secondary
% axis on the other dim.
if priSlot == 1
    H = Hc; dom = domC; flg = flgC;
else
    H = Hc.'; dom = domC.'; flg = flgC.';
end
SU.H = H; SU.dom = dom; SU.flg = flg; SU.nDom = sum(dom(:));

% Per-slot angle vectors: primary nodes on priSlot; the other slot gets
% the secondary nodes when it models them, else zeros (pinned at 0 --
% constant robot surface along that dim, matching the hand-tuned row).
secOn1 = ~isempty(secSlot) && secSlot == 1;
secOn2 = ~isempty(secSlot) && secSlot == 2;
if priSlot == 1
    th1 = pNodes;
    th2 = sNodes * double(secOn2);   % 0-vector when the secondary is
else
    th1 = sNodes * double(secOn1);   % unmodeled (slot pinned at 0)
    th2 = pNodes;
end
SU.th1 = th1;
SU.th2 = th2;
if isempty(secSlot)
    note = sprintf(['secondary %s not a class hinge: remaining joint ' ...
        'pinned at 0, robot surface constant along that axis (same as ' ...
        'the hand-tuned row); %d mm old-code BPA'], SU.secDof, dia);
else
    note = sprintf(['both surface dofs modeled (slots %d+%d); %d mm ' ...
        'old-code BPA'], priSlot, secSlot, dia);
end
end

%% ------------ [S2-4] Old-code diameter (Size.m rule verbatim) ---------
function [dia, why] = oldCodeDiameter(mif)
% Size.m's own thresholds: Fc = 0.9 on the Festo theoretical maxima
% F1 = 630 N (10 mm), F2 = 1500 N (20 mm), F3 = 6000 N (40 mm).
F1 = 630; F2 = 1500; F3 = 6000; Fc = 0.9;
if mif < Fc * F1
    dia = 10;
    why = '';
elseif mif >= Fc * F1 && mif < Fc * F2
    dia = 20;
    why = '';
elseif mif >= Fc * F2 && mif <= Fc * F3
    dia = 40;
    why = '';
else
    dia = NaN;
    why = sprintf('mif %g beyond 40 mm (Use additional muscles)', mif);
end
end

%% ------------ Small helpers ------------------------------------------
function sites = iAppendSite(sites, row)
% Struct-array growth from an EMPTY struct([]) errors in MATLAB; branch.
if isempty(sites)
    sites = row;
else
    sites(end + 1) = row;   %#ok<AGROW>
end
end

function iAppendSurveyRow(surveyCsv, S)
% One checkpoint row per actuator ([S2-8]): appended the moment the site
% reaches a terminal status. built:false rows carry the reason in the
% note (empty or escalated per the ask). Built as a cell array and
% joined so field counting cannot drift; master_name / ladderChain /
% note are double-quoted (master names contain commas).
if ~strcmp(S.status, 'ready') || isempty(S.best)
    if ~isempty(S.surf)
        pri = S.surf.priDof; sec = S.surf.secDof;
        nDomS = sprintf('%d', S.surf.nDom);
        nExS = sprintf('%d', numel(S.surf.dom) - S.surf.nDom);
    else
        pri = ''; sec = ''; nDomS = ''; nExS = '';
    end
    if isnan(S.dia), diaS = ''; else, diaS = sprintf('%d', S.dia); end
    if isempty(S.note), nt = S.reason; else, nt = S.note; end
    c = {sprintf('%d', S.id), ['"', iCsvSafe(S.spec.masterName), '"'], ...
        S.spec_name, 'false', pri, sec, diaS, nDomS, nExS, '0', '0', ...
        '0', '', '', '', '', '', '', '', '', '0', '0.00', '0', ...
        S.status, '', '', ['"', iCsvSafe(nt), '"']};
    fid = fopen(surveyCsv, 'a');
    fprintf(fid, '%s\n', strjoin(c, ','));
    fclose(fid);
    return
end
b = S.best;
if isempty(S.note), nt = ''; else, nt = S.note; end
if isnan(b.meetPctU), mU = ''; else, mU = sprintf('%.4f', b.meetPctU); end
if isnan(b.activePctU), aU = ''; else, aU = sprintf('%.4f', b.activePctU); end
if isnan(b.signAgree), sA = ''; else, sA = sprintf('%.4f', b.signAgree); end
c = {sprintf('%d', S.id), ['"', iCsvSafe(S.spec.masterName), '"'], ...
    S.spec_name, 'true', S.surf.priDof, S.surf.secDof, ...
    sprintf('%d', S.dia), sprintf('%d', S.surf.nDom), ...
    sprintf('%d', numel(S.surf.dom) - S.surf.nDom), ...
    sprintf('%d', sum(S.chain == '(')), ...   % rungs = visited configs
    sprintf('%d', b.nB), sprintf('%d', b.G), ...
    sprintf('%.6g', b.xBest(1)), sprintf('%.6g', b.xBest(2)), ...
    sprintf('%.4f', b.meetPct), sprintf('%.4f', b.activePct), ...
    mU, aU, sA, sprintf('%.4f', b.score), ...
    sprintf('%d', S.evalsUsed), sprintf('%.2f', S.wallMin), ...
    sprintf('%d', S.convergedSite), S.stopReason, ...
    ['"', iCsvSafe(S.chain), '"'], b.matFile, ['"', iCsvSafe(nt), '"']};
fid = fopen(surveyCsv, 'a');
fprintf(fid, '%s\n', strjoin(c, ','));
fclose(fid);
end

function s = iCsvSafe(txt)
s = strrep(txt, '"', '''');
end

function s = iIf(cond, a, b)
if cond
    s = a;
else
    s = b;
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
    mode = 'DRY-RUN (buildability + surfaces verification; no optimization)';
elseif useParallel
    mode = 'SURVEY2 (all sites, parpool <= 8)';
else
    mode = 'SURVEY2 (all sites, serial)';
end
end

%% ------------ Rigid route lengths (10-22 driver, verbatim logic) ------
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
