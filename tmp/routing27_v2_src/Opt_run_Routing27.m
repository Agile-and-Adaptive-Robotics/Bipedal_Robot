% Opt Run Routing27 (v2, penalty restructure — Ben, 2026-10-03)
% Campaign: routing27_v2. Identical machinery and modes to routing27
% (built on the committed BiPulley batch pattern, stock-Gait2392 targets),
% with the constraint/penalty structure Ben specified:
%
%   HARD CONSTRAINTS (inside the TASK-RoM envelope only):
%     - strain floor, pulley feasibility, bone clearance >= inflated BPA
%       radius + 5 mm, lateral footprint cap, muscle-to-muscle spacing.
%     Task envelopes per primary DOF from the functional-task literature
%     (deep squat, jump, sprint, walking/side lunge, grapevine; envelope
%     table and provenance in buildRoutingContext.m).
%
%   PENALTIES (objective, never fatal):
%     - RoM shortfall: mean-square torque deficit over the OUTER poses
%       (outside the task envelope), same normalization and scale as the
%       mean-square envelope term, so partial-RoM solutions are scored.
%     - Anatomy / femoral-neck guard: mean over primary-DOF cells of the
%     - Too-close-to-other-BPAs penalty: soft, squared, threshold at the
%       BPA-diameter comfort distance (the BPAs will squish each other
%       physically; the optimizer should just prefer non-squishing routes).
%       squared excess of the route's off-bone distance beyond 20 mm;
%       routes ballooning far above the femoral neck are far from every
%       bone cloud and pay quadratically (anatW = 1.0).
%     - Pulley-axis internal reaction: unchanged from v1 (5e-2 x rho^2).
%     - Outer-pose violation terms (bone/spacing/strain beyond the
%       envelope) at outerW = 1e2, soft.
%
%   EFFICIENCY (Ben: efficient eval counts):
%     - surrogateopt: 500 evals OR 10 min wall (MaxTime), whichever first.
%     - patternsearch: 1200 evals OR 8 min wall, FunctionTolerance 1e-8.
%     - Search-time pose thinning: the expensive geometry checks
%       (bone clearance, spacing) run on every other envelope pose during
%       the solvers; the FINAL reported evaluation uses the full grid.
%     - Env overrides: ROUTING27_SURROGATE / ROUTING27_PATTERN /
%       ROUTING27_MAXTIME_S / ROUTING27_MAXTIME_P (seconds).
%
%   NOT yet in this version (documented TODO): the improvement term
%   against the 2019 routing baseline (needs the 2019 torque/RoM tables
%   plumbed in); the bowden-cable extensor variant and two-phase Xi
%   scheme belong to the follow-up spec (ROUTING27_V2_SPEC).
%
% MODES (env ROUTING27_MODE): plan / validate / full / smoke, as v1.
% ARTIFACTS (Results\ under this folder): same schema as v1 (CSV columns
% unchanged for comparability); the penalty values are logged to the
% diary per config. Resume via checkpoint, as v1.
% TODO(v2.1): improvement term vs the 2019 routing baseline (needs the
% 2019 per-actuator torque/RoM tables plumbed in; documented, not guessed).

mode = getenv('ROUTING27_MODE');
if isempty(mode)
    mode = 'plan';
end
if strcmp(getenv('ROUTING27_SMOKE'), '1')
    mode = 'smoke';
end
if ~any(strcmp(mode, {'plan', 'validate', 'full', 'smoke'}))
    error('Opt_run_Routing27:mode', 'Unknown ROUTING27_MODE %s', mode)
end

%% Path setup (repo pattern: this folder first, then the Matlab tree)
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:10
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break;
    end
    if strcmp(parent, root) || isempty(parent)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(scriptDir);
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

resDir = fullfile(scriptDir, 'Results');
if ~exist(resDir, 'dir')
    mkdir(resDir)
end
diaryFile = fullfile(resDir, sprintf('routing27_%s.log', mode));
diary(diaryFile);
fprintf('==== ROUTING27 CAMPAIGN v2 (penalty restructure) ====\n');
fprintf('Started %s | mode = %s\n', char(string(datetime('now'))), mode);

%% Map, specs, targets
selected = [];
switch mode
    case 'smoke'
        selected = 12;
    case 'validate'
        selected = [11, 24, 12];
end

mapFile = fullfile(scriptDir, 'actuator_map_27.json');
mapStruct = jsondecode(fileread(mapFile));

[specs, report, chain] = routingSpecsFromRobotbody(mapStruct);

targetsDir = fullfile(scriptDir, 'Human_Torques_27');
if ~exist(fullfile(targetsDir, 'targets27_manifest.json'), 'file')
    error('Opt_run_Routing27:noTargets', ...
        ['Human targets missing. Run gen_gait2392_torque_targets.py ' ...
         'first (D:\\Anaconda\\envs\\opensim\\python.exe).'])
end
for a = 1:numel(mapStruct.actuators)
    if ~isempty(selected) && ...
            ~any(selected == mapStruct.actuators(a).id)
        continue;
    end
    pat = sprintf('targets27_%02d_*_%s_primary.csv', ...
        mapStruct.actuators(a).id, mapStruct.actuators(a).primary_dof);
    if isempty(dir(fullfile(targetsDir, pat)))
        error('Opt_run_Routing27:noTarget', 'Missing target CSV %s', pat)
    end
end

boneDir = fullfile(root, 'Code', 'Matlab', 'Bone_Mesh_Plots', ...
    'Open_Sim_Bone_Geometry');

fid = fopen(fullfile(resDir, 'routing27_spec_report.txt'), 'w');
for r = 1:numel(report)
    fprintf(fid, '[%02d] %s -> %s\n', report(r).actuatorId, ...
        report(r).masterName, report(r).seedChosen);
    for t = 1:numel(report(r).seedTried)
        fprintf(fid, '    seed: %s\n', report(r).seedTried{t});
    end
    fprintf(fid, '    %s\n', report(r).note);
end
fclose(fid);

%% Mode configuration
CONFIGS = struct( ...
    'name', {'single', 'par2', 'par3', 'pulley2', 'pulley2par2'}, ...
    'nBPA', {1, 2, 3, 1, 2}, ...
    'gain', {1, 1, 1, 2, 2});

switch mode
    case 'plan'
        selected = [];
        validateMap = [];
        budgetS = 0;
        budgetP = 0;
        useParallel = false;
    case 'smoke'
        selected = 12;
        validateMap = struct('ids', {12}, 'configs', {{'pulley2'}});
        budgetS = 0;
        budgetP = 40;
        useParallel = false;
    case 'validate'
        selected = [11, 24, 12];
        validateMap = struct('ids', {11, 24, 12}, ...
            'configs', {{'single'}, {'single'}, {'pulley2'}});
        budgetS = 60;
        budgetP = 200;
        useParallel = false;
    case 'full'
        selected = [];
        validateMap = [];
        budgetS = 500;      % v2: 400 -> 500 (surrogate is the global phase)
        budgetP = 1200;     % v2: 2000 -> 1200 (refinement, time-capped too)
        useParallel = true;
end
maxTimeS = 600;             % surrogateopt wall cap, seconds
maxTimeP = 480;             % patternsearch wall cap, seconds
if strcmp(mode, 'validate') || strcmp(mode, 'smoke')
    maxTimeS = 120;
    maxTimeP = 120;
end
envS = str2double(getenv('ROUTING27_SURROGATE'));
envP = str2double(getenv('ROUTING27_PATTERN'));
envTS = str2double(getenv('ROUTING27_MAXTIME_S'));
envTP = str2double(getenv('ROUTING27_MAXTIME_P'));
if ~isnan(envS)
    budgetS = envS;
end
if ~isnan(envP)
    budgetP = envP;
end
if ~isnan(envTS)
    maxTimeS = envTS;
end
if ~isnan(envTP)
    maxTimeP = envTP;
end

Nangles = 13;

%% Plan print
fprintf('\n---- PLAN (v2: task-envelope hard set + penalty outer set) ----\n');
fprintf('Grid %d primary angles | budgets surrogate %d (max %ds) + pattern %d (max %ds) | parallel %d\n', ...
    Nangles, budgetS, maxTimeS, budgetP, maxTimeP, useParallel);
fprintf('%-4s %-32s %-12s %-18s %-16s %s\n', 'id', 'actuator', 'seed', ...
    'primary DOF', 'task range', 'configs');
for s = 1:numel(specs)
    sp = specs(s);
    if ~isempty(selected) && ~any(selected == sp.actuatorId)
        continue;
    end
    cfgNames = configNamesFor(sp.actuatorId, selected, validateMap, CONFIGS);
    trTxt = '(plan defers to ctx)';
    fprintf('%-4d %-32s %-12s %-18s %-16s %s\n', sp.actuatorId, sp.masterName, ...
        sp.name, sp.primaryDof, trTxt, strjoin(cfgNames, '/'));
end
fprintf('--------------\n\n');

if strcmp(mode, 'plan')
    fprintf('PLAN complete. Set ROUTING27_MODE=validate|full to execute.\n');
    diary off;
    return;
end

%% Checkpoint (resume support; per-mode file)
ckptFile = fullfile(resDir, sprintf('routing27_checkpoint_%s.mat', mode));
doneKeys = {};
finalizedRoutes = struct('name', {}, 'Spelvis', {}, 'latExtent', {}, ...
    'config', {}, 'margin', {}, 'actuatorId', {});
if isfile(ckptFile)
    Sck = load(ckptFile);
    doneKeys = Sck.doneKeys;
    if isfield(Sck, 'finalizedRoutes')
        finalizedRoutes = Sck.finalizedRoutes;
    end
    fprintf('Checkpoint: %d config(s) done, %d route(s) finalized.\n', ...
        numel(doneKeys), numel(finalizedRoutes));
end

summaryFile = fullfile(resDir, sprintf('routing27_summary_%s.csv', mode));
if ~isfile(summaryFile)
    fid = fopen(summaryFile, 'a');
    fprintf(fid, ['actuatorId,masterName,seed,config,nBPA,gain,exitflagS,' ...
        'exitflagP,objective,worstC,marginWorstRel,minBoneMM,' ...
        'minSpacingMM,latExtentMM,rho,feasible,humanPeakNm,' ...
        'tauRobotPeakNm,chosen,timestamp\n']);
    fclose(fid);
end

if useParallel && isempty(gcp('nocreate'))
    parpool(10);    % easteregg2 (10 cores, 128 GB)
end

%% Contexts for the selected actuators (finalized list refreshes as
%% routes finalize: contexts are rebuilt from the checkpointed list)
t0 = tic;
for s = 1:numel(specs)
    sp = specs(s);
    if ~isempty(selected) && ~any(selected == sp.actuatorId)
        continue;
    end
    cfgIdx = configIndicesFor(sp.actuatorId, selected, validateMap, CONFIGS);
    if ~any(~isdone(doneKeys, sp.actuatorId, cfgIdx, CONFIGS))
        continue;
    end
    fprintf('\n===== [%02d] %s (seed %s; %s) =====\n', sp.actuatorId, ...
        sp.masterName, sp.name, sp.primaryDof);
    ctx = buildRoutingContext(sp, chain, targetsDir, boneDir, ...
        finalizedRoutes, struct('Nangles', Nangles));
    fprintf('ctx %02d: %d pts, cross [%d %d], %s at crossing %d, Lchar %.3f m, peak human %.1f N*m\n', ...
        sp.actuatorId, size(sp.points, 1), sp.cross(1), sp.cross(2), ...
        sp.primaryDof, ctx.targetCrossing, ctx.Lchar, ctx.humanPeakNm);
    fprintf('task envelope [%.0f, %.0f] deg: %d of %d poses hard\n', ...
        ctx.taskRangeDeg(1), ctx.taskRangeDeg(2), ...
        numel(ctx.envCells), numel(ctx.thetaPrimaryDeg));

    best = struct('f', inf, 'res', [], 'cfgName', '');
    for ci = cfgIdx
        cfg = CONFIGS(ci);
        key = sprintf('%02d_%s', sp.actuatorId, cfg.name);
        if any(strcmp(doneKeys, key))
            fprintf('  %s: already done (checkpoint), skipping.\n', cfg.name);
            continue;
        end
        fprintf('  --- config %s (nBPA=%d, G=%d) at %s ---\n', cfg.name, ...
            cfg.nBPA, cfg.gain, char(string(datetime('now'))));
        try
            res = runConfig(ctx, cfg, budgetS, budgetP, useParallel, ...
                maxTimeS, maxTimeP);
            doneKeys{end + 1} = key; %#ok<AGROW>
            save(ckptFile, 'doneKeys', 'finalizedRoutes');
            writeSummary(summaryFile, sp, cfg, res, false);
            feas = res.worstC <= 1e-6;
            if feas && res.f < best.f
                best = struct('f', res.f, 'res', res, 'cfgName', cfg.name);
            elseif isempty(best.res) && ~feas && strcmp(best.cfgName, '')
                best = struct('f', res.f, 'res', res, 'cfgName', cfg.name);
            end
        catch err
            fprintf('  CONFIG FAILED (%s): %s\n%s\n', cfg.name, err.message, ...
                getReport(err, 'extended', 'hyperlinks', 'off'));
            writeSummaryError(summaryFile, sp, cfg);
        end
    end

    if isempty(best.res)
        fprintf('  no new feasible-or-not result this pass; route not (re)finalized.\n');
        continue;
    end

    fprintf('  chosen: %s | objective %.4g | worstC %.3g | marginRel %+.4f | Jrom %.3g | Janat %.3g\n', ...
        best.cfgName, best.res.f, best.res.worstC, best.res.marginRel, ...
        best.res.Jrom, best.res.Janat);
    fr = struct();
    fr.name = sprintf('%02d_%s', sp.actuatorId, sp.masterName);
    fr.Spelvis = best.res.Spelvis;
    fr.latExtent = best.res.latExtent;
    fr.config = best.cfgName;
    fr.margin = best.res.marginRel;
    fr.actuatorId = sp.actuatorId;
    finalizedRoutes(end + 1) = fr; %#ok<AGROW>
    save(ckptFile, 'doneKeys', 'finalizedRoutes');
    ci = find(strcmp({CONFIGS.name}, best.cfgName), 1);
    writeSummary(summaryFile, sp, CONFIGS(ci), best.res, true);

    % Dated full-workspace capture (bare save behind liveRun, Ben
    % directive; the flag is cleared immediately before the save).
    liveRun = true;
    stamp = char(string(datetime('now'), 'yyyyMMdd_HHmm'));
    outFile = fullfile(resDir, sprintf('Routing27_%s_%s.mat', ...
        ctx.label, stamp));
    clear liveRun
    save(outFile);
    fprintf('  saved %s\n', outFile);
end

fprintf('\nCAMPAIGN %s done: %d config(s) done, %d route(s) finalized, %.1f min elapsed.\n', ...
    mode, numel(doneKeys), numel(finalizedRoutes), toc(t0) / 60);
fprintf('Summary: %s\n', summaryFile);
diary off;

%% =====================================================================
%% One actuator x one config (v2: envelope-hard / outer-soft evaluate)
%% =====================================================================
function res = runConfig(ctx, cfg, budgetS, budgetP, useParallel, maxTimeS, maxTimeP)
spec = ctx.spec;
k = ctx.targetCrossing;
pul = struct('nPulleyBPA', cfg.nBPA, ...
    'tackleLineParts', [1, cfg.gain], ...    % one tackle, distal crossing
    'routingMode', 'moving_via');
nFin = numel(ctx.finalizedSegs);
nA = numel(ctx.thetaPrimaryDeg);
nSeg0 = size(spec.points, 1) - 1;

    function points = applyDeltas(x)
    x = x(:);                          % surrogateopt passes rows
    points = spec.points;
    points(1, :) = points(1, :) + x(1:3)';
    points(end, :) = points(end, :) + x(4:6)';
    for j = 1:numel(ctx.movableViaRows)
        r = ctx.movableViaRows(j);
        points(r, :) = points(r, :) + x(6 + 3 * (j - 1) + (1:3))';
    end
    end

    function [st, ok] = evaluate(x, thin)
    if nargin < 2
        thin = true;
    end
    x = x(:);
    points = applyDeltas(x);
    rest = x(ctx.iRest);
    tendon = x(ctx.iTend);
    kmax = 0.75 * rest;             % 25 percent maximum contraction (KMAX)
    try
        obj = BiPam_pulley(spec.name, points, spec.cross, ctx.diameter, ...
            ctx.T, rest, kmax, tendon, ctx.fitting, ctx.pressure, ...
            ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.wraps, pul);
    catch err
        ok = false; %#ok<NASGU>
        st = struct('f', 1e12, 'worstC', 1e3, 'marginRel', -1e3, ...
            'minBone', 0, 'minSpacing', inf, 'latExtent', 1e3, ...
            'rho', 0, 'Jrom', 0, 'Janat', 0, ...
            'tauEff', zeros(nA, 1), 'Sp', []);
        return
    end

    % Joint torque on the primary crossing: the insertion moment
    % projected on the hinge axis (invariant across a hinge's frames).
    tq3 = obj.Torque_ins(:, :, :, k);           % N1 x 3 x N2
    [n1, ~, n2] = size(tq3);
    tauGrid = zeros(n1, n2);
    for i1 = 1:n1
        for i2 = 1:n2
            tauGrid(i1, i2) = dot(reshape(tq3(i1, :, i2), 1, 3), ctx.axisHat);
        end
    end
    strainGrid = reshape(obj.strain_p, n1, n2);
    if ctx.primaryDim == 3
        tau = abs(tauGrid(:, 1));
        infeas = obj.PulleyInfeasible(:, 1);
        strainV = strainGrid(:, 1);
    else
        tau = abs(tauGrid(1, :))';
        infeas = obj.PulleyInfeasible(1, :)';
        strainV = strainGrid(1, :)';
    end

    envCells = ctx.envCells;
    outerCells = ctx.outerCells;
    % Search-time thinning of the expensive geometry poses (full grid at
    % the final evaluation): every other envelope and outer cell.
    if thin
        if numel(envCells) > 6
            envCells = envCells(1:2:end);
        end
        if numel(outerCells) > 4
            outerCells = outerCells(1:2:end);
        end
    end

    % ---- Hard set (envelope poses only) ----
    c1 = obj.strainFloor() - min(strainV(envCells));
    c2 = double(sum(infeas(envCells)));

    % Torque meet-or-beat over the FULL envelope (cheap: grid is done).
    required = (1 + ctx.requiredMargin) * ctx.humanAbsGrid;
    tauEff = tau;
    tauEff(isnan(tauEff)) = 0;
    deficit = max(0, required - tauEff) / ctx.torqueScale;
    Jworst = max(deficit(ctx.envCells));
    Jtorque = mean(deficit(ctx.envCells) .^ 2);
    overshoot = max(0, tauEff(ctx.envCells) - ctx.humanAbsGrid(ctx.envCells));
    Jover = 1e-3 * mean((overshoot / ctx.torqueScale) .^ 2);

    % ---- PENALTY: RoM shortfall outside the task envelope ----
    if ~isempty(ctx.outerCells)
        Jrom = ctx.romW * mean(deficit(ctx.outerCells) .^ 2);
    else
        Jrom = 0;
    end

    % ---- Pulley-axis reaction penalty (v1 convention kept) ----
    rho = 0;
    Jreact = 0;
    if cfg.gain > 1
        Rv = obj.ReactionFmag(~isnan(obj.ReactionFmag));
        if ~isempty(Rv)
            rho = max(Rv) / max(cfg.nBPA * max(obj.Fmax, 1), 1);
            Jreact = 5e-2 * rho .^ 2;
        end
    end

    % Rigid-route geometry: segments in frame 1 per primary cell.
    S = routingChainLib('routeSegments', points, spec.cross, ctx.T);
    nRows = size(S, 1);
    cellOf = zeros(nRows, 1);
    for r = 1:nRows
        lin = floor((r - 1) / nSeg0) + 1;
        if ctx.primaryDim == 3
            cellOf(r) = mod(lin - 1, size(ctx.T, 3)) + 1;
        else
            cellOf(r) = floor((lin - 1) / size(ctx.T, 3)) + 1;
        end
    end

    % ---- Bone clearance per evaluated cell + anatomy penalty ----
    % Hard: envelope cells. Anatomy: mean squared excess of per-cell
    % off-bone distance beyond anatD0 over ALL evaluated cells (the
    % femoral-neck guard: a route ballooning above the neck stands far
    % from every cloud and pays quadratically). Outer violations soft.
    geoCells = [envCells, outerCells];
    dminCell = inf(1, numel(geoCells));
    for gi = 1:numel(geoCells)
        a = geoCells(gi);
        SA = S(cellOf == a, :);
        samp = [SA(:, 1:3); SA(:, 4:6); (SA(:, 1:3) + SA(:, 4:6)) / 2];
        dminCell(gi) = routingChainLib('cloudMinDist', samp, ctx.cellClouds{a});
    end
    nEnvG = numel(envCells);
    dEnv = dminCell(1:nEnvG);
    dOut = dminCell(nEnvG + 1:end);
    minBone = min(dminCell);
    c3 = ctx.boneRequiredClear - min(dEnv);
    excess = max(0, dminCell - ctx.anatD0) / ctx.anatD0;
    Janat = ctx.anatW * mean(excess .^ 2);
    if ~isempty(dOut)
        JboneOut = ctx.outerW * mean((max(0, ...
            ctx.boneRequiredClear - dOut) / ctx.boneRequiredClear) .^ 2);
    else
        JboneOut = 0;
    end

    % ---- Lateral footprint (hard cap at envelope extent; the gradient
    % penalty uses the full extent so ballooning is always discouraged) ----
    mids = (S(:, 1:3) + S(:, 4:6)) / 2;
    zAll = [S(:, 3); mids(:, 3)];
    latExtent = max(abs(zAll)) + ctx.bpaRadius ...
        + ((cfg.nBPA - 1) * ctx.bundleSpacing) / 2;
    c4 = latExtent - ctx.latCap;
    Jfoot = max(0, latExtent / ctx.latCap - 1) .^ 2;

    % ---- Spacing vs finalized routes (envelope hard, outer soft) ----
    envSegMask = ismember(cellOf, envCells);
    SpE = xformSegs(S(envSegMask, :), ctx.TtoPelvis);
    minSpacing = inf;
    for r = 1:nFin
        Sf = ctx.finalizedSegs{r};
        for a = 1:size(SpE, 1)
            for b = 1:size(Sf, 1)
                dd = routingChainLib('segSegDist', SpE(a, 1:3), ...
                    SpE(a, 4:6), Sf(b, 1:3), Sf(b, 4:6));
                minSpacing = min(minSpacing, dd);
            end
        end
    end
    c5 = ctx.spacingRequired - minSpacing;
    outerSegMask = ismember(cellOf, outerCells);
    if any(outerSegMask) && nFin > 0
        SpO = xformSegs(S(outerSegMask, :), ctx.TtoPelvis);
        minSpacingOut = inf;
        for r = 1:nFin
            Sf = ctx.finalizedSegs{r};
            for a = 1:size(SpO, 1)
                for b = 1:size(Sf, 1)
                    dd = routingChainLib('segSegDist', SpO(a, 1:3), ...
                        SpO(a, 4:6), Sf(b, 1:3), Sf(b, 4:6));
                    minSpacingOut = min(minSpacingOut, dd);
                end
            end
        end
        JspaceOut = ctx.outerW * (max(0, ...
            ctx.spacingRequired - minSpacingOut) / ctx.spacingRequired) .^ 2;
    else
        JspaceOut = 0;
    end

    % Too-close-to-other-BPAs soft penalty (Ben 2026-10-03: 'penalize if the
    % routes are too close'; the BPAs will squish each other physically, but
    % the optimizer should prefer routes that don't force squishing).
    JspaceClose = 0;
    if nFin > 0 && minSpacing < ctx.spacingComfort
        JspaceClose = 1e-1 * ((ctx.spacingComfort - minSpacing) / ...
            ctx.spacingComfort) .^ 2;
    end

    % Outer strain-floor violations, soft.
    if ~isempty(ctx.outerCells)
        so = max(0, obj.strainFloor() - min([strainV(ctx.outerCells); inf]));
        JstrainOut = ctx.outerW * (so / max(abs(obj.strainFloor()), 1)) .^ 2;
    else
        JstrainOut = 0;
    end

    worstC = max([c1, c2, c3, c4, c5, 0]);

    % Design-change penalty (moved-point deltas + rest/tendon drift).
    dx = x(1:ctx.iRest - 1) - ctx.x0(1:ctx.iRest - 1);
    J = 1e5 * Jworst + 1e3 * Jtorque + Jover + Jreact + Jfoot + JspaceClose ...
        + Jrom + Janat + JboneOut + JspaceOut + JstrainOut ...
        + 1e-2 * sum(dx .^ 2) / (0.060 .^ 2) ...
        + 1e-3 * ((x(ctx.iRest) - ctx.x0(ctx.iRest)) / 0.04) .^ 2 ...
        + 1e-3 * ((x(ctx.iTend) - ctx.x0(ctx.iTend)) / 0.04) .^ 2 ...
        + 10 * max(0, worstC);
    if ~isfinite(J)
        J = 1e12;
    end

    nz = required(ctx.envCells) > 1e-9;
    if any(nz)
        marginRel = min((tauEff(ctx.envCells(nz)) - required(ctx.envCells(nz))) ...
            ./ required(ctx.envCells(nz)));
    else
        marginRel = 0;
    end

    st = struct('f', J, 'worstC', worstC, 'marginRel', marginRel, ...
        'minBone', minBone, 'minSpacing', minSpacing, ...
        'latExtent', latExtent, 'rho', rho, ...
        'Jrom', Jrom, 'Janat', Janat, 'JspaceClose', JspaceClose, ...
        'tauEff', tauEff, 'Sp', SpE);
    ok = true;
    end

    function f = objF(xq)
    [st, okk] = evaluate(xq, true);
    if ~okk
        f = 1e12;
    else
        f = st.f;
        if ~isscalar(f)
            error('objF:notScalar', 'objective has size %s', ...
                mat2str(size(f)));
        end
    end
    end

    function [c, ceq] = nonlconF(xq)
    [st, okk] = evaluate(xq, true);
    if ~okk
        c = 1e3;
    else
        c = st.worstC;
    end
    ceq = [];
    end

    function [fv, c, ceq] = objconstrF(xq)
    fv = objF(xq);
    c = nonlconF(xq);
    ceq = [];
    end

st0 = evaluate(ctx.x0, false);
fprintf('  initial: f = %.4g, worstC = %.3g, marginRel = %+.4f, Jrom %.3g, Janat %.3g\n', ...
    st0.f, st0.worstC, st0.marginRel, st0.Jrom, st0.Janat);

x_best = ctx.x0;
f_best = st0.f;
exitflagS = NaN;
exitflagP = NaN;

if budgetS > 0
    optsS = optimoptions('surrogateopt', 'Display', 'iter', ...
        'UseParallel', useParallel, 'MaxFunctionEvaluations', budgetS, ...
        'MaxTime', maxTimeS, 'ConstraintTolerance', 1e-6);
    [xS, fS, exitflagS] = surrogateopt(@objconstrF, ctx.lb, ctx.ub, optsS);
    if fS < f_best
        x_best = xS;
        f_best = fS;
    end
end
optsP = optimoptions('patternsearch', 'Display', 'iter', ...
    'UseParallel', useParallel, 'MaxFunctionEvaluations', budgetP, ...
    'MaxTime', maxTimeP, 'MeshTolerance', 1e-4, 'StepTolerance', 1e-4, ...
    'FunctionTolerance', 1e-8, 'ConstraintTolerance', 1e-6);
[xB, fB, exitflagP] = patternsearch(@objF, x_best, [], [], [], [], ...
    ctx.lb, ctx.ub, @nonlconF, optsP);
if fB < f_best
    x_best = xB;
    f_best = fB;
end

stF = evaluate(x_best, false);    % full grid: the honest reported numbers
fprintf('  final: f = %.4g, worstC = %.3g, marginRel = %+.4f, minBone %.1f mm, lat %.1f mm, Jrom %.3g, Janat %.3g\n', ...
    stF.f, stF.worstC, stF.marginRel, 1000 * stF.minBone, ...
    1000 * stF.latExtent, stF.Jrom, stF.Janat);

res = struct('x', x_best, 'f', f_best, 'worstC', stF.worstC, ...
    'marginRel', stF.marginRel, 'minBone', stF.minBone, ...
    'minSpacing', stF.minSpacing, 'latExtent', stF.latExtent, ...
    'rho', stF.rho, 'Jrom', stF.Jrom, 'Janat', stF.Janat, 'JspaceClose', stF.JspaceClose, ...
    'tau', stF.tauEff, 'Spelvis', stF.Sp, ...
    'nBPA', cfg.nBPA, 'gain', cfg.gain, 'exitflagS', exitflagS, ...
    'exitflagP', exitflagP, 'humanPeakNm', ctx.humanPeakNm, ...
    'tauRobotPeakNm', max(stF.tauEff));
end

%% ---------------------------------------------------------------------
function Sp = xformSegs(S, T)
if isempty(S)
    return
end
Sp = S;
Sp(:, 1:3) = (T(1:3, 1:3) * S(:, 1:3).' + T(1:3, 4)).';
Sp(:, 4:6) = (T(1:3, 1:3) * S(:, 4:6).' + T(1:3, 4)).';
end

function names = configNamesFor(id, selected, validateMap, CONFIGS)
idx = configIndicesFor(id, selected, validateMap, CONFIGS);
names = {CONFIGS(idx).name};
end

function idx = configIndicesFor(id, selected, validateMap, CONFIGS) %#ok<INUSL>
if isstruct(validateMap) && ~isempty(validateMap)
    hit = find([validateMap.ids] == id, 1);
    if ~isempty(hit)
        idx = find(ismember({CONFIGS.name}, validateMap(hit).configs{1}));
        return;
    end
end
idx = 1:numel(CONFIGS);
end

function tf = isdone(doneKeys, id, cfgIdx, CONFIGS)
tf = false(1, numel(cfgIdx));
for j = 1:numel(cfgIdx)
    tf(j) = any(strcmp(doneKeys, sprintf('%02d_%s', id, CONFIGS(cfgIdx(j)).name)));
end
end

function writeSummary(summaryFile, sp, cfg, res, chosen)
fid = fopen(summaryFile, 'a');
fprintf(fid, ['%d,%s,%s,%s,%d,%d,%g,%g,%.6g,%.6g,%.6f,%.3f,%.3f,' ...
    '%.3f,%.4f,%d,%.2f,%.2f,%d,%s\n'], ...
    sp.actuatorId, sp.masterName, sp.name, cfg.name, cfg.nBPA, cfg.gain, ...
    res.exitflagS, res.exitflagP, res.f, res.worstC, res.marginRel, ...
    1000 * res.minBone, 1000 * min(res.minSpacing, 999), ...
    1000 * res.latExtent, res.rho, (res.worstC <= 1e-6), ...
    res.humanPeakNm, res.tauRobotPeakNm, chosen, ...
    char(string(datetime('now'), 'yyyyMMdd_HHmm:ss')));
fclose(fid);
end

function writeSummaryError(summaryFile, sp, cfg)
fid = fopen(summaryFile, 'a');
fprintf(fid, ['%d,%s,%s,%s,%d,%d,error,error,nan,nan,nan,nan,nan,' ...
    'nan,nan,0,nan,nan,0,%s\n'], ...
    sp.actuatorId, sp.masterName, sp.name, cfg.name, cfg.nBPA, cfg.gain, ...
    char(string(datetime('now'), 'yyyyMMdd_HHmm:ss')));
fclose(fid);
end
