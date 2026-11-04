function fr_launch_run1()
%FR_LAUNCH_RUN1 run_1_corrected campaign driver (Ben 2026-11-04 directive).
%
% Corrected formulation relative to run_0:
%   - SIGNED per-pose human torque targets (targets_v2 via fr_targets v2;
%     force-length modulated, sign-verified) - Ben item 2;
%   - Xi3 per-pose inside the equilibrium (fr_mechBuild v2 + derived
%     class, fr_xi3_check PASSED required) - Ben item 3;
%   - sizing winners overrides (fr_sizingSearch output; contract counts
%     were centers, never mandatory) - Ben item 4;
%   - packaging: rigid hardware hard, BPA corridors soft (fr_obstacles v2)
%     and the packaging_v2 obstacle model when registered - Ben items 5/6;
%   - corrected feasibility labeling: 'solved' only when every hard
%     constraint is within its stated tolerance; 'best-infeasible'
%     otherwise, top violator named - Ben item 7;
%   - warm starts from run_0 checkpoints (STARTING POINTS, never results);
%   - SERIAL compliant smoke on one group BEFORE any parfor (the run_0
%     crash guard: all 10 workers died inside their first compliant
%     evaluations on 2026-11-04 ~06:49).
%
% Launched detached:
%   "D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd( ...
%    '.../Functional_Routing_20261027'); addpath('campaign_20261103'); ...
%    fr_launch_run1"
% A run_1/STOP file stops cleanly at the next layer/round boundary.

C = fr_setup();
runDir = fullfile(C.camp, 'run_1');
ckDir = fullfile(runDir, 'checkpoints');
logDir = fullfile(runDir, 'logs');
if ~exist(ckDir, 'dir'); mkdir(ckDir); end
if ~exist(logDir, 'dir'); mkdir(logDir); end
masterLog = fullfile(logDir, 'master.log');
diary(masterLog);
cleanupObj = onCleanup(@() diary('off'));

tAll = tic;
startedAt = datetime('now');
fprintf('=== fr_launch_run1 (corrected formulation) | %s | machine %s ===\n', ...
    startedAt, getenv('COMPUTERNAME'));

% ---- PRECONDITIONS (fail loudly, before any compute) ----------------------
tgtMat = fullfile(C.camp, 'targets_v2', 'signed_targets.mat');
assert(exist(tgtMat, 'file') == 2, ...
    'fr_launch_run1:precondition', ...
    'targets_v2/signed_targets.mat missing - signed targets required');
xiLog = fullfile(C.camp, 'tmp', 'xi3_check_final.log');
assert(exist(xiLog, 'file') == 2 && ~isempty(regexpi( ...
    fileread(xiLog), 'PASSED', 'once')), ...
    'fr_launch_run1:precondition', ...
    'fr_xi3_check.log with PASSED missing - Xi3 per-pose required');
sizingFile = fullfile(C.camp, 'sizing_winners.json');
haveSizing = exist(sizingFile, 'file') == 2;
if haveSizing
    sz = jsondecode(fileread(sizingFile));
    overrides = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if iscell(sz.groups)
        groupsList = sz.groups;
    else
        groupsList = arrayfun(@(g) g, sz.groups, 'UniformOutput', false);
    end
    for i = 1:numel(groupsList)
        g = groupsList{i};
        if ~isfield(g, 'diaStr'); continue; end
        ov = struct('diaStr', g.diaStr, 'nBPA', g.nBPA, ...
            'family', g.family, 'gain', g.gain);
        overrides(g.name) = ov;
    end
    fprintf('sizing winners: %d group config(s) applied\n', overrides.Count);
else
    overrides = containers.Map('KeyType', 'char', 'ValueType', 'any');
    fprintf(['WARNING: sizing_winners.json missing - running on the ' ...
        'contract seed configs (documented deviation)\n']);
end

rng(C.rngSeed);
kin = fr_kinematics(C.osimStock);
CASES = fr_cases27(C, kin, overrides);
nG = numel(CASES);
byId = containers.Map('KeyType', 'double', 'ValueType', 'double');
byName = containers.Map('KeyType', 'char', 'ValueType', 'double');
for c = 1:nG
    byId(CASES(c).groupMapId) = c;
    byName(CASES(c).name) = c;
end

tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
[layerIds, regionName] = fr_layerIds(tc, CASES);
allIds = [];
for r = 1:numel(layerIds)
    allIds = [allIds, layerIds{r}]; %#ok<AGROW>
end
fprintf('layer order: %s\n', mat2str(allIds));

% ---- run_0 warm starts (STARTING POINTS, never results) -------------------
% v2: pack the checkpoint's DESIGN STRUCT into the x-vector the solver
% expects (gate r1 of run_1: struct warm starts crashed every warm-started
% (group, branch) with MATLAB:badsubscript in fr_unpackX)
wRigid = containers.Map('KeyType', 'char', 'ValueType', 'any');
wComp = containers.Map('KeyType', 'char', 'ValueType', 'any');
ck0 = fullfile(C.camp, 'run_0', 'checkpoints');
d0 = dir(fullfile(ck0, '*.mat'));
for i = 1:numel(d0)
    S = load(fullfile(ck0, d0(i).name));
    if ~isfield(S, 'csName') || ~isfield(S, 'design'); continue; end
    if ~byName.isKey(S.csName); continue; end
    cs0 = CASES(byName(S.csName));
    if S.nX ~= cs0.nX; continue; end   % dimension guard for x0Override
    % 10 mm guard (run_1 resume evidence): every run_0 checkpoint predates
    % the sizing search, i.e. it is a 20 mm-era design. Warm-starting a
    % 10 mm winner from a 20 mm design poisons the solve (grac/sar came
    % out 100% tx-infeasible: the 10 mm sleeve cannot close from a 20 mm
    % length world). 10 mm groups start from their SEED instead.
    if overrides.isKey(S.csName)
        ovw = overrides(S.csName);
        if isfield(ovw, 'diaStr') && strcmpi(ovw.diaStr, '10')
            continue
        end
    end
    if S.branch(1) == 'r'
        wRigid(S.csName) = fr_packX(cs0, S.design);
    else
        wComp(S.csName) = fr_packX(cs0, S.design);
    end
end
fprintf('warm starts from run_0: %d rigid, %d compliant\n', ...
    wRigid.Count, wComp.Count);

% ---- SERIAL compliant smoke (run_0 crash guard) ---------------------------
fprintf('\n--- COMPLIANT SMOKE (serial, psoas, budget 100) ---\n');
csSmoke = CASES(byName('psoas'));
smokeOpts = struct('quiet', true, 'budget', 100, 'envelopes', {{}});
try
    [x0Smoke, ~, infoSmoke] = fr_solveCase(C, csSmoke, kin, ...
        'compliant', smokeOpts);
    if isempty(x0Smoke)
        error('fr_launch_run1:smokeNoDesign', ...
            'compliant smoke returned no design');
    end
    fprintf('compliant smoke PASSED (f %.3f) - parallel solves unlocked\n', ...
        infoSmoke.fBest);
catch smokeErr
    fprintf(2, 'COMPLIANT SMOKE FAILED: %s | %s\n', ...
        smokeErr.identifier, smokeErr.message);
    fr_writeStatus(runDir, startedAt, true, [], [], containers.Map( ...
        'KeyType', 'double', 'ValueType', 'any'));
    error('fr_launch_run1:smokeFailed', ...
        'serial compliant smoke failed - no parallel campaign (run_0 crash guard)');
end

% ---- manifest + resume scan ------------------------------------------------
resArr = fr_scanCheckpoints(ckDir, CASES);
fprintf(['resume scan (run_1): %d/%d rigid, %d/%d compliant, %d/%d ' ...
    'refine found\n'], sum([resArr.hasRigid]), nG, ...
    sum([resArr.hasCompliant]), nG, sum([resArr.hasRefine]), nG);
manFile = fullfile(runDir, 'manifest.json');
if ~exist(manFile, 'file')
    fr_writeManifest(manFile, C, tc, CASES, layerIds, regionName);
    fprintf('manifest written: %s\n', manFile);
end

% ---- envelope-of-record (none at fresh launch; resume-aware) ---------------
envRecord = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(resArr)
    ra = resArr(i);
    if ra.hasRigid
        S = load(ra.ckRigid);
        envRecord(S.mapId) = S.envelope;
    elseif ra.hasCompliant
        S = load(ra.ckCompliant);
        envRecord(S.mapId) = S.envelope;
    end
end
fprintf('envelope-of-record restored for %d placed group(s)\n', ...
    envRecord.Count);

% ---- pool: env vars BEFORE parpool, capped at 10 ---------------------------
setenv('OMP_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS', '1');
if isempty(gcp('nocreate'))
    pool = parpool('Processes', C.parallelMaxWorkers);
    fprintf('parpool: %d workers\n', pool.NumWorkers);
end

stopFile = fullfile(runDir, 'STOP');

% ---- PASS 1: sequential layers, parfor within layer ------------------------
layerResults = cell(1, numel(layerIds));
solvedEnv = containers.Map('KeyType', 'double', 'ValueType', 'any');
for L = 1:numel(layerIds)
    if exist(stopFile, 'file') == 2
        fprintf(['STOP detected before layer %d (%s) - stopping ' ...
            'cleanly; checkpoints intact and resumable.\n'], L, ...
            regionName{L});
        fr_writeStatus(runDir, startedAt, true, layerResults, resArr, ...
            envRecord);
        return
    end
    ids = layerIds{L};
    fprintf('\n--- LAYER %d (%s): groups %s ---\n', L, regionName{L}, ...
        mat2str(ids));
    % v2 layering correction: a layer carries ONLY the envelopes of
    % STRICTLY EARLIER layers (run_1 pass 2 carried every restored
    % envelope into layer 1, inverting the layering: psoas re-solved
    % AGAINST superficial envelopes and regressed from solved to
    % best-infeasible). Deepest layer runs against bone only.
    envsInner = {};
    for LL = 1:L - 1
        for id2 = layerIds{LL}
            if isKey(envRecord, id2)
                envsInner{end + 1} = envRecord(id2); %#ok<AGROW>
            end
        end
    end
    fprintf('obstacle envelopes carried into this layer: %d (from %d earlier layers)\n', ...
        numel(envsInner), L - 1);
    nL = numel(ids);
    templ = fr_pairTemplate();
    pairOut = repmat(templ, 1, nL);
    parfor ii = 1:nL
        id = ids(ii);
        cs = CASES(byId(id));
        x0r = [];
        x0c = [];
        nm = cs.name;
        if wRigid.isKey(nm); x0r = wRigid(nm); end
        if wComp.isKey(nm); x0c = wComp(nm); end
        pairOut(ii) = fr_solvePair(C, cs, kin, resArr, ckDir, logDir, ...
            envsInner, C.solver, x0r, x0c);
    end
    for ii = 1:nL
        id = ids(ii);
        po = pairOut(ii);
        if strcmp(po.rigid.status, 'solved')
            envRecord(id) = po.rigid.envelope;
            solvedEnv(id) = po.rigid.envelope;
        elseif strcmp(po.compliant.status, 'solved')
            envRecord(id) = po.compliant.envelope;
            solvedEnv(id) = po.compliant.envelope;
        end
        fprintf(['  [%2d] %-9s rigid %-16s f %8.3f violEff %.1e | ' ...
            'compliant %-16s f %8.3f violEff %.1e\n'], id, po.name, ...
            po.rigid.status, po.rigid.fBest, po.rigid.violation, ...
            po.compliant.status, po.compliant.fBest, ...
            po.compliant.violation);
        if isfield(po.rigid, 'topViolator') ...
                && strcmp(po.rigid.status, 'best-infeasible')
            fprintf('       rigid top violator: %s\n', po.rigid.topViolator);
        end
        if isfield(po.compliant, 'topViolator') ...
                && strcmp(po.compliant.status, 'best-infeasible')
            fprintf('       compliant top violator: %s\n', ...
                po.compliant.topViolator);
        end
    end
    layerResults{L} = pairOut;
end

% ---- REFINEMENT ROUND 1 -----------------------------------------------------
if exist(stopFile, 'file') == 2
    fprintf(['STOP detected before refinement - stopping cleanly; ' ...
        'checkpoints intact and resumable.\n']);
    fr_writeStatus(runDir, startedAt, true, layerResults, resArr, ...
        envRecord);
    return
end
fprintf(['\n--- REFINEMENT ROUND 1 (all %d groups, snapshot ' ...
    'envelopes) ---\n'], nG);
envSnapshot = cell(1, nG);
for i = 1:nG
    id = allIds(i);
    if isKey(envRecord, id)
        envSnapshot{i} = envRecord(id);
    end
end
refineOut = fr_refineTemplate();
parfor i = 1:nG
    id = allIds(i);
    cs = CASES(byId(id));
    ra = resArr(byName(cs.name));
    others = envSnapshot;
    others(i) = [];
    refineOut(i) = fr_refineOne(C, cs, kin, ra, ckDir, logDir, others, ...
        C.solver);
end
for i = 1:nG
    id = allIds(i);
    ro = refineOut(i);
    if ro.accepted
        envRecord(id) = ro.envelope;
    end
    fprintf('  [%2d] %-9s refine %s f %.3f -> %.3f\n', id, ro.name, ...
        ro.verdict, ro.fPrev, ro.fNew);
end

fprintf('\n=== run_1 complete | wall %.2f h | %s ===\n', ...
    toc(tAll) / 3600, datetime('now'));
fr_writeStatus(runDir, startedAt, false, layerResults, resArr, ...
    envRecord, CASES, refineOut, allIds);
end
