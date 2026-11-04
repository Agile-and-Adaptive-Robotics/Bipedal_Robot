function fr_launch_run0()
%FR_LAUNCH_RUN0 run_0_weights_placeholder campaign driver (ONE process).
%
% Launched detached exactly as:
%   "D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd( ...
%    '.../Functional_Routing_20261027'); addpath('campaign_20261103'); ...
%    fr_launch_run0"
%
% Structure (the launcher ask, verbatim):
%   - parpool capped at 10 (AGENTS.md easteregg2; env vars set BEFORE
%     parpool so workers do not oversubscribe BLAS threads);
%   - layer loop SEQUENTIAL in the contract's layerOrder (inner layers'
%     occupied envelopes are HARD obstacles for outer layers, fr_obstacles
%     constraint 7 via opts.envelopes);
%   - parfor over groups WITHIN each layer (within-layer groups are
%     placed in parallel - no mutual obstacles in pass 1; the refinement
%     round coordinates every envelope);
%   - per group the RIGID branch then the COMPLIANT branch (Xi0-Xi3,
%     fr_mechBuild), both full re-optimizations with matched budgets
%     (fr_solvePair);
%   - refinement round 1 after the first deep-outward pass: rigid re-solve
%     (budget 150) against ALL OTHER groups' round-start envelope
%     snapshot, accepted only if the retained objective is not degraded
%     beyond C.refineTolerance (revert rule, fr_refineOne);
%   - checkpoints: per group per branch .mat under run_0/checkpoints/
%     written as each completes; resume skips completed (group, branch)
%     pairs and a completed refinement (fr_scanCheckpoints);
%   - run_0/manifest.json at startup (weights snapshot + checksum, seeds,
%     budgets, group list, layer order, machine); run_0/logs/master.log +
%     per-group logs; run_0/STATUS.json at the end;
%   - a run_0/STOP file stops cleanly at the next layer/round boundary
%     (checkpoints stay intact and resumable). Delete STOP before
%     relaunching.

C = fr_setup();
runDir = fullfile(C.camp, 'run_0');
ckDir = fullfile(runDir, 'checkpoints');
logDir = fullfile(runDir, 'logs');
if ~exist(ckDir, 'dir'); mkdir(ckDir); end
if ~exist(logDir, 'dir'); mkdir(logDir); end
masterLog = fullfile(logDir, 'master.log');
diary(masterLog);
cleanupObj = onCleanup(@() diary('off'));

tAll = tic;
startedAt = datetime('now');
fprintf('=== fr_launch_run0 | %s | machine %s ===\n', startedAt, ...
    getenv('COMPUTERNAME'));
fprintf('run dir: %s\n', runDir);
stopFile = fullfile(runDir, 'STOP');
if exist(stopFile, 'file') == 2
    fprintf(['STOP file present at launch - exiting without work ' ...
        '(delete run_0/STOP to relaunch).\n']);
    fr_writeStatus(runDir, startedAt, true, [], []);
    return
end

rng(C.rngSeed);
kin = fr_kinematics(C.osimStock);
CASES = fr_cases27(C, kin);
nG = numel(CASES);
byId = containers.Map('KeyType', 'double', 'ValueType', 'double');
for c = 1:nG
    byId(CASES(c).groupMapId) = c;
end
byName = containers.Map('KeyType', 'char', 'ValueType', 'double');
for c = 1:nG
    byName(CASES(c).name) = c;
end

% ---- contract layer order (regions in order, ids within region) ----------
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
[layerIds, regionName] = fr_layerIds(tc, CASES);
allIds = [];
for r = 1:numel(layerIds)
    allIds = [allIds, layerIds{r}]; %#ok<AGROW>
end
fprintf('layer order: %s\n', mat2str(allIds));

% ---- resume scan ----------------------------------------------------------
resArr = fr_scanCheckpoints(ckDir, CASES);
fprintf(['resume scan: %d/%d rigid, %d/%d compliant, %d/%d refine ' ...
    'checkpoints found\n'], sum([resArr.hasRigid]), nG, ...
    sum([resArr.hasCompliant]), nG, sum([resArr.hasRefine]), nG);

% ---- manifest (written once; a resumed run keeps the original) -----------
manFile = fullfile(runDir, 'manifest.json');
if ~exist(manFile, 'file')
    fr_writeManifest(manFile, C, tc, CASES, layerIds, regionName);
    fprintf('manifest written: %s\n', manFile);
else
    fprintf('manifest kept (resumed run): %s\n', manFile);
end

% ---- envelope-of-record restore (rigid; refine1 accepted overrides) ------
envRecord = containers.Map('KeyType', 'double', 'ValueType', 'any');
for i = 1:numel(resArr)
    ra = resArr(i);
    if ra.hasRefine
        S = load(ra.ckRefine);
        if S.accepted
            envRecord(S.mapId) = S.envelope;
            continue
        end
    end
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

% ---- pool: env vars BEFORE parpool, capped at 10 --------------------------
setenv('OMP_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS', '1');
if isempty(gcp('nocreate'))
    pool = parpool('Processes', C.parallelMaxWorkers);
    fprintf('parpool: %d workers\n', pool.NumWorkers);
else
    fprintf('parpool already running: %d workers\n', ...
        gcp('nocreate').NumWorkers);
end

% ---- PASS 1: sequential layers, parfor within layer ------------------------
layerResults = cell(1, numel(layerIds));
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
    envsInner = envRecord.values();
    fprintf('obstacle envelopes carried into this layer: %d\n', ...
        numel(envsInner));
    nL = numel(ids);
    templ = fr_pairTemplate();
    pairOut = repmat(templ, 1, nL);
    parfor ii = 1:nL
        id = ids(ii);
        cs = CASES(byId(id));
        pairOut(ii) = fr_solvePair(C, cs, kin, resArr, ckDir, logDir, ...
            envsInner, C.solver);
    end
    for ii = 1:nL
        id = ids(ii);
        po = pairOut(ii);
        if strcmp(po.rigid.status, 'solved')
            envRecord(id) = po.rigid.envelope;
        elseif strcmp(po.compliant.status, 'solved')
            envRecord(id) = po.compliant.envelope;
        end
        fprintf(['  [%2d] %-9s rigid %-7s f %8.3f viol %.1e | ' ...
            'compliant %-7s f %8.3f viol %.1e\n'], id, po.name, ...
            po.rigid.status, po.rigid.fBest, po.rigid.violation, ...
            po.compliant.status, po.compliant.fBest, ...
            po.compliant.violation);
    end
    layerResults{L} = pairOut;
end

% ---- REFINEMENT ROUND 1 (round-start envelope snapshot) --------------------
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

% ---- final status ----------------------------------------------------------
fprintf('\n=== run_0 complete | wall %.2f h | %s ===\n', ...
    toc(tAll) / 3600, datetime('now'));
fr_writeStatus(runDir, startedAt, false, layerResults, resArr, ...
    envRecord, CASES, refineOut, allIds);
end
