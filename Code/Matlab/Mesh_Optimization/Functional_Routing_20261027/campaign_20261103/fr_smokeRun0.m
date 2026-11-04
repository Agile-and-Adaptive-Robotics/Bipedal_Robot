function ok = fr_smokeRun0()
%FR_SMOKERUN0 Launch-plumbing smoke for run_0 (foreground, one group).
%
% Exercises exactly the plumbing the detached driver uses, on ONE real
% group (psoas, full budgets so its checkpoints are legitimately reused
% by the campaign via resume): manifest write -> fr_solvePair (rigid then
% compliant) -> resume scan sees the checkpoints -> envelope-of-record
% restore from the checkpoint. NO parfor/pool here; the pool path is
% verified by the detached launch's first-layer evidence.

C = fr_setup();
runDir = fullfile(C.camp, 'run_0');
ckDir = fullfile(runDir, 'checkpoints');
logDir = fullfile(runDir, 'logs');
if ~exist(ckDir, 'dir'); mkdir(ckDir); end
if ~exist(logDir, 'dir'); mkdir(logDir); end

kin = fr_kinematics(C.osimStock);
CASES = fr_cases27(C, kin);
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
[layerIds, regionName] = fr_layerIds(tc, CASES);
manFile = fullfile(runDir, 'manifest.json');
if ~exist(manFile, 'file')
    fr_writeManifest(manFile, C, tc, CASES, layerIds, regionName);
    fprintf('smoke: manifest written (%s)\n', manFile);
else
    fprintf('smoke: manifest already present\n');
end

cs = CASES(strcmp({CASES.name}, 'psoas'));
cs = cs(1);
resArr = fr_scanCheckpoints(ckDir, CASES);
t0 = tic;
po = fr_solvePair(C, cs, kin, resArr, ckDir, logDir, {}, C.solver);
fprintf(['smoke solvePair wall %.1f s | rigid %s f %.3f viol %.1e | ' ...
    'compliant %s f %.3f viol %.1e\n'], toc(t0), po.rigid.status, ...
    po.rigid.fBest, po.rigid.violation, po.compliant.status, ...
    po.compliant.fBest, po.compliant.violation);

resArr2 = fr_scanCheckpoints(ckDir, CASES);
idx = find(strcmp({resArr2.name}, 'psoas'), 1);
ok = resArr2(idx).hasRigid && resArr2(idx).hasCompliant && ...
    strcmp(po.rigid.status, 'solved') && ...
    strcmp(po.compliant.status, 'solved');
fprintf('smoke resume scan: hasRigid %d hasCompliant %d\n', ...
    resArr2(idx).hasRigid, resArr2(idx).hasCompliant);

S = load(resArr2(idx).ckRigid);
fprintf(['smoke envelope restore: capsules %d poses x %d, spheres %d ' ...
    'poses, fBest %.3f\n'], numel(S.envelope.capsules), ...
    numel(S.envelope.capsules{1}), numel(S.envelope.spheres), S.fBest);
ok = ok && ~isempty(S.envelope.capsules);
if ok
    fprintf('SMOKE PASSED\n');
else
    fprintf('SMOKE FAILED\n');
end
end
