function ro = fr_refineOne(C, cs, kin, ra, ckDir, logDir, envOthers, ...
    solver)
%FR_REFINEONE Coordinated refinement round 1 for one group (worker-side).
%
% Rigid re-solve (budget C.solver.refineMaxFunctionEvaluations) against
% ALL OTHER groups' ROUND-START envelope snapshot (documented semantics:
% the sequential fr_runLayered round updates envelopes mid-round; the
% parallel round uses the round-start snapshot - the per-group revert
% rule protects either way). ACCEPTED only if the new objective does not
% degrade the retained objective beyond C.refineTolerance; otherwise the
% pass-1 design is KEPT (revert rule) and the rejection is recorded.
% Checkpoint <name>_refine1.mat written either way (resume skips a
% completed refinement). fPrev = the rigid pass-1 checkpoint's fBest.

ro = struct('name', cs.name, 'mapId', cs.groupMapId, ...
    'verdict', 'pending', 'accepted', false, 'fPrev', NaN, 'fNew', NaN, ...
    'envelope', [], 'design', [], 'info', []);
glog = fullfile(logDir, ['group_' cs.name '.log']);

if ra.hasRefine
    S = load(ra.ckRefine);
    ro.verdict = S.verdict;
    ro.accepted = S.accepted;
    ro.fPrev = S.fPrev;
    ro.fNew = S.fNew;
    if S.accepted
        ro.envelope = S.envelope;
        ro.design = S.design;
        ro.info = S.info;
    end
    fr_rlog(glog, sprintf('RESUME skip refine1 (%s): f %.4f', ...
        S.verdict, S.fNew));
    return
end
if ~ra.hasRigid
    ro.verdict = 'skipped-no-rigid';
    fr_rlog(glog, 'REFINE skipped: no rigid checkpoint to refine');
    return
end

S = load(ra.ckRigid);
fPrev = S.fBest;
t0 = tic;
opts = struct('quiet', true, 'envelopes', {envOthers}, ...
    'budget', solver.refineMaxFunctionEvaluations);
try
    [bestNew, repNew, infoNew] = fr_solveCase(C, cs, kin, 'rigid', opts);
    ro.fPrev = fPrev;
    ro.fNew = infoNew.fBest;
    if infoNew.fBest <= fPrev + C.refineTolerance
        ro.verdict = 'accepted';
        ro.accepted = true;
        ro.design = bestNew;
        ro.info = infoNew;
        ro.envelope = fr_envelope(C, cs, bestNew, kin, ...
            cs.passes{cs.envelopePassIdx}.qStructs, 'rigid');
    else
        ro.verdict = 'rejected';
        ro.accepted = false;
        ro.design = S.design;
        ro.info = S.info;
        ro.envelope = S.envelope;
    end
    fr_rlog(glog, sprintf( ...
        ['REFINE1 %s: f %.4f -> %.4f (%s, others %d, %.1f s)'], ...
        cs.name, fPrev, ro.fNew, ro.verdict, numel(envOthers), ...
        toc(t0)));
    if exist(ra.ckRefine, 'file')
        delete(ra.ckRefine);
    end
    ck = struct();
    ck.csName = cs.name;
    ck.mapId = cs.groupMapId;
    ck.branch = 'rigid';
    ck.phase = 'refine1';
    ck.verdict = ro.verdict;
    ck.accepted = ro.accepted;
    ck.fPrev = ro.fPrev;
    ck.fNew = ro.fNew;
    ck.envelope = ro.envelope;
    ck.design = ro.design;
    ck.info = ro.info;
    ck.writtenAt = datetime('now');
    ck.machine = getenv('COMPUTERNAME');
    save(ra.ckRefine, '-struct', 'ck');
catch err
    ro.verdict = 'failed';
    fr_rlog(glog, sprintf('REFINE1 FAILED %s: %s | %s', cs.name, ...
        err.identifier, err.message));
end
end

function fr_rlog(logFile, line)
fid = fopen(logFile, 'a');
if fid == -1
    return
end
fprintf(fid, '%s\n', line);
fclose(fid);
end
