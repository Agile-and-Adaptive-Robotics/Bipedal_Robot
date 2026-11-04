function po = fr_solvePair(C, cs, kin, resArr, ckDir, logDir, envsInner, ...
    solver)
%FR_SOLVEPAIR RIGID branch then COMPLIANT branch for one group.
%
% Resume: completed (group, branch) checkpoints are skipped (restored via
% fr_outFromCk) unless the checkpoint's case signature (nX, xSeed sum)
% disagrees with the current case, in which case the group solves fresh.
% Each branch gets the FULL solve budget; opts.envelopes carries the
% retained deeper-layer obstacles for both branches.

po = fr_pairTemplate();
po.name = cs.name;
po.mapId = cs.groupMapId;
glog = fullfile(logDir, ['group_' cs.name '.log']);
ra = fr_resumeRow(resArr, cs.name);

% ---- rigid ---------------------------------------------------------------
if ra.hasRigid
    S = load(ra.ckRigid);
    if S.nX ~= cs.nX || abs(S.xSeedSum - sum(cs.xSeed)) > 1e-9
        fr_glog2(glog, sprintf(['RESUME %s rigid: checkpoint case ' ...
            'signature differs (nX %d vs %d) - solving fresh'], ...
            cs.name, S.nX, cs.nX));
    else
        po.rigid = fr_outFromCk(S, 'solved');
        fr_glog2(glog, sprintf('RESUME skip (checkpoint): rigid f %.4f', ...
            S.fBest));
    end
end
if ~strcmp(po.rigid.status, 'solved')
    opts = struct('quiet', true, 'envelopes', {envsInner}, ...
        'budget', solver.MaxFunctionEvaluations);
    po.rigid = fr_groupSolve(C, cs, kin, 'rigid', opts, ...
        fr_ckPath(cs.name, 'rigid', ckDir), glog, 'pass1');
end

% ---- compliant (Xi0-Xi3) - full re-optimization, matched budget -----------
if ra.hasCompliant
    S = load(ra.ckCompliant);
    if S.nX ~= cs.nX || abs(S.xSeedSum - sum(cs.xSeed)) > 1e-9
        fr_glog2(glog, sprintf(['RESUME %s compliant: checkpoint case ' ...
            'signature differs (nX %d vs %d) - solving fresh'], ...
            cs.name, S.nX, cs.nX));
    else
        po.compliant = fr_outFromCk(S, 'solved');
        fr_glog2(glog, sprintf( ...
            'RESUME skip (checkpoint): compliant f %.4f', S.fBest));
    end
end
if ~strcmp(po.compliant.status, 'solved')
    opts = struct('quiet', true, 'envelopes', {envsInner}, ...
        'budget', solver.MaxFunctionEvaluations);
    po.compliant = fr_groupSolve(C, cs, kin, 'compliant', opts, ...
        fr_ckPath(cs.name, 'compliant', ckDir), glog, 'pass1');
end
end
