function po = fr_solvePair(C, cs, kin, resArr, ckDir, logDir, envsInner, ...
    solver, x0Rigid, x0Compliant)
%FR_SOLVEPAIR RIGID branch then COMPLIANT branch for one group (v2).
%
% Resume: a checkpoint is skipped ONLY when it is FEASIBLE-of-record
% (ck.status == 'solved', i.e. every hard constraint within its stated
% tolerance) and its case signature (nX, xSeed sum) matches the current
% case. 'best-infeasible' checkpoints are NOT treated as completed work:
% the group re-solves. Each branch gets the FULL solve budget;
% opts.envelopes carries the retained deeper-layer obstacles
% (rigid-hard + BPA-corridor-soft per fr_obstacles v2) for both branches.
% x0Rigid / x0Compliant (optional design structs from the run_0 warm-start
% layer) are threaded as opts.x0Override - previous bests as STARTING
% POINTS, never as results of record (Ben 2026-11-04 item 1).

if nargin < 9; x0Rigid = []; end
if nargin < 10; x0Compliant = []; end
po = fr_pairTemplate();
po.name = cs.name;
po.mapId = cs.groupMapId;
glog = fullfile(logDir, ['group_' cs.name '.log']);
ra = fr_resumeRow(resArr, cs.name);

% ---- rigid ---------------------------------------------------------------
if ra.hasRigid
    S = load(ra.ckRigid);
    ckStatus = 'solved';
    if isfield(S, 'status'); ckStatus = S.status; end
    if S.nX ~= cs.nX || abs(S.xSeedSum - sum(cs.xSeed)) > 1e-9
        fr_glog2(glog, sprintf(['RESUME %s rigid: checkpoint case ' ...
            'signature differs (nX %d vs %d) - solving fresh'], ...
            cs.name, S.nX, cs.nX));
    elseif ~strcmp(ckStatus, 'solved')
        fr_glog2(glog, sprintf(['RESUME %s rigid: checkpoint is %s - ' ...
            're-solving (kept as warm start, not as a result)'], ...
            cs.name, ckStatus));
    else
        po.rigid = fr_outFromCk(S, 'solved');
        fr_glog2(glog, sprintf('RESUME skip (feasible checkpoint): rigid f %.4f', ...
            S.fBest));
    end
end
if ~strcmp(po.rigid.status, 'solved')
    opts = struct('quiet', true, 'envelopes', {envsInner}, ...
        'budget', solver.MaxFunctionEvaluations);
    if ~isempty(x0Rigid); opts.x0Override = x0Rigid; end
    po.rigid = fr_groupSolve(C, cs, kin, 'rigid', opts, ...
        fr_ckPath(cs.name, 'rigid', ckDir), glog, 'pass1');
end

% ---- compliant (Xi0-Xi3) - full re-optimization, matched budget -----------
if ra.hasCompliant
    S = load(ra.ckCompliant);
    ckStatus = 'solved';
    if isfield(S, 'status'); ckStatus = S.status; end
    if S.nX ~= cs.nX || abs(S.xSeedSum - sum(cs.xSeed)) > 1e-9
        fr_glog2(glog, sprintf(['RESUME %s compliant: checkpoint case ' ...
            'signature differs (nX %d vs %d) - solving fresh'], ...
            cs.name, S.nX, cs.nX));
    elseif ~strcmp(ckStatus, 'solved')
        fr_glog2(glog, sprintf(['RESUME %s compliant: checkpoint is %s - ' ...
            're-solving (kept as warm start, not as a result)'], ...
            cs.name, ckStatus));
    else
        po.compliant = fr_outFromCk(S, 'solved');
        fr_glog2(glog, sprintf( ...
            'RESUME skip (feasible checkpoint): compliant f %.4f', S.fBest));
    end
end
if ~strcmp(po.compliant.status, 'solved')
    opts = struct('quiet', true, 'envelopes', {envsInner}, ...
        'budget', solver.MaxFunctionEvaluations);
    if ~isempty(x0Compliant); opts.x0Override = x0Compliant; end
    po.compliant = fr_groupSolve(C, cs, kin, 'compliant', opts, ...
        fr_ckPath(cs.name, 'compliant', ckDir), glog, 'pass1');
end
end
