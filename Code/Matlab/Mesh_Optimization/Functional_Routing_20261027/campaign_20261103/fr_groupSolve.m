function out = fr_groupSolve(C, cs, kin, branch, opts, ckFile, logFile, ...
    phaseLabel)
%FR_GROUPSOLVE Solve one (group, branch) and write its checkpoint + log.
%
% Worker-side helper for fr_launch_run0's parfor. Solves with
% fr_solveCase (patternsearch, opts.budget, opts.envelopes = retained
% deeper-layer obstacles), builds the pose-dependent occupied envelope of
% the BEST design (fr_envelope, same branch), writes the checkpoint mat
% (explicit variable names - never a bare save after a struct load) and
% appends structured lines to the group's log file (ONE writer per file:
% safe under parfor because each group owns exactly one file).
%
% On failure: no checkpoint is written (so a relaunch RETRIES the group)
% and the failure is logged. out.status is one of
%   'solved' | 'failed'
% out fields are FIXED so parfor can slice the struct array.

t0 = tic;
out = fr_emptyOut();
out.name = cs.name;
out.mapId = cs.groupMapId;
out.branch = branch;
out.phase = phaseLabel;

fr_glog(logFile, sprintf('=== %s %s | %s | budget %d | envelopes %d ===', ...
    phaseLabel, branch, datetime('now'), opts.budget, ...
    numel(opts.envelopes)));
try
    [best, rep, info] = fr_solveCase(C, cs, kin, branch, opts);
    env = fr_envelope(C, cs, best, kin, ...
        cs.passes{cs.envelopePassIdx}.qStructs, branch);
    out.status = 'solved';
    out.fBest = info.fBest;
    out.violation = rep.violation;
    out.fObj = rep.fObj;
    out.envelope = env;
    out.design = best;
    out.rep = rep;
    out.info = info;
    out.nX = cs.nX;
    out.xSeedSum = sum(cs.xSeed);
    out.elapsedSec = toc(t0);
    if exist(ckFile, 'file')
        delete(ckFile);
    end
    ck = struct();
    ck.csName = cs.name;
    ck.mapId = cs.groupMapId;
    ck.branch = branch;
    ck.phase = phaseLabel;
    ck.design = best;
    ck.rep = rep;
    ck.info = info;
    ck.envelope = env;
    ck.fBest = info.fBest;
    ck.violation = rep.violation;
    ck.nX = cs.nX;
    ck.xSeedSum = sum(cs.xSeed);
    ck.elapsedSec = out.elapsedSec;
    ck.writtenAt = datetime('now');
    ck.machine = getenv('COMPUTERNAME');
    save(ckFile, '-struct', 'ck');
    fr_glog(logFile, sprintf( ...
        ['SOLVED %s %s: fSeed %.4f -> fBest %.4f (obj %.4f, viol ' ...
        '%.2e) | capsules(pose1) %d | %.1f s'], cs.name, branch, ...
        info.fSeed, info.fBest, rep.fObj, rep.violation, ...
        numel(env.capsules{1}), out.elapsedSec));
catch err
    out.status = 'failed';
    out.elapsedSec = toc(t0);
    fr_glog(logFile, sprintf('FAILED %s %s: %s | %s', cs.name, branch, ...
        err.identifier, err.message));
end
end

% ==========================================================================
function fr_glog(logFile, line)
fid = fopen(logFile, 'a');
if fid == -1
    return
end
fprintf(fid, '%s\n', line);
fclose(fid);
end
