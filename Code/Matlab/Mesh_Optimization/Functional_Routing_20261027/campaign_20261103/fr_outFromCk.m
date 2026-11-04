function out = fr_outFromCk(S, status)
%FR_OUTFROMCK Restore a fr_groupSolve-shaped output from a checkpoint.
out = fr_emptyOut();
out.name = S.csName;
out.mapId = S.mapId;
out.branch = S.branch;
out.phase = S.phase;
out.status = status;
out.fBest = S.fBest;
out.violation = S.violation;
out.fObj = S.rep.fObj;
out.envelope = S.envelope;
out.design = S.design;
out.rep = S.rep;
out.info = S.info;
out.nX = S.nX;
out.xSeedSum = S.xSeedSum;
out.elapsedSec = S.elapsedSec;
end
