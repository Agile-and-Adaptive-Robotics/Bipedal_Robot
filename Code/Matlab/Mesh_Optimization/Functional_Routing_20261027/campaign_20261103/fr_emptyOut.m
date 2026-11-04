function out = fr_emptyOut()
%FR_EMPTYOUT Fixed field set for fr_groupSolve outputs (parfor struct
% array slicing requires uniform fields in a fixed order).
out = struct();
out.name = '';
out.mapId = NaN;
out.branch = '';
out.phase = '';
out.status = '';
out.fBest = NaN;
out.violation = NaN;
out.fObj = NaN;
out.envelope = [];
out.design = [];
out.rep = [];
out.info = [];
out.nX = NaN;
out.xSeedSum = NaN;
out.elapsedSec = NaN;
end
