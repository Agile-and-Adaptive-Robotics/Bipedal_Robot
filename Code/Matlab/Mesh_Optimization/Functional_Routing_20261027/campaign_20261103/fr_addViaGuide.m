function candidate = fr_addViaGuide(cs, afterRow, point, beforeCrossing)
% Insert a body-fixed guide between consecutive rows on the same body.
% Existing attachment design variables keep their positions in the vector.
next = afterRow + 1;
assert(afterRow >= 1 && next <= numel(cs.rows), 'Invalid guide segment');
body = cs.rowOwning{afterRow};
if nargin<4; beforeCrossing=false; end
assert(beforeCrossing || strcmp(body, cs.rowOwning{next}), 'Guide segment must share a body');
assert(isequal(size(point), [1 3]) && all(isfinite(point)), 'Invalid guide point');
candidate = cs;
newRow = struct('body', body, 'loc', point);
candidate.rows = insert(cs.rows, afterRow, {newRow});
candidate.rowOwning = insert(cs.rowOwning, afterRow, {body});
candidate.seedLoc = [cs.seedLoc(1:afterRow,:); point; cs.seedLoc(next:end,:)];
for m = 1:numel(cs.movable)
    if cs.movable(m).row > afterRow
        candidate.movable(m).row = cs.movable(m).row + 1;
    end
end
for p = 1:numel(cs.passes)
    mask = cs.rowUse{p};
    used = mask(afterRow) && mask(next);
    candidate.rowUse{p} = insert(mask, afterRow, used);
    candidate.rowExpr{p} = insert(cs.rowExpr{p}, afterRow, cs.rowExpr{p}(afterRow));
    candidate.rowFrame{p} = insert(cs.rowFrame{p}, afterRow, cs.rowFrame{p}(afterRow));
    if used
        assert(beforeCrossing || strcmp(cs.rowFrame{p}{afterRow}, cs.rowFrame{p}{next}), ...
            'Guide endpoints must share an assigned frame');
        cross = cs.crossFor{p};
        candidate.crossFor{p} = cross + (cross > nnz(mask(1:afterRow)));
    end
end
candidate.cross = candidate.crossFor{1};
end

function out = insert(in, k, element)
out = [in(1:k), element, in(k+1:end)];
end
