function [pts, bodies, absRows] = fr_locForPass(~, cs, design, ~, passIdx)
%FR_LOCFORPASS Route locations for one pass: nUsed x 3 in each row's
% ASSIGNED frame (cs.rowFrame{k}{absRow}), with movable attachments placed
% by their surface patches (fr_patch.project) in their OWNING body frame
% and re-expressed by the pass's constant matrix where one exists.
%
% master directive item 3: the mount projection keeps every candidate on
% the bone-surface band (nearest cloud point) plus the explicit mount
% standoff in millimeters, so a Cartesian search step cannot carry a mount
% far off the bone.
mask = cs.rowUse{passIdx};
absRows = find(mask);
nU = numel(absRows);
pts = zeros(nU, 3);
bodies = cell(1, nU);
for u = 1:nU
    rAbs = absRows(u);
    p = cs.seedLoc(rAbs, :);
    mvIdx = find([cs.movable.row] == rAbs, 1);
    if ~isempty(mvIdx)
        mv = cs.movable(mvIdx);
        uv = design.uv(mvIdx, :);
        if isfield(mv, 'patch') && ~isempty(mv.patch)
            p = mv.patch.project(uv(1), uv(2), uv(3));
        end
        % (an empty patch = the finish-time seed call before patches are
        % attached: the seed location itself is the right answer)
    end
    ex = cs.rowExpr{passIdx}{rAbs};
    if ~isempty(ex)
        p = (ex(1:3, 1:3) * p(:) + ex(1:3, 4))';
    end
    pts(u, :) = p;
    bodies{u} = cs.rowFrame{passIdx}{rAbs};
end
end
