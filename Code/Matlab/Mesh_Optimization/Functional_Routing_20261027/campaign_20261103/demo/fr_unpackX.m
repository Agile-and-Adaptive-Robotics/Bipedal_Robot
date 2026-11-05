function design = fr_unpackX(cs, x)
%FR_UNPACKX Decision vector -> design struct.
%
% x layout (per case, see fr_finishCase): for each movable attachment
% [u; v; standoff_mm], then resting length L0 (dm), tendon length (dm),
% series split fraction. Scaling keeps patternsearch conditioning sane.
design = struct();
nMov = numel(cs.movable);
design.uv = zeros(nMov, 3);
for m = 1:nMov
    design.uv(m, :) = x(3 * m - 2:3 * m);
end
design.L0 = x(3 * nMov + 1) / 10;
design.tendon = x(3 * nMov + 2) / 10;
design.seriesSplit = min(1, max(1e-3, x(end)));
end

function x = fr_packX(cs, design) %#ok<DEFNU>
%FR_PACKX design struct -> decision vector (inverse of fr_unpackX)
nMov = numel(cs.movable);
x = zeros(cs.nX, 1);
for m = 1:nMov
    x(3 * m - 2:3 * m) = design.uv(m, :);
end
x(3 * nMov + 1) = design.L0 * 10;
x(3 * nMov + 2) = design.tendon * 10;
x(end) = design.seriesSplit;
end
