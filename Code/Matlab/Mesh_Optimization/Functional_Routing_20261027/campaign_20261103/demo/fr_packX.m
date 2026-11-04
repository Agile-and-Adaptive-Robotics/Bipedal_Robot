function x = fr_packX(cs, design)
%FR_PACKX design struct -> decision vector (inverse of fr_unpackX).
nMov = numel(cs.movable);
x = zeros(cs.nX, 1);
for m = 1:nMov
    x(3 * m - 2:3 * m) = design.uv(m, :);
end
x(3 * nMov + 1) = design.L0 * 10;
x(3 * nMov + 2) = design.tendon * 10;
x(end) = design.seriesSplit;
end
