function y = fr_splinemm(x, yv, xq)
%FR_SPLINEMM Natural cubic spline (OpenSim SimmSpline convention).
[xs, iu] = sort(x(:));
ys = yv(:); ys = ys(iu);
n = numel(xs);
h = diff(xs);
A = zeros(n);
b = zeros(n, 1);
A(1, 1) = 1; A(n, n) = 1;
for i = 2:n-1
    A(i, i-1) = h(i-1);
    A(i, i) = 2 * (h(i-1) + h(i));
    A(i, i+1) = h(i);
    b(i) = 6 * ((ys(i+1) - ys(i)) / h(i) - (ys(i) - ys(i-1)) / h(i-1));
end
M = A \ b;
xq = xq(:);
xq = min(max(xq, xs(1)), xs(end));   % clamp to the spline domain
idx = discretize(xq, xs, 'IncludedEdge', 'right');
idx = min(max(idx, 1), n-1);   % segment containing xq
y = zeros(size(xq));
for k = 1:numel(xq)
    i = idx(k);
    hi = h(i);
    A1 = (xs(i+1) - xq(k)) / hi;
    A2 = (xq(k) - xs(i)) / hi;
    y(k) = A1 * ys(i) + A2 * ys(i+1) + ...
        ((A1^3 - A1) * M(i) + (A2^3 - A2) * M(i+1)) * hi^2 / 6;
end
y = reshape(y, 1, []);
end
