function w = taperKernel_20261022(x, lo, hi, band, wMin)
%TAPERKERNEL_20261022 Shared soft-taper kernel for the 2026-10-22 guards.
%   w = TAPERKERNEL_20261022(X, LO, HI, BAND, WMIN) returns 1 inside
%   [LO, HI] (inclusive), a raised-cosine taper down to WMIN across BAND
%   beyond either edge, and WMIN farther out. Elementwise over X; LO/HI
%   may be scalars or arrays the same size as X (the Item 5 hip-coupling
%   mask passes per-pose limit pairs). Pure function, deterministic, and
%   dependent on no file outside this campaign folder
%   (tmp/paper_bridge/opt_improve_20261022).
%
%   Shared-kernel contract: hipAdductionEnvelope_20261022 (Item 5 mask)
%   and softCornerWeights_20261022 (Item 6 corner weights) MUST both use
%   THIS kernel so the two guards compose without mixed tapers.
%
%   Review-blocker guardrail (plan review 2026-10-22, Item 2): weights
%   from this kernel only ever DOWN-WEIGHT requirement/shortfall terms.
%   They must never rescue slack cells: at slack the insertion force is
%   zeroed by construction (MonoPam_pulley.m:595; BiPam_pulley.m:671), so
%   Torque_ins/tauTarget - 1 == -1 identically there; consumers must
%   EXCLUDE slack cells from any torque test by mask and route them to
%   the RoM/slack side, not re-weight them.

if nargin ~= 5
    error('taperKernel_20261022:nargin', 'Need exactly 5 arguments.');
end
if ~(isscalar(band) && isscalar(wMin))
    error('taperKernel_20261022:args', 'BAND and WMIN must be scalars.');
end
if band <= 0
    error('taperKernel_20261022:args', 'Need BAND > 0.');
end
if wMin < 0 || wMin > 1
    error('taperKernel_20261022:args', 'Need 0 <= WMIN <= 1.');
end

x = double(x);
lo = double(lo);
hi = double(hi);
sz = size(x);
if ~isequal(size(lo), sz) || ~isequal(size(hi), sz)
    if isscalar(lo) && isscalar(hi)
        lo = repmat(lo, sz);
        hi = repmat(hi, sz);
    else
        error('taperKernel_20261022:args', ...
            'LO/HI must be scalars or the same size as X.');
    end
end
if any(hi(:) < lo(:))
    error('taperKernel_20261022:args', 'Every HI must be >= its LO.');
end

w = ones(sz);
inside = (x >= lo) & (x <= hi);
w(inside) = 1;

% Normalized distance beyond the nearer edge, clipped to [0, 1].
tHigh = (x - hi) / band;   % > 0 beyond the upper edge
tLow = (lo - x) / band;    % > 0 below the lower edge
t = min(max(max(tHigh, tLow), 0), 1);

out = ~inside;
w(out) = wMin + (1 - wMin) * (1 + cos(pi * t(out))) / 2;
% At t == 1 this evaluates to wMin exactly (cos(pi) == -1), so the floor
% is bit-exact for poses beyond the taper band.
end
