function W = softCornerWeights_20261022(theta1Deg, theta2Deg, env1, env2, opts)
%SOFTCORNERWEIGHTS_20261022 Soft extreme-pose penalty weights (plan Item 6).
%   W = SOFTCORNERWEIGHTS_20261022(THETA1DEG, THETA2DEG, ENV1, ENV2, OPTS)
%   returns per-pose weights for the requirement envelope on a two-DOF
%   pose grid, per the pulley sweep directive §5: "poses like full flexion
%   combined with abduction may be achievable by humans but the solver
%   must NOT be over-penalized for failing to reach them — soft-penalize
%   or cap the requirement envelope at extreme combined poses rather than
%   enforcing full coverage there." This implements the DOWN-WEIGHT branch
%   of that order.
%
%   Weight law: W = 1 inside BOTH per-joint task envelopes (bit-exact 1,
%   so a mid-envelope failure can never be softened — anti-gaming), then a
%   raised-cosine taper per axis across OPTS.band beyond either edge,
%   with the PRODUCT clamped at the global floor OPTS.wMin (poses beyond
%   both envelopes rest AT wMin, not wMin^2 — the plan's "taper to a
%   floor wMin", and the softest allowed down-weight; humans may reach
%   the corners, so the requirement is down-weighted, not removed).
%   Kernel = taperKernel_20261022, shared with
%   hipAdductionEnvelope_20261022 (Item 5) so the guards compose.
%   OPTS.hipMask (same size as the grid, e.g. the Item 5 coupling mask)
%   is AND-ed in via min(W, hipMask) — that is how the combined
%   flexion-x-adduction corner is caught even when each axis alone sits
%   inside its envelope.
%
%   Consumption rule (driver copies only, never Ben's objectives):
%   weighted margin = min over cells of W .* margin; report BOTH the raw
%   and the softened worst margin so no softened number can be mistaken
%   for a raw one. W == 1 exactly mid-envelope, so the softening can
%   never hide a functional-RoM failure.
%
%   Review-blocker guardrail (plan review 2026-10-22, Item 2): these
%   weights multiply the REQUIREMENT side only. They must never be used
%   to rescue slack cells: at slack, insertion force is zeroed by
%   construction (MonoPam_pulley.m:595; BiPam_pulley.m:671), so
%   |Torque_ins|/tauTarget - 1 == -1 identically there. Consumers must
%   define torque tests over TAUT, FEASIBLE cells only (mask slack cells
%   out and route them to the RoM/slack side), then apply W.
%
%   PENDING-BEN defaults: band = 10 deg, wMin = 0.1 (plan Item 6); both
%   are printed when OPTS.verbose and must be revisited by Ben.

if nargin < 5
    opts = struct();
end
band = getDef(opts, 'band', 10);    % PENDING-BEN
wMin = getDef(opts, 'wMin', 0.1);   % PENDING-BEN
verbose = getDef(opts, 'verbose', false);
hipMask = getDef(opts, 'hipMask', []);

t1 = double(theta1Deg);
t2 = double(theta2Deg);
if ~isequal(size(t1), size(t2))
    error('softCornerWeights_20261022:args', ...
        'theta1Deg and theta2Deg must be the same size.');
end
env1 = double(env1(:).');
env2 = double(env2(:).');
if numel(env1) ~= 2 || numel(env2) ~= 2
    error('softCornerWeights_20261022:args', ...
        'env1/env2 must each be [lo hi].');
end

W = taperKernel_20261022(t1, env1(1), env1(2), band, wMin) .* ...
    taperKernel_20261022(t2, env2(1), env2(2), band, wMin);
W = max(wMin, W);   % global floor: double-extreme corners rest AT wMin

if ~isempty(hipMask)
    if ~isequal(size(hipMask), size(W))
        error('softCornerWeights_20261022:args', ...
            'opts.hipMask must be the same size as the pose grid.');
    end
    W = min(W, hipMask);   % AND semantics: the tighter guard wins
end

if verbose
    fprintf('softCornerWeights_20261022: band=%g deg (PENDING-BEN), ' ...
        + 'wMin=%g (PENDING-BEN), hipMask=%d\n', band, wMin, ...
        ~isempty(hipMask));
end
end

function v = getDef(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
