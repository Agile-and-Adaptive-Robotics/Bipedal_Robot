function [best, rep, info] = fr_solveCase(C, cs, kin, xiBranch, opts)
%FR_SOLVECASE Optimize one demonstrator case with patternsearch.
%
% Decision variables (per case, fr_unpackX): each movable attachment's
% surface-patch coordinates (u, v in [-1,1]^2) and mount standoff
% (millimeters), then resting length, tendon length, and the series split
% fraction (initialized at equal halves, master item 2). Bounds have
% NONZERO extent by construction (checked by the gate).
%
% The solver is patternsearch (deterministic given the seed; available on
% this machine per AGENTS.md) with an internal penalty objective (see
% fr_objective); retained deeper-layer envelopes arrive inside opts.
% opts.budget overrides MaxFunctionEvaluations (gate uses a short budget).
%
% best = best design; rep = full evaluation report at best; info = solver
% diagnostics (f including violations at seed/best, counts).
if nargin < 5
    opts = struct();
end
if ~isfield(opts, 'budget') || isempty(opts.budget)
    opts.budget = C.solver.MaxFunctionEvaluations;
end
if ~isfield(opts, 'quiet') || isempty(opts.quiet)
    opts.quiet = true;
end
oq = opts;
oq.quiet = true;

% seed evaluation (x0Override = warm start, e.g. from an unconstrained
% solve that the constrained solve must repair)
x0 = cs.xSeed;
if isfield(opts, 'x0Override') && ~isempty(opts.x0Override)
    x0 = opts.x0Override;
end
if isfield(opts, 'penaltyScale') && ~isempty(opts.penaltyScale)
    oq.penaltyScale = opts.penaltyScale;
end
[f0, rep0] = fr_objective(C, cs, kin, x0, xiBranch, oq);

obj = @(x) fr_objective(C, cs, kin, x, xiBranch, oq);
optsPs = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', opts.budget, ...
    'MeshTolerance', C.solver.MeshTolerance, ...
    'StepTolerance', C.solver.StepTolerance, ...
    'Display', 'none', ...
    'UseCompleteSearch', false);
% deterministic poll/search order
optsPs = optimoptions(optsPs, 'PollMethod', 'GPSPositiveBasis2N');
[xb, fb] = patternsearch(obj, x0, [], [], [], [], cs.lb, cs.ub, [], ...
    optsPs);
[~, repB] = fr_objective(C, cs, kin, xb, xiBranch, oq);
% keep the SEED if the search failed to improve f (penalty-inclusive)
if fb < f0
    best = fr_unpackX(cs, xb);
    rep = repB;
    info.fBest = fb;
    info.xBest = xb;
else
    best = fr_unpackX(cs, x0);
    rep = rep0;
    info.fBest = f0;
    info.xBest = x0;
end
info.fSeed = f0;
info.xSeed = x0;
info.budget = opts.budget;
info.branch = xiBranch;
info.improved = info.fBest < info.fSeed;
if ~opts.quiet
    fprintf(['SOLVE %-16s [%s] budget %d: fSeed %.4f -> fBest %.4f ' ...
        '(improved %d)\n'], cs.name, xiBranch, opts.budget, ...
        info.fSeed, info.fBest, info.improved);
end
end
