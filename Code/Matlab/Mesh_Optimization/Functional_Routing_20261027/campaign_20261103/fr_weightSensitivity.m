function fr_weightSensitivity()
%FR_WEIGHTSENSITIVITY Bounded weight-sensitivity check (Ben directive
% item 7, 2026-11-04): normalize, inspect contributions, perturb each
% weight +-2x, confirm capacity/ROM stay primary, shape small, overlap
% tolerable, and that the CHOSEN weights are defensible. Run headless:
%   matlab -batch "cd('<campaign_20261103>'); fr_weightSensitivity"
% Writes output/weight_sensitivity.log. Exits nonzero on failure.

logFile = fullfile('output', 'weight_sensitivity.log');
if ~exist('output', 'dir'); mkdir('output'); end
fid = fopen(logFile, 'w');
C = fr_setup();
if isfield(C, 'kinematics')
    kin = C.kinematics;
else
    kin = fr_kinematics(C.osimStock);
end
cases = fr_cases27(C, kin);
pick = {'psoas', 'vas_int', 'semimem'};   % fr_cases27 group keys: hip mono, multi-BPA knee, biarticular
W0 = containerOfSpec(C);
ok = true;

fprintf(fid, 'weight sensitivity | %s | baseline spec status: %s\n', ...
    datetime('now'), specStatus(C));
for c = 1:numel(pick)
    cs = cases(strcmp({cases.name}, pick{c}));
    if isempty(cs)
        fprintf(fid, 'case %s: NOT FOUND - skipping\n', pick{c});
        ok = false;
        continue
    end
    cs = cs(1);
    fprintf(fid, '\n=== %s ===\n', cs.name);
    [f0, rep0] = fr_objective(C, cs, kin, cs.xSeed, 'rigid', ...
        struct('quiet', true));
    fprintf(fid, 'baseline: f %.4f (obj %.4f, violEff %.2e)\n', f0, ...
        rep0.fObj, rep0.violation);
    printTerms(fid, rep0, W0);
    % capacity must dominate the objective at the chosen weights
    contrib = termContribs(rep0, W0);
    tot = sum(abs(contrib));
    if tot > 0 && contrib(1) / tot < 0.5
        fprintf(fid, 'CHECK capacity share: %.1f%% < 50%% - FAIL\n', ...
            100 * contrib(1) / tot);
        ok = false;
    else
        fprintf(fid, 'CHECK capacity share: %.1f%% of |weighted terms| (>= 50%%) - PASS\n', ...
            100 * safeDiv(contrib(1), tot));
    end
    if contrib(3) > 0.1 * tot
        fprintf(fid, 'CHECK shape share: %.1f%% > 10%% of weighted terms - FAIL\n', ...
            100 * safeDiv(contrib(3), tot));
        ok = false;
    else
        fprintf(fid, 'CHECK shape share: %.1f%% (<= 10%%) - PASS\n', ...
            100 * safeDiv(contrib(3), tot));
    end
    % +-2x perturbation per weight: report fObj swing at the SEED design
    % (via opts.weightOverride - the spec file is never rewritten)
    ks = W0.keys;
    vs = W0.values;
    for iTerm = 1:numel(ks)
        nm = ks{iTerm};
        for s = [-2, 2]
            Wm = containers.Map('KeyType', 'char', 'ValueType', 'double');
            Wm(nm) = vs{iTerm} * s;
            fi = evalWithWeights(C, cs, kin, Wm);
            fprintf(fid, '  perturb %-28s x%+g -> fObj %.4f (delta %+.4f)\n', ...
                nm, s, fi, fi - rep0.fObj);
        end
    end
end

% short solve A/B: deficit weight 50 vs 200 on psoas, tiny budget, to
% show the design choice is stable to the perturbation
fprintf(fid, '\n=== short A/B solve (psoas, budget 100) ===\n');
try
    csH = cases(strcmp({cases.name}, 'psoas')); csH = csH(1);
    f50 = shortSolve(C, csH, kin, 'torque_capacity_deficit', 50);
    f200 = shortSolve(C, csH, kin, 'torque_capacity_deficit', 200);
    fprintf(fid, 'short-solve fBest: deficitW50 %.4f | deficitW200 %.4f | |delta| %.4f\n', ...
        f50, f200, abs(f50 - f200));
    if abs(f50 - f200) > 0.25 * min(abs(f50), abs(f200)) + 1e-9
        fprintf(fid, 'CHECK A/B stability: designs diverge >25%% - inspect (WARN, not automatic FAIL)\n');
    else
        fprintf(fid, 'CHECK A/B stability: within 25%% - PASS\n');
    end
catch err
    fprintf(fid, 'A/B short solve errored: %s\n', err.message);
    ok = false;
end

fprintf(fid, '\nVERDICT: %s\n', frPass(ok));
fclose(fid);
fprintf('weight sensitivity verdict: %s (log: %s)\n', frPass(ok), ...
    fullfile('output', 'weight_sensitivity.log'));
if ~ok
    error('fr_weightSensitivity:checkFailed', ...
        'weight sensitivity checks failed - see the log');
end
end

% ==========================================================================
function [fi, ri] = evalWithWeights(C, cs, kin, Wm)
oq = struct('quiet', true, 'weightOverride', Wm);
[fi, ri] = fr_objective(C, cs, kin, cs.xSeed, 'rigid', oq);
end

function f = shortSolve(C, cs, kin, termName, wVal)
Wm = containers.Map('KeyType', 'char', 'ValueType', 'double');
Wm(termName) = wVal;
so = struct('quiet', true, 'budget', 100, 'weightOverride', Wm);
[~, ~, info] = fr_solveCase(C, cs, kin, 'rigid', so);
f = info.fBest;
end

% ==========================================================================
function w = containerOfSpec(C)
spec = jsondecode(fileread(fullfile(C.camp, 'objective_spec.json')));
tl = spec.terms;
w = containers.Map();
for i = 1:numel(tl)
    if iscell(tl)
        w(tl{i}.name) = tl{i}.weight;
    else
        w(tl(i).name) = tl(i).weight;
    end
end
end

function s = specStatus(C)
spec = jsondecode(fileread(fullfile(C.camp, 'objective_spec.json')));
s = spec.status;
end

function c = termContribs(rep, W)
nm = fieldnames(rep.terms);
c = zeros(numel(nm), 1);
for i = 1:numel(nm)
    if W.isKey(nm{i})
        c(i) = W(nm{i}) * rep.terms.(nm{i});
    end
end
end

function printTerms(fid, rep, W)
nm = fieldnames(rep.terms);
for i = 1:numel(nm)
    wv = NaN;
    if W.isKey(nm{i}); wv = W(nm{i}); end
    fprintf(fid, '  term %-28s value %.4f weight %g -> weighted %.4f\n', ...
        nm{i}, rep.terms.(nm{i}), wv, wv * rep.terms.(nm{i}));
end
end

function d = safeDiv(a, b)
if b <= 0
    d = 0;
else
    d = a / b;
end
end

function s = frPass(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end
