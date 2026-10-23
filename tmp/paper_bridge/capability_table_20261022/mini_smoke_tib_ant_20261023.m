% Mini Smoke Tib Ant 20261023
% GRADING-ONLY mini smoke (NO optimizer, NO ladder -- the ask prohibits
% solves): the established smoke site tib_ant_r (the 10-22 survey's own
% smoke choice), graded on its EXACT human-CSV 25-point primary grid with
% the subtalar pinned at its default, using the survey's exact meet/active
% definitions (Opt_run_BiPulleySpecs27_20261022 grade()), on the extended
% masters builder biPulleySpecsFromMasters_20261023. Sizing: the best of
% the same 5-seed box by taut count (deterministic, no search).
% Purpose: give the campaign-shaped meet/active numbers for the result
% schema with zero optimization anywhere.

function out = mini_smoke_tib_ant_20261023()

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(scriptDir);

pres = 620; xi0 = 0.005; xi1 = 5e5; xi2 = 1e4; wraps = 6; fitn = 0.025;

sp = biPulleySpecsFromMasters_20261023({'tib_ant_r'});
sp = sp(1);

% Human graded curve (single source of truth), actuator 16 / ankle_angle_r
curves = loadHuman_(fullfile(scriptDir, 'human_torque_curves_20261022.csv'));
HC = [];
es = curves(16);
for e = 1:numel(es)
    if strcmp(es(e).dof, 'ankle_angle_r')
        HC = es(e);
    end
end
assert(~isempty(HC), 'no human curve for actuator 16 ankle_angle_r');
anglesDeg = HC.anglesDeg(:).';
nA = numel(anglesDeg);

% Graded transforms: primary slot 1 swept over the EXACT CSV angles, slot 2
% pinned at 0 and duplicated (N2 = 2), the survey's buildGradedTransforms
% slot-1 layout.
T = zeros(4, 4, nA, 2);
for ii = 1:nA
    T(:, :, ii, 1) = hingeT_(sp, 1, anglesDeg(ii));
end
T(:, :, 1, 2) = hingeT_(sp, 2, 0);
T(:, :, 2, 2) = hingeT_(sp, 2, 0);

% Lchar on the graded grid (survey convention)
Lrig = rigidLens_(sp.points, sp.cross, T, nA, 2);
Lchar = max(Lrig(:));

% 5-seed box (survey box seeds), best by taut count -- deterministic
ubT = max(0.30 * Lchar, 0.06);
seeds = [0.60 * Lchar, 0.05; 0.30 * Lchar, 0.02; 0.90 * Lchar, 0.02; ...
         0.30 * Lchar, ubT; 0.90 * Lchar, ubT];
slot = sp.kin.primarySlot;
axisU = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));
tauHuman = repmat(reshape(HC.tauNm, nA, 1), 1, 2);
best = struct('meetPct', 0, 'activePct', 0, 'seed', [NaN NaN], 'score', -inf);
for sIdx = 1:size(seeds, 1)
    rest = seeds(sIdx, 1); tendon = seeds(sIdx, 2); kmax = 0.75 * rest;
    obj = BiPam_pulley2_20261022(sp.name, sp.points, sp.cross, 20, T, ...
        rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, wraps, sp.pulley);
    T3 = obj.Torque_ins(:, :, :, slot);
    tauB = reshape(sum(T3 .* reshape(axisU, 1, 3, 1), 2), nA, 2);
    strain = reshape(obj.strain_p, nA, 2);
    slack = logical(obj.PulleySlack);
    infeas = logical(obj.PulleyInfeasible);
    taut = ~slack & ~infeas & isfinite(tauB);
    meet = taut & isfinite(tauHuman) & ...
        (sign(tauB) == sign(tauHuman) | abs(tauHuman) <= 1e-9) & ...
        abs(tauB) >= abs(tauHuman) - 1e-12;
    active = taut & isfinite(strain) & strain > 0 & abs(tauB) > 0;
    meetPct = 100 * sum(meet(:)) / numel(meet);
    activePct = 100 * sum(active(:)) / numel(active);
    score = sum(taut(:));
    fprintf('seed [%.4f %.4f]: taut %d, meet %.1f%%, active %.1f%%\n', ...
        rest, tendon, score, meetPct, activePct);
    if score > best.score
        best = struct('meetPct', meetPct, 'activePct', activePct, ...
            'seed', seeds(sIdx, :), 'score', score);
    end
end
fprintf(['MINI SMOKE tib_ant_r (grading only, no optimizer): meet %.1f%%, ' ...
    'active %.1f%% at seed [%.4f %.4f]\n'], best.meetPct, ...
    best.activePct, best.seed(1), best.seed(2));
out = best;
end

function T = hingeT_(sp, slot, thetaDeg)
th = thetaDeg * pi / 180;
if strcmpi(strtrim(sp.kin.types{slot}), 'fixed')
    R = eye(3);
    p = [0, 0, 0];
else
    a = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));
    K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
    R = eye(3) + sin(th) * K + (1 - cos(th)) * (K * K);
    p = sp.kin.pivots(slot, :);
    if sp.kin.kneeSlot == slot && ~isempty(sp.kin.kneeFit.thetaY)
        thK = sp.kin.kneeFit.thetaY;
        xK = interp1(sp.kin.kneeFit.thetaX, sp.kin.kneeFit.x, thK, ...
            'linear', 'extrap');
        p = [interp1(thK, xK, th, 'linear', 'extrap'), ...
            interp1(thK, sp.kin.kneeFit.y, th, 'linear', 'extrap'), p(3)];
    end
end
T = [R, p(:); 0, 0, 0, 1];
end

function L = rigidLens_(points, cross, T, N1, N2)
nPts = size(points, 1);
L = zeros(N1, N2);
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for iii = 1:N2
        T2 = T(:, :, iii, 2);
        tot = 0;
        for r = 1:nPts - 1
            pA = mapRow_(points, cross, r, T1, T2);
            pB = mapRow_(points, cross, r + 1, T1, T2);
            tot = tot + norm(pA - pB);
        end
        L(ii, iii) = tot;
    end
end
end

function p = mapRow_(points, cross, r, T1, T2)
if r < cross(1)
    f = 1;
elseif cross(1) == cross(2)
    f = 3;
elseif r < cross(2)
    f = 2;
else
    f = 3;
end
p = points(r, :);
while f > 1
    if f == 3
        rr = [p, 1] * T2.';  p = rr(1:3);  f = 2;
    else
        rr = [p, 1] * T1.';  p = rr(1:3);  f = 1;
    end
end
end

function curves = loadHuman_(csvPath)
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
cols = raw.Properties.VariableNames;
angCols = orderCols_(cols, 'theta_deg_');
tauCols = orderCols_(cols, 'tau_h_Nm_');
curves = containers.Map('KeyType', 'double', 'ValueType', 'any');
for r = 1:height(raw)
    idn = double(raw.actuator_id(r));
    if isnan(idn); continue; end
    angles = cellfun(@(c) raw.(c)(r), angCols);
    taus = cellfun(@(c) raw.(c)(r), tauCols);
    e = struct('dof', strtrim(raw.primary_dof{r}), 'anglesDeg', angles, ...
        'tauNm', taus);
    if isKey(curves, idn)
        curves(idn) = [curves(idn), e];
    else
        curves(idn) = e;
    end
end
end

function picked = orderCols_(cols, prefix)
idx = [];
for c = 1:numel(cols)
    if strncmpi(cols{c}, prefix, numel(prefix))
        n = str2double(strtrim(cols{c}(numel(prefix) + 1:end)));
        if ~isnan(n)
            idx(end + 1, :) = [n, c];   %#ok<AGROW>
        end
    end
end
idx = sortrows(idx, 1);
picked = cols(idx(:, 2));
end
