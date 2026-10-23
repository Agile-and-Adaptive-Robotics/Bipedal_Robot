% Sanity Bi Pulley Specs 27 20261023
% Per-actuator BUILD + TORQUE-CURVE SANITY PASS (NO SOLVES) for
% biPulleySpecsFromMasters_20261023 -- the extended masters builder
% (FIX_SPEC_20261023.md section 3). "No solves" = no design optimization:
% every actuator is built, its two-hinge 5x5 grid (the machinery's grading
% convention) is transformed, the BiPam_pulley2_20261022 class (machinery
% of record) is constructed at the survey's default sizing (rest = 0.6*Lchar,
% tendon = 0.05 m, 620 kPa, xi 0.005/5e5/1e4, wraps 6, config
% (1,1,'moving_via')), and the primary-slot torque curve is inspected. The
% class's internal per-cell equilibrium root find is part of any torque
% evaluation (the survey's grade()), not a design solve.
%
% CHECKS per actuator:
%   1. BUILD: the builder returns a spec (it also self-asserts that the
%      primary dof is one of the two modeled hinges).
%   2. NEUTRAL-POSE LENGTH vs THE MASTERS' OWN FORMULAS: an independent
%      reference path recomputes the route length at theta = 0 straight
%      from the raw AttachPoints_Robot.mat variables plus the masters'
%      published pivot offsets (per route family, derived from the
%      masters' VecTrans conventions); the class's MuscleLength at the
%      neutral cell must agree to 1e-9 m. This validates the frame/pivot/
%      fold bookkeeping end to end.
%   3. TORQUE-CURVE SANITY: on the 5x5 grid, across the survey's 5-seed
%      size box (rest in [0.30, 0.90]*Lchar, tendon in [0.02, max(0.30*
%      Lchar, 0.06)]; plain evaluations, no optimizer), the primary-axis
%      torque must have >= 1 taut finite cell with |tau| > 0.01 N*m at SOME
%      seed. A muscle with an empty taut window at every seed is still
%      buildable/gradeable (the honest slack/infeasible flags are the
%      machinery's own output, cf. med_gas in the 10-22 survey) but is
%      reported loudly as capability-flagged for the ladder stage.
%   4. SIGN vs HUMAN CURVE (reported, not pass/fail): sign agreement over
%      taut cells where the human curve exceeds 1 N*m in magnitude, with
%      human_torque_curves_20261022.csv at matched primary angles and the
%      secondary pinned at its default; < 50% prints a WARNING (a
%      systematic 0% is the flag for an axis-sign mismatch).
% Output: printed table + sanity_biPulleySpecs27_20261023.csv (dated new
% artifact). Exit: error if any actuator fails a hard check (1-3).

function results = sanity_biPulleySpecs27_20261023()

scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));   % festo4/Spr/etc + class
addpath(scriptDir);                                    % the dated builder

%% Survey-default physical constants (Opt_run_BiPulley_survey_20261022)
pres = 620; xi0 = 0.005; xi1 = 5e5; xi2 = 1e4; wraps = 6; fitn = 0.025;

%% Build ALL 27
specs = biPulleySpecsFromMasters_20261023();
nA = numel(specs);
fprintf('Built %d specs from biPulleySpecsFromMasters_20261023\n', nA);
assert(nA == 27, 'expected 27 specs, got %d', nA);

%% Human curves (sign reference only)
curves = loadHumanCurves_(fullfile(scriptDir, 'human_torque_curves_20261022.csv'));

%% Masters-side neutral-pose offsets per route family. Derived from the
%% masters' VecTrans chain at theta = 0 (every hinge is then a pure
%% translation by its pivot): frame2 rows shift by the hinge-1 pivot (plus
%% any builder fold), frame3 rows by hinge-2 pivot + hinge-1 pivot (plus
%% folds). Independent of spec.kin -- raw data + published constants only.
pHip = [-0.0707, -0.0661, 0.0835];
pAnkle = [0, -0.43, 0];
pSub = [-0.04877, -0.04195, 0.00792];
pMtp = [0.1788, -0.002, 0.00108];
pBack = [-0.1007, 0.0815, 0];
knee0 = [interp1(specs(1).kin.kneeFit.thetaX, specs(1).kin.kneeFit.x, 0), ...
    interp1(specs(1).kin.kneeFit.thetaY, specs(1).kin.kneeFit.y, 0), 0];
fprintf('Rolling-knee neutral center (femur coords): [%.4f %.4f %.4f] m\n', ...
    knee0);

fam = containers.Map('KeyType', 'char', 'ValueType', 'any');
fam('pelvis|torsoA|torsoB')    = struct('f2', [],           'f3', pBack);
fam('pelvis|femur|tibia')      = struct('f2', pHip,        'f3', knee0 + pHip);
fam('femur|tibia|tibia_rigid') = struct('f2', knee0,       'f3', knee0);
fam('femur|tibia|calcn')       = struct('f2', [],          'f3', pAnkle + knee0);
fam('tibia|calcnA|calcnSub')   = struct('f2', [],          'f3', pSub + pAnkle);
fam('tibia|calcnA|toes')       = struct('f2', pSub + pAnkle, ...
                                       'f3', pMtp + pSub + pAnkle);
fam('pelvis|femurA|femurB')    = struct('f2', [],          'f3', pHip);
fam('torso|pelvisShift|femur') = struct('f2', -pBack,      'f3', pHip - pBack);

%% Per-actuator pass
rows = cell(nA, 1);
nFail = 0;
fprintf(['\n%-3s %-12s %-17s %-16s %-16s %8s %9s %9s %10s %6s %5s\n'], ...
    'id', 'spec', 'primary', 'hinge1', 'hinge2', 'Lchar', 'tauPkNm', ...
    'neutral', 'ref', 'sign%', 'ok');
for a = 1:nA
    sp = specs(a);
    key = strjoin(sp.frames, '|');
    if ~isKey(fam, key)
        error('sanity:noFamily', 'no neutral-reference family for %s', key)
    end
    FR = fam(key);

    % -- grid: 5x5, both hinges swept over their ranges (machinery conv.) --
    nGrid = 5;
    th1 = linspace(sp.kin.thetaRangesDeg(1, 1), sp.kin.thetaRangesDeg(1, 2), nGrid);
    th2 = linspace(sp.kin.thetaRangesDeg(2, 1), sp.kin.thetaRangesDeg(2, 2), nGrid);
    T = gridT_(sp, th1, th2, nGrid);

    % -- rigid lengths -> characteristic length (survey convention) -------
    points = sp.points;   % folds already applied by the builder
    Lrig = rigidLens_(points, sp.cross, T, nGrid);
    Lchar = max(Lrig(:));

    % -- neutral reference (masters formulas, raw data, independent path) -
    raw = points - sp.kin.rowFolds;    % undo the builder folds
    c1 = sp.cross(1);  c2 = sp.cross(2);
    nRows = size(raw, 1);
    fOf = ones(nRows, 1);
    if c1 == c2
        fOf(c1:end) = 3;
    else
        fOf(c1:c2 - 1) = 2;
        fOf(c2:end) = 3;
    end
    ref = raw;
    if ~isempty(FR.f2)
        ref(fOf == 2, :) = ref(fOf == 2, :) + FR.f2;
    end
    ref(fOf == 3, :) = ref(fOf == 3, :) + FR.f3;
    Lref = 0;
    for r = 1:nRows - 1
        Lref = Lref + norm(ref(r + 1, :) - ref(r, :));
    end

    % -- construct the class (torque evaluation; no design solve) ---------
    % Seed-box scan (the survey's own 5-point box multi-start convention,
    % Opt_run_BiPulley_survey_20261022 runRung): 5 plain evaluations of the
    % class across the legal size box -- no optimizer, no ladder. A taut
    % window anywhere in the box passes; the best seed (most taut cells,
    % then peak torque) is the reported curve.
    ok = true;  failReason = '';
    neutralLen = NaN; neutralErr = NaN; tauPk = NaN; nTaut = 0; ...
        nFinite = 0; signPct = NaN;
    bestSeed = [NaN, NaN];
    ubT = max(0.30 * Lchar, 0.06);
    seeds = [0.60 * Lchar, 0.05; ...
             0.30 * Lchar, 0.02; ...
             0.90 * Lchar, 0.02; ...
             0.30 * Lchar, ubT; ...
             0.90 * Lchar, ubT];
    try
        slot = sp.kin.primarySlot;
        axisU = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));
        bestScore = -inf;
        tauBest = [];
        tautBest = [];
        for sIdx = 1:size(seeds, 1)
            restS = seeds(sIdx, 1);
            tendonS = seeds(sIdx, 2);
            kmaxS = 0.75 * restS;
            obj = BiPam_pulley2_20261022(sp.name, points, sp.cross, 20, ...
                T, restS, kmaxS, tendonS, fitn, pres, xi0, xi1, xi2, ...
                wraps, sp.pulley);
            if sIdx == 1
                % neutral-pose length check (rigid MuscleLength is sizing-
                % independent; theta = 0 on both hinges; dim 3 duplicated
                % per the class's T layout contract)
                T0 = zeros(4, 4, 2, 2);
                T0(:, :, 1, 1) = hingeT_(sp, 1, 0);
                T0(:, :, 2, 1) = hingeT_(sp, 1, 0);
                T0(:, :, 1, 2) = hingeT_(sp, 2, 0);
                T0(:, :, 2, 2) = hingeT_(sp, 2, 0);
                obj0 = BiPam_pulley2_20261022(sp.name, points, sp.cross, ...
                    20, T0, restS, kmaxS, tendonS, fitn, pres, xi0, xi1, ...
                    xi2, wraps, sp.pulley);
                neutralLen = obj0.MuscleLength(1, 1);
                neutralErr = abs(neutralLen - Lref);
            end
            T3 = obj.Torque_ins(:, :, :, slot);
            tauB = reshape(sum(T3 .* reshape(axisU, 1, 3, 1), 2), ...
                nGrid, nGrid);
            taut = ~logical(obj.PulleySlack) & ...
                ~logical(obj.PulleyInfeasible) & isfinite(tauB);
            nT = sum(taut(:));
            pk = 0;
            if any(taut(:))
                pk = max(abs(tauB(taut)));
            end
            score = nT + pk;   % most taut cells first, then peak torque
            if score > bestScore
                bestScore = score;
                bestSeed = seeds(sIdx, :);
                nTaut = nT;
                tauPk = pk;
                tauBest = tauB;
                tautBest = taut;
                nFinite = sum(isfinite(tauB(:)));
            end
        end

        % sign vs human curve on the primary sweep, secondary pinned at 0
        % (row orientation forced: slot-1 slices come out as columns)
        HC = curvesFor_(curves, sp.actuatorId, sp.kin.primaryDof);
        if ~isempty(HC) && ~isempty(tauBest)
            if slot == 1
                j0 = iNearest0_(th2);
                tauSlice = reshape(tauBest(:, j0), 1, []);
                tautSlice = reshape(tautBest(:, j0), 1, []);
                thP = th1;
            else
                i0 = iNearest0_(th1);
                tauSlice = reshape(tauBest(i0, :), 1, []);
                tautSlice = reshape(tautBest(i0, :), 1, []);
                thP = th2;
            end
            tauH = interp1(HC.anglesDeg, HC.tauNm, thP, 'linear', NaN);
            m = tautSlice & isfinite(tauH) & abs(tauH) > 1 & ...
                abs(tauSlice) > 1e-9;
            if any(m)
                signPct = 100 * mean(sign(tauSlice(m)) == sign(tauH(m)));
                if signPct < 50
                    fprintf(['  WARNING %s: sign agreement %.0f%% (%d ' ...
                        'cells) -- check axis convention\n'], sp.name, ...
                        signPct, sum(m));
                end
            end
        end
    catch err
        ok = false;
        failReason = sprintf('runtime: %s (%s:%d)', err.message, ...
            err.stack(1).name, err.stack(1).line);
    end

    % hard checks
    if ok && neutralErr > 1e-9
        ok = false;
        failReason = sprintf('neutral length mismatch %.3g m', neutralErr);
    end
    if ok && (nTaut < 1 || tauPk <= 0.01)
        ok = false;
        failReason = sprintf('no taut torque (taut %d, peak %.4g N*m)', ...
            nTaut, tauPk);
    end
    if ~ok
        nFail = nFail + 1;
    end

    h1 = sp.kin.dofs{1};  h2 = sp.kin.dofs{2};
    if isempty(h1); h1 = '(rigid)'; end
    if isempty(h2); h2 = '(rigid)'; end
    fprintf(['%-3d %-12s %-17s %-16s %-16s %8.4f %9.3f %9.5f %9.5f ' ...
        '%6.0f %5s\n'], sp.actuatorId, sp.name, sp.kin.primaryDof, h1, h2, ...
        Lchar, tauPk, neutralLen, Lref, signPct, iIf(ok, 'PASS', 'FAIL'));
    if ~ok
        fprintf('      FAIL: %s\n', failReason);
    end

    rows{a} = { ...
        sp.actuatorId, sp.masterName, sp.name, sp.kin.primaryDof, ...
        sp.kin.primarySlot, iIf(ok, 'PASS', 'FAIL'), failReason, ...
        h1, mat2str(sp.kin.axes(1, :)), ...
        sprintf('[%g %g]', sp.kin.thetaRangesDeg(1, :)), ...
        h2, mat2str(sp.kin.axes(2, :)), ...
        sprintf('[%g %g]', sp.kin.thetaRangesDeg(2, :)), ...
        sprintf('[%d %d]', sp.cross(1), sp.cross(2)), ...
        nRows, Lchar, bestSeed(1), bestSeed(2), neutralLen, Lref, ...
        neutralErr, nTaut, nFinite, tauPk, signPct};
end

%% CSV (dated new artifact)
hdr = {'actuator_id', 'master_name', 'spec', 'primary_dof', 'primary_slot', ...
    'status', 'fail_reason', 'hinge1_dof', 'hinge1_axis', ...
    'hinge1_range_deg', 'hinge2_dof', 'hinge2_axis', 'hinge2_range_deg', ...
    'cross', 'n_points', 'Lchar_m', 'best_seed_rest_m', ...
    'best_seed_tendon_m', 'neutral_len_m', 'neutral_ref_m', ...
    'neutral_abs_err_m', 'taut_cells_of_25', 'finite_tau_cells', ...
    'peak_abs_tau_Nm', 'sign_agree_pct'};
csvFile = fullfile(scriptDir, 'sanity_biPulleySpecs27_20261023.csv');
fid = fopen(csvFile, 'w');
fprintf(fid, ['# sanity_biPulleySpecs27_20261023 - per-actuator build + ' ...
    'torque-curve sanity pass (no solves)\n']);
fprintf(fid, ['# builder biPulleySpecsFromMasters_20261023; class ' ...
    'BiPam_pulley2_20261022; 5x5 grid; 620 kPa; taut window scanned over ' ...
    'the survey 5-seed size box (no optimizer)\n']);
fprintf(fid, '%s\n', strjoin(hdr, ','));
for a = 1:nA
    fprintf(fid, '%s\n', joinRow_(rows{a}));
end
fclose(fid);

fprintf('\nSANITY PASS COMPLETE: %d/27 PASS, %d FAIL\n', 27 - nFail, nFail);
fprintf('CSV: %s\n', csvFile);
if nFail > 0
    error('sanity_biPulleySpecs27_20261023:failures', ...
        '%d actuator(s) failed the sanity pass', nFail)
end
results = rows;
end

%% =====================================================================
%% Local helpers
%% =====================================================================
function T = gridT_(sp, th1, th2, nGrid)
% 4x4xnGridxnGrid: T(:,:,ii,1) = hinge 1 at th1(ii); T(:,:,iii,2) =
% hinge 2 at th2(iii) (both hinge angle grids indexed along dim 3, dim 4
% selects the joint -- the class's T layout contract).
T = zeros(4, 4, nGrid, nGrid);
for ii = 1:nGrid
    T(:, :, ii, 1) = hingeT_(sp, 1, th1(ii));
end
for iii = 1:nGrid
    T(:, :, iii, 2) = hingeT_(sp, 2, th2(iii));
end
end

function T = hingeT_(sp, slot, thetaDeg)
th = thetaDeg * pi / 180;
if strcmpi(strtrim(sp.kin.types{slot}), 'fixed')
    R = eye(3);
    p = [0, 0, 0];
else
    R = rotAxis_(sp.kin.axes(slot, :), th);
    p = sp.kin.pivots(slot, :);
    if sp.kin.kneeSlot == slot && isfield(sp.kin, 'kneeFit') && ...
            ~isempty(sp.kin.kneeFit.thetaY)
        thK = sp.kin.kneeFit.thetaY;
        xK = interp1(sp.kin.kneeFit.thetaX, sp.kin.kneeFit.x, thK, ...
            'linear', 'extrap');
        p = [interp1(thK, xK, th, 'linear', 'extrap'), ...
            interp1(thK, sp.kin.kneeFit.y, th, 'linear', 'extrap'), p(3)];
    end
end
T = [R, p(:); 0, 0, 0, 1];
end

function R = rotAxis_(axis, theta)
a = axis(:) / norm(axis);
K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end

function L = rigidLens_(points, cross, T, nGrid)
nPts = size(points, 1);
L = zeros(nGrid, nGrid);
for ii = 1:nGrid
    T1 = T(:, :, ii, 1);
    for iii = 1:nGrid
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
% Route row r expressed in the proximal (frame-1) grid cell (the class's
% mapVectorToFrame chain, replicated for the rigid length reference).
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
        p = rowTrans_(T2, p);
        f = 2;
    else
        p = rowTrans_(T1, p);
        f = 1;
    end
end
end

function v = rowTrans_(T, v)
r = [v, 1] * T.';
v = r(1:3);
end

function ix = iNearest0_(v)
[~, ix] = min(abs(v));
end

function s = iIf(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end

function s = joinRow_(cellsRow)
parts = cell(size(cellsRow));
for i = 1:numel(cellsRow)
    v = cellsRow{i};
    if isnumeric(v)
        if isscalar(v) && mod(v, 1) == 0
            parts{i} = sprintf('%d', v);
        else
            parts{i} = sprintf('%.10g', v);
        end
    elseif ischar(v)
        parts{i} = ['"', strrep(v, '"', ''), '"'];
    else
        parts{i} = sprintf('"%s"', char(v));
    end
end
s = strjoin(parts, ',');
end

function curves = loadHumanCurves_(csvPath)
curves = [];
if ~isfile(csvPath)
    warning('sanity:noHumanCsv', 'human CSV not found: %s', csvPath)
    return
end
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

function HC = curvesFor_(curves, actId, dof)
HC = [];
if isempty(curves) || ~isKey(curves, double(actId))
    return
end
es = curves(double(actId));
for e = 1:numel(es)
    if strcmp(es(e).dof, dof)
        HC = es(e);
        return
    end
end
end
