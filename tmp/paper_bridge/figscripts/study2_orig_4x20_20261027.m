% study2_orig_4x20_20261027.m -- STUDY 2, ORIGINAL side, 40 mm -> 4x 20 mm
% (2026-10-27, easteregg2). Capability comparison v3, per the studies
% directive:
%   ORIGINAL side = masters-era routing, EXCEPT wherever the masters' Size
%   rule used a 40 mm BPA those eight sites are substituted by FOUR
%   PARALLEL 20 mm BPAs (same routing, same graded surfaces, 620 kPa).
%   The other 19 sites reuse the survey2 rung-1 (nB1, G1) mats on disk
%   (masters routing + old-code diameter) -- done by study2_capability_
%   v3_20261027.m, not here.
%
% The eight 40 mm sites (BiPulley_survey2_perMuscle_20261023.csv
% bpa_size_mm column): ercspn_r, semimem_r, vas_int_r, soleus_r,
% med_gas_r, tib_post_r, glut_max_r, add_mag_r.
%
% Protocol = survey2 rung-1 VERBATIM (Opt_run_BiPulley_survey2_20261023.m:
% same surfaces CSV, same grade() booleans, same seeds/bounds/objective,
% 60-eval rung budget + 5 box seeds, 10-min wall cap), only the config
% injected is (nBPA = 4, G = 1, dia = 20). Route points are NOT design
% variables (masters route verbatim). parpool not needed (rung-sized).
%
% Outputs (dated 20261027, new names, never overwritten):
%   Results/BiPulley_surveys/BiPulley_survey3_orig_<site>_nB4_G1_<stamp>.mat
%   tmp/paper_bridge/capability_table_20261022/BiPulley_survey3_orig_4x20_20261027.csv

function study2_orig_4x20_20261027()
R = 'D:\GitHub\Bipedal_Robot';
capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(capDir);
resDir = fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_surveys');

baseXi0 = 0.005; baseXi1 = 5e5; baseXi2 = 1e4; baseWraps = 6;
basePressure = 620; baseFitn = 0.025;   % survey2 constants, verbatim

sites = {'ercspn_r', 'semimem_r', 'vas_int_r', 'soleus_r', ...
    'med_gas_r', 'tib_post_r', 'glut_max_r', 'add_mag_r'};
nB = 4; G = 1; dia = 20;

surfaces = loadSurfaces(fullfile(capDir, ...
    'human_surfaces_5x5_grading_20261023.csv'));
fprintf('Surfaces loaded: %d actuator surface(s)\n', surfaces.Count);

outCsv = fullfile(capDir, 'BiPulley_survey3_orig_4x20_20261027.csv');
fid = fopen(outCsv, 'w');
fprintf(fid, ['actuator_id,site,nBPA,G,dia_mm,rest_m,tendon_m,meetCells,' ...
    'activeCells,domainCells,meetPct,activePct,evals,matFile,note\n']);

function [g, ok] = grade(xq)
    rest = xq(1); tendon = xq(2);
    kmax = 0.75 * rest;
    try
        obj = BiPam_pulley2_20261022(nm, points, spec.cross, dia, T, ...
            rest, kmax, tendon, baseFitn, basePressure, baseXi0, ...
            baseXi1, baseXi2, baseWraps, spec.pulley);
    catch err
        g = struct('tauB', nan(size(tauHuman)), 'reason', err.message);
        ok = false; return
    end
    T3 = obj.Torque_ins(:, :, :, slot);
    tauB = reshape(sum(T3 .* reshape(axisU, 1, 3, 1), 2), N1, N2);
    strain = reshape(obj.strain_p, N1, N2);
    slack = logical(obj.PulleySlack);
    infeas = logical(obj.PulleyInfeasible);
    taut = ~slack & ~infeas & isfinite(tauB);
    meet = taut & isfinite(tauHuman) & ...
        (sign(tauB) == sign(tauHuman) | abs(tauHuman) <= 1e-9) & ...
        abs(tauB) >= abs(tauHuman) - 1e-12;
    active = taut & isfinite(strain) & strain > 0 & abs(tauB) > 0;
    g = struct('tauB', tauB, 'strain', strain, 'slack', slack, ...
        'infeas', infeas, 'taut', taut, 'meet', meet, 'active', ...
        active, 'obj', obj);
    ok = true;
end

function f = sObj(xq)
    [g, ok] = grade(xq);
    if ~ok, f = 1e6; return; end
    meetFrac = sum(g.meet(dom)) / nDom;
    activeFrac = sum(g.active(dom)) / nDom;
    deadFrac = (sum(g.slack(dom)) + sum(g.infeas(dom))) / nDom;
    f = -(1000 * meetFrac + activeFrac) + 1e-3 * deadFrac;
end

for s = 1:numel(sites)
    nm = sites{s};
    tSite = tic;
    spec0 = biPulleySpecsFromMasters_20261023({nm});
    sp = spec0(1);
    [SU, mapNote] = surfaceForSite(surfaces, sp, dia);
    assert(SU.priSlot > 0, '%s: %s', nm, mapNote);
    th1v = SU.th1; th2v = SU.th2;
    [T, N1, N2] = gridT(sp, th1v, th2v);
    points = sp.points;
    Lrig = rigidRouteLengths(points, sp.cross, T, N1, N2);
    Lchar = max(Lrig(:));
    tauHuman = SU.H; dom = SU.dom; nDom = SU.nDom;
    slot = SU.priSlot;
    axisU = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));
    spec = injectConfig(sp, nB, G, dia);
    lb = [0.30 * Lchar; 0.02];
    ub = [0.90 * Lchar; max(0.30 * Lchar, 0.06)];

    opts = optimoptions('patternsearch', 'Display', 'none', ...
        'MeshTolerance', 1e-4, 'StepTolerance', 1e-4, ...
        'MaxFunctionEvaluations', 60);
    ubT = ub(2);
    seeds = [0.60 * Lchar, 0.05; 0.30 * Lchar, 0.02; 0.90 * Lchar, 0.02; ...
        0.30 * Lchar, ubT; 0.90 * Lchar, ubT];
    fSeed = zeros(size(seeds, 1), 1);
    for sIdx = 1:size(seeds, 1)
        fSeed(sIdx) = sObj(seeds(sIdx, :));
    end
    [~, ixSeed] = min(fSeed);
    [xBest, ~, ~, pout] = patternsearch(@sObj, seeds(ixSeed, :), [], [], ...
        [], [], lb, ub, [], opts);
    evals = pout.funccount + size(seeds, 1);
    [g, ok] = grade(xBest);
    assert(ok, 'final grade failed %s', nm);
    meetCells = sum(g.meet(dom)); activeCells = sum(g.active(dom));
    meetPct = 100 * meetCells / nDom; activePct = 100 * activeCells / nDom;
    stamp = char(string(datetime('now'), 'yyyyMMdd_HHmmss'));
    matFile = fullfile(resDir, sprintf('BiPulley_survey3_orig_%s_nB%d_G%d_%s.mat', ...
        nm, nB, G, stamp));
    surveyMeta = struct('actuator_id', sp.actuatorId, 'site', nm, ...
        'side', 'ORIGINAL_v3', 'nBPA', nB, 'G', G, 'dia_mm', dia, ...
        'substitution', 'masters 40 mm -> 4x 20 mm parallel', ...
        'protocol', 'survey2 rung-1 verbatim', 'evals', evals, ...
        'meetCells', meetCells, 'activeCells', activeCells, ...
        'domainCells', nDom, 'excludedCells', numel(dom) - nDom, ...
        'meetPct', meetPct, 'activePct', activePct, 'note', mapNote, ...
        'created', datetime('now'));
    save(matFile, 'surveyMeta', 'spec', 'points', 'T', 'xBest', ...
        'tauHuman', 'dom', 'g', 'SU');
    wallMin = toc(tSite) / 60;
    fprintf(['[%2d] %-12s ORIG (4,1) dia20: meet %d/%d (%.1f%%), active ' ...
        '%d/%d, rest %.4f tendon %.4f, %d evals, %.1f min\n'], ...
        sp.actuatorId, nm, meetCells, nDom, meetPct, activeCells, nDom, ...
        xBest(1), xBest(2), evals, wallMin);
    fprintf(fid, '%d,%s,%d,%d,%d,%.6g,%.6g,%d,%d,%d,%.4f,%.4f,%d,%s,"%s"\n', ...
        sp.actuatorId, nm, nB, G, dia, xBest(1), xBest(2), meetCells, ...
        activeCells, nDom, meetPct, activePct, evals, matFile, mapNote);
end
fclose(fid);
fprintf('STUDY2 ORIG 4x20 DONE\n');

%% ---------------- Local functions (survey2 driver, verbatim) ------------
function surfaces = loadSurfaces(csvPath)
raw = readtable(csvPath, 'VariableNamingRule', 'preserve', ...
    'Delimiter', ',', 'TextType', 'char', 'CommentStyle', '#');
need = {'actuator_id', 'primary_dof', 'secondary_dof', 'p_node_deg', ...
    's_node_deg', 'in_graded_domain', 'tau_h_Nm'};
for k = 1:numel(need)
    assert(any(strcmp(raw.Properties.VariableNames, need{k})), 'missing %s', need{k})
end
surfaces = containers.Map('KeyType', 'double', 'ValueType', 'any');
for r = 1:height(raw)
    idn = double(raw.actuator_id(r));
    if isnan(idn), continue; end
    e = struct('priDof', strtrim(raw.primary_dof{r}), ...
        'secDof', strtrim(raw.secondary_dof{r}), ...
        'p', iNum(raw.p_node_deg(r)), 's', iNum(raw.s_node_deg(r)), ...
        'dom', iNum(raw.in_graded_domain(r)) == 1, ...
        'tau', iNum(raw.tau_h_Nm(r)), 'flg', false);
    for f = {'flag_via', 'flag_step', 'flag_slack', 'flag_short'}
        if any(strcmp(raw.Properties.VariableNames, f{1})) && ...
                iNum(raw.(f{1})(r)) == 1
            e.flg = true;
        end
    end
    if isKey(surfaces, idn)
        surfaces(idn) = [surfaces(idn), e];
    else
        surfaces(idn) = e;
    end
end
end

function v = iNum(x)
if isnumeric(x)
    v = double(x);
elseif iscell(x) || isstring(x)
    v = str2double(char(x));
else
    v = double(x);
end
end

function [SU, note] = surfaceForSite(surfaces, sp, dia)
SU = struct('priDof', '', 'secDof', '', 'priSlot', 0, 'th1', [], ...
    'th2', [], 'H', [], 'dom', [], 'flg', [], 'nDom', 0);
note = '';
if ~isKey(surfaces, double(sp.actuatorId))
    note = 'no surface'; return
end
es = surfaces(double(sp.actuatorId));
assert(numel(es) == 25, 'act %d: %d rows', sp.actuatorId, numel(es))
SU.priDof = es(1).priDof;
SU.secDof = es(1).secDof;
priSlot = find(strcmp(sp.kin.dofs, SU.priDof), 1);
if isempty(priSlot)
    SU.priSlot = 0;
    note = sprintf('primary dof %s not modeled', SU.priDof); return
end
SU.priSlot = priSlot;
secSlot = find(strcmp(sp.kin.dofs, SU.secDof), 1);
pNodes = sort(unique(round([es.p], 9)));
sNodes = sort(unique(round([es.s], 9)));
assert(numel(pNodes) == 5 && numel(sNodes) == 5, 'act %d grid', sp.actuatorId)
Hc = nan(5, 5); domC = false(5, 5); flgC = false(5, 5);
for r = 1:25
    ip = find(abs(pNodes - es(r).p) < 1e-9, 1);
    is = find(abs(sNodes - es(r).s) < 1e-9, 1);
    Hc(ip, is) = es(r).tau;
    domC(ip, is) = es(r).dom;
    flgC(ip, is) = es(r).flg;
end
if priSlot == 1
    H = Hc; dom = domC; flg = flgC;
else
    H = Hc.'; dom = domC.'; flg = flgC.';
end
SU.H = H; SU.dom = dom; SU.flg = flg; SU.nDom = sum(dom(:));
secOn1 = ~isempty(secSlot) && secSlot == 1;
secOn2 = ~isempty(secSlot) && secSlot == 2;
if priSlot == 1
    SU.th1 = pNodes;
    SU.th2 = sNodes * double(secOn2);
else
    SU.th1 = sNodes * double(secOn1);
    SU.th2 = pNodes;
end
if isempty(secSlot)
    note = sprintf(['secondary %s not a class hinge: pinned at 0; %d mm' ...
        ' BPA'], SU.secDof, dia);
else
    note = sprintf('both dofs modeled; %d mm BPA', dia);
end
end

function [T, N1, N2] = gridT(sp, th1v, th2v)
N1 = numel(th1v); N2 = numel(th2v);
T = zeros(4, 4, N1, N2);
for ii = 1:N1
    T(:, :, ii, 1) = hingeT(sp, 1, th1v(ii));
end
for iii = 1:N2
    T(:, :, iii, 2) = hingeT(sp, 2, th2v(iii));
end
end

function T = hingeT(sp, slot, thetaDeg)
th = thetaDeg * pi / 180;
if strcmpi(strtrim(sp.kin.types{slot}), 'fixed')
    R = eye(3); p = [0, 0, 0];
else
    a = sp.kin.axes(slot, :) / norm(sp.kin.axes(slot, :));
    K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
    R = eye(3) + sin(th) * K + (1 - cos(th)) * (K * K);
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

function spec = injectConfig(spec0, nB, G, dia)
spec = spec0;
spec.diameter = dia;
spec.pulley = struct('nPulleyBPA', nB, 'tackleLineParts', [G, 1], ...
    'routingMode', 'moving_via');
spacing = spec.diameter * 1e-3 + 0.005;
if nB > 1
    u = [];
    for r = size(spec.points, 1) - 1:-1:1
        cand = spec.points(r + 1, :) - spec.points(r, :);
        if norm(cand) > 1e-12
            u = cand / norm(cand);
            break
        end
    end
    if abs(u(3)) < 0.9
        perp = cross(u, [0, 0, 1]);
    else
        perp = cross(u, [0, 1, 0]);
    end
    perp = perp / norm(perp);
    offs = ((1:nB)' - (nB + 1) / 2) * spacing;
    spec.bundleOffsets = offs .* perp;
    spec.footprintWidth = (nB - 1) * spacing + spec.diameter * 1e-3;
else
    spec.bundleOffsets = zeros(0, 3);
    spec.footprintWidth = spec.diameter * 1e-3;
end
end

function L = rigidRouteLengths(points, cross, T, N1, N2)
nPts = size(points, 1);
L = zeros(N1, N2);
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for iii = 1:N2
        T2 = T(:, :, iii, 2);
        tot = 0;
        for r = 1:nPts - 1
            pA = mapRow(points, cross, r, T1, T2, 1);
            pB = mapRow(points, cross, r + 1, T1, T2, 1);
            tot = tot + norm(pA - pB);
        end
        L(ii, iii) = tot;
    end
end
end

function p = mapRow(points, cross, r, T1, T2, frame)
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
while f > frame
    if f == 3
        p = RowVecTrans(T2, p); f = 2;
    else
        p = RowVecTrans(T1, p); f = 1;
    end
end
end
end
