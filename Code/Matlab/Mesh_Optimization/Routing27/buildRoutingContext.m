function ctx = buildRoutingContext(spec, chain, targetsDir, boneDir, ...
    finalizedRoutes, opts)
% Build Routing Context
% Campaign: routing27 (2026-10-03)
% Assembles the per-actuator optimization context: the two-slab transform
% grid on the PRIMARY DOF (the other crossing holds its default pose),
% the stock-human target curve interpolated onto that grid (pchip, the
% objective_KneeExt20mm convention), the Xi-of-record constants, the
% per-cell bone clouds expressed in the proximal frame (the no-penetration
% constraint input), the lateral-footprint cap (bone-cloud z extent + BPA
% diameter + 20 mm), and the design-variable seed/bounds.
%
% Everything map-based is resolved HERE so the evaluation context is
% plain arrays/cells (serializable to parallel workers; containers.Map
% does not serialize).
%
% Xi of record (flexor 2brk CV pick 77 with its locked pair; scout assets
% list 2026-10-03): Xi0 +3.94 mm, Xi1 3.998e4 N/m, Xi2 1.4734e4 N/m,
% applied to every actuator mount (documented choice: one advisor-
% consistent triple across configurations; Xi1 sits on a flat likelihood
% valley per AGENTS, and Opt_run_BiPulley's hand-tune neighborhood is
% superseded by the pick of record).

if nargin < 6
    opts = struct();
end
Nangles = getOrDefault(opts, 'Nangles', 13);
requiredMargin = getOrDefault(opts, 'requiredMargin', 0.05);
pressure = getOrDefault(opts, 'pressure', 620);
wraps = getOrDefault(opts, 'wraps', 6);
fitting = getOrDefault(opts, 'fitting', 0.025);
bpaDia_mm = getOrDefault(opts, 'diameter_mm', 20);
boneMargin = getOrDefault(opts, 'boneMargin_mm', 5) * 1e-3;
spacingMargin = getOrDefault(opts, 'spacingMargin_mm', 5) * 1e-3;
viaDelta = getOrDefault(opts, 'viaDelta_mm', 30) * 1e-3;
endDelta = getOrDefault(opts, 'endDelta_mm', 30) * 1e-3;

ctx = struct();
ctx.spec = spec;
ctx.label = sprintf('a%02d_%s', spec.actuatorId, ...
    strrep(strrep(strrep(strrep(spec.masterName, ' ', '_'), ',', ''), ...
    '(', ''), ')', ''));
ctx.diameter = bpaDia_mm;
ctx.bpaRadius = 0.0385 / 2;          % fully inflated 20 mm BPA radius
ctx.pressure = pressure;
ctx.wraps = wraps;
ctx.fitting = fitting;
ctx.Xi0 = 0.00394;
ctx.Xi1 = 3.998e4;
ctx.Xi2 = 1.4734e4;
ctx.requiredMargin = requiredMargin;

%% Primary crossing split and the transform grid
k = spec.targetCrossing;
ctx.targetCrossing = k;
ctx.axisHat = spec.kin.axes(k, :) / norm(spec.kin.axes(k, :));
rPrim = spec.kin.thetaRangesDeg(k, :);
ctx.thetaPrimaryDeg = linspace(rPrim(1), rPrim(2), Nangles)';
otherSlot = 3 - k;
thetaOtherDeg = 0;                    % frozen crossing at its default

% BiPamData slab convention (verified against the committed batch): BOTH
% slabs are indexed along dim3 (the angle axis); dim4 is the slab label,
% so T(:,:,ii,1) is crossing-1's transform at angle ii and T(:,:,iii,2)
% crossing-2's at angle iii. When the primary DOF is crossing 1 the
% frozen joint-2 slab is replicated, so dim4 = 2 suffices (26 cells);
% when the primary is crossing 2 the class's second grid dimension must
% span the primary angles, so the grid is Nang x Nang.
if k == 1
    ctx.primaryDim = 3;
    nDim4 = 2;
    th1v = (pi / 180) * ctx.thetaPrimaryDeg(:);
    th2v = (pi / 180) * thetaOtherDeg * ones(Nangles, 1);
else
    ctx.primaryDim = 4;
    nDim4 = Nangles;
    th1v = (pi / 180) * thetaOtherDeg * ones(Nangles, 1);
    th2v = (pi / 180) * ctx.thetaPrimaryDeg(:);
end
T = zeros(4, 4, Nangles, nDim4);
for i = 1:Nangles
    p1 = spec.kin.pivots(1, :);
    if spec.kin.kneeSlot == 1 && rollingSpec(spec.kin)
        p1 = rollPivot(spec.kin, th1v(i), p1);
    end
    T(:, :, i, 1) = localHinge(spec.kin.types{1}, spec.kin.axes(1, :), ...
        p1, th1v(i));
    p2 = spec.kin.pivots(2, :);
    if spec.kin.kneeSlot == 2 && rollingSpec(spec.kin)
        p2 = rollPivot(spec.kin, th2v(i), p2);
    end
    T(:, :, i, 2) = localHinge(spec.kin.types{2}, spec.kin.axes(2, :), ...
        p2, th2v(i));
end
ctx.T = T;

%% Characteristic route length (sizing, the batch convention): the max
%% rigid route length over the primary sweep (frozen joint at default).
nPts = size(spec.points, 1);
Lrig = zeros(Nangles, 1);
for i1 = 1:Nangles
    tot = 0;
    for r = 1:nPts - 1
        pA = segMapRow(spec.points, spec.cross, r, T(:, :, i1, 1), ...
            T(:, :, i1, 2), 1);
        pB = segMapRow(spec.points, spec.cross, r + 1, ...
            T(:, :, i1, 1), T(:, :, i1, 2), 1);
        tot = tot + norm(pA - pB);
    end
    Lrig(i1) = tot;
end
ctx.Lchar = max(Lrig(:));

%% Human target: stock CSV -> grid (pchip interpolation, extrap)
csvName = sprintf('targets27_%02d_%s_%s_primary.csv', spec.actuatorId, ...
    safeLabel(spec.masterName), spec.primaryDof);
csvPath = fullfile(targetsDir, csvName);
if ~exist(csvPath, 'file')
    error('buildRoutingContext:noTarget', ...
        'Missing human target %s (run gen_gait2392_torque_targets.py)', ...
        csvPath)
end
H = readmatrix(csvPath, 'NumHeaderLines', 1);
ctx.humanAngleDeg = H(:, 1);
tauHumanSigned = H(:, end);            % tau_group_Nm is the last column
ctx.humanAbsGrid = abs(interp1(ctx.humanAngleDeg, abs(tauHumanSigned), ...
    ctx.thetaPrimaryDeg, 'pchip', 'extrap'));
ctx.torqueScale = max(1, max(ctx.humanAbsGrid));
ctx.humanPeakNm = max(ctx.humanAbsGrid);

%% Design variables: origin/end/via bounded deltas + rest + tendon
viaRows = [];
if nPts >= 4
    viaRows = (2:min(3, nPts - 1))';
end
ctx.movableViaRows = viaRows;
nD = 6 + 3 * numel(viaRows) + 2;
x0 = zeros(nD, 1);
lb = x0;
ub = x0;
lb(1:3) = -endDelta;  ub(1:3) = endDelta;    % origin row 1 (frame 1)
lb(4:6) = -endDelta;  ub(4:6) = endDelta;    % insertion row end (distal)
for j = 1:numel(viaRows)
    lb(6 + 3 * (j - 1) + (1:3)) = -viaDelta;
    ub(6 + 3 * (j - 1) + (1:3)) = viaDelta;
end
iRest = nD - 1;
iTend = nD;
x0(iRest) = 0.60 * ctx.Lchar;
x0(iTend) = 0.05;
lb(iRest) = 0.30 * ctx.Lchar;
ub(iRest) = 0.90 * ctx.Lchar;
lb(iTend) = 0.02;
ub(iTend) = max(0.30 * ctx.Lchar, 0.06);
ctx.x0 = x0;
ctx.lb = lb;
ctx.ub = ub;
ctx.iRest = iRest;
ctx.iTend = iTend;

%% Bone clouds: read once per body, transform per cell into frame 1
routeBodies = cellfun(@stripNote, spec.frames, 'UniformOutput', false);
proxBody = routeBodies{1};
distalBody = routeBodies{3};
order = {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r', 'toes_r'};
bodies = {};
if any(strcmp(routeBodies, 'torso'))
    bodies = [bodies, {'torso', 'pelvis'}];
end
i1 = find(strcmp(order, proxBody), 1);
i2 = find(strcmp(order, distalBody), 1);
if ~isempty(i1) && ~isempty(i2)
    bodies = [bodies, order(min(i1, i2):max(i1, i2))]; %#ok<AGROW>
elseif strcmp(proxBody, 'pelvis')
    bodies = [bodies, 'pelvis']; %#ok<AGROW>
end
bodies = unique(bodies, 'stable');

cloudFiles = struct( ...
    'pelvis', {{'Pelvis_R'}}, ...
    'torso', {{'Spine', 'Sacrum'}}, ...
    'femur_r', {{'Femur'}}, ...
    'tibia_r', {{'Tibia'}}, ...
    'talus_r', {{'Talus'}}, ...
    'calcn_r', {{'Calcaneus', 'Talus'}}, ...
    'toes_r', {{'Toes', 'Calcaneus'}});
bodyClouds = cell(1, numel(bodies));
maxZ = 0;
for b = 1:numel(bodies)
    fl = cloudFiles.(bodies{b});
    pts = zeros(0, 3);
    for f = 1:numel(fl)
        p = fullfile(boneDir, sprintf('%s_Mesh_Points.txt', fl{f}));
        if exist(p, 'file')
            pts = [pts; readmatrix(p)]; %#ok<AGROW>
        else
            warning('buildRoutingContext:noCloud', 'Missing cloud %s', p)
        end
    end
    bodyClouds{b} = pts;
    maxZ = max(maxZ, max(abs(pts(:, 3))));
end

primKey = jointEdgeKey(spec.kin.coordNames{k});
cellClouds = cell(Nangles, 1);
for a = 1:Nangles
    X = bodyTransforms(chain, ctx.thetaPrimaryDeg(a), primKey);
    Xf1inv = inv(X.(safeField(proxBody)));
    clouds1 = cell(1, numel(bodies));
    for b = 1:numel(bodies)
        M = Xf1inv * X.(safeField(bodies{b}));
        C = bodyClouds{b};
        clouds1{b} = (M(1:3, 1:3) * C.' + M(1:3, 4)).';
    end
    cellClouds{a} = clouds1;
end
ctx.cellClouds = cellClouds;
ctx.boneRequiredClear = ctx.bpaRadius + boneMargin;
ctx.latCap = maxZ + ctx.bpaRadius + 0.020;   % bone z extent + BPA + 20 mm
ctx.spacingRequired = (bpaDia_mm * 1e-3) + spacingMargin;
ctx.bundleSpacing = bpaDia_mm * 1e-3 + spacingMargin;

% Static default-pose transform of the proximal frame into the GLOBAL
% pelvis frame (cross-actuator spacing metric; envelope approximation).
ctx.frame1Body = proxBody;
Xdef = bodyTransforms(chain, 0, '');
ctx.TtoPelvis = Xdef.(safeField(proxBody));

%% Finalized routes (spacing + set footprint), pelvis-frame segments
ctx.finalizedNames = {};
ctx.finalizedSegs = {};
for r = 1:numel(finalizedRoutes)
    ctx.finalizedNames{end + 1} = finalizedRoutes(r).name; %#ok<AGROW>
    ctx.finalizedSegs{end + 1} = finalizedRoutes(r).Spelvis; %#ok<AGROW>
end
end

%% ---------------------------------------------------------------------
function v = getOrDefault(s, name, def)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = def;
end
end

function s = stripNote(s)
s = regexprep(strtrim(s), '\s*\(.*$', '');
end

function lbl = safeLabel(masterName)
lbl = strrep(masterName, ' ', '_');
lbl = strrep(lbl, ',', '');
lbl = strrep(lbl, '(', '');
lbl = strrep(lbl, ')', '');
lbl = strrep(lbl, '.', '');
end

function f = safeField(body)
% femur_r -> femur etc. (valid struct field names)
f = ['b_', regexprep(body, '_r$', '')];
end

function tf = rollingSpec(kin)
tf = (isfield(kin, 'transThetaX') && ~isempty(kin.transThetaX)) || ...
    (isfield(kin, 'transThetaY') && ~isempty(kin.transThetaY));
end

function p = rollPivot(kin, theta, p)
px = 0;
py = 0;
if ~isempty(kin.transThetaX)
    px = interp1(kin.transThetaX, kin.transX, theta, 'linear', 'extrap');
end
if ~isempty(kin.transThetaY)
    py = interp1(kin.transThetaY, kin.transY, theta, 'linear', 'extrap');
end
p = [px, py, p(3)];
end

function T = localHinge(jtype, axis, pivot, theta)
switch jtype
    case 'fixed'
        R = eye(3);
    otherwise
        a = axis(:) / norm(axis);
        Km = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
        R = eye(3) + sin(theta) * Km + (1 - cos(theta)) * (Km * Km);
end
T = [R, pivot(:); 0 0 0 1];
end

function p = segMapRow(points, cross, r, T1, T2, frame)
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
        p = (T2(1:3, 1:3) * p.' + T2(1:3, 4)).';
        f = 2;
    else
        p = (T1(1:3, 1:3) * p.' + T1(1:3, 4)).';
        f = 1;
    end
end
end

function key = jointEdgeKey(coordName)
switch coordName
    case 'hip_flexion_r'
        key = 'pelvis>femur_r';
    case 'knee_angle_r'
        key = 'femur_r>tibia_r';
    case 'ankle_angle_r'
        key = 'tibia_r>talus_r';
    case 'subtalar_angle_r'
        key = 'talus_r>calcn_r';
    case 'mtp_angle_r'
        key = 'calcn_r>toes_r';
    case 'lumbar_extension'
        key = 'pelvis>torso';
    otherwise
        key = '';
end
end

function X = bodyTransforms(chain, thDegPrimary, primKey)
% Per-body 4x4 transforms in the pelvis-anchored TRUE chain at one pose
% (primary joint at thDegPrimary, every other joint at its default 0).
edges = {'pelvis>femur_r', 'femur_r>tibia_r', 'tibia_r>talus_r', ...
    'talus_r>calcn_r', 'calcn_r>toes_r'};
parents = {'pelvis', 'femur_r', 'tibia_r', 'talus_r', 'calcn_r'};
X = struct();
X.b_pelvis = eye(4);
X.b_torso = chain.edges('pelvis>torso').Tdefault;
for e = 1:numel(edges)
    th = 0;
    if strcmp(edges{e}, primKey)
        th = thDegPrimary * pi / 180;
    end
    E = routingChainLib('edgeTrans', chain.edges(edges{e}), th, true);
    child = edges{e}(strfind(edges{e}, '>') + 1:end);
    X.(safeField(child)) = X.(safeField(parents{e})) * E;
end
end
