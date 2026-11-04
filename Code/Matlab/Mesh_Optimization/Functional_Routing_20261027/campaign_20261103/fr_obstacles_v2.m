function [d, bid, perBody] = fr_obstacles_v2(C, ptsWorld, kin, q)
%FR_OBSTACLES_V2 Registered CAD obstacle distance model (task 1).
%
% Per-body SDF grids built from the ACTUAL robot CAD meshes
% (Code/MuJoCo_SNS/robot_torso split meshes, built by build_ah001.py and
% registered into the OpenSim gait2392_robotbody frames), plus cloud-only
% bodies (Bone_Mesh_Plots) where the robot CAD has no part. Built by
% tmp/functional_routing_20261103/packaging_v2_build/build_sdf.py; the
% grids live in packaging_v2/sdf/*.mat with the frame chain recorded in
% each file. Registration quality (mm numbers, uncertainty bounds):
% packaging_v2/registration_report.json.
%
% CALL FORMS:
%   model = fr_obstacles_v2(C)            % manifest: bodies, kinds, chains
%   [d, bid, perBody] = fr_obstacles_v2(C, ptsWorld, kin, q)
%       ptsWorld  Nx3 ground-frame query points at pose q (meters)
%       d         Nx1 min signed distance over all bodies (m; negative
%                 = inside the CAD surface; clouds return >= 0)
%       bid       Nx1 nearest-body id (index into model.bodies)
%       perBody(i).body/.kind/.dist  per-body distances (N x 1)
%
% TRANSFORM CHAIN (recorded per .mat, composed here per pose):
%   robot body frames are ISB-oriented (x anterior, y superior, z right)
%   at the pelvis datum (build_ah001.py cad_to_isb pelvis-part ruling;
%   AGENTS.md 2026-10-27). OpenSim leg frames at neutral are ground-aligned
%   and the OpenSim femur_r origin IS the hip center (validated in the
%   build against geometry_20261103 bodyTransformsInGround.neutral), so:
%     thigh  : M_ground<-robotThigh(q) = fr_bodyT(kin,q,'femur_r')
%              (both origins are the hip; rotations are the same hip delta)
%     back_r : thigh * trans(b1) * Rot(axBack, thB(qk))
%     shank_r: thigh * trans(b1)*Rot(axBack,thB) * trans(dbl)*Rot(axShank,thS)
%     front_r: thigh * trans(f1) * Rot(axFront, thF(qk))
%     foot_r : shank * trans(footOff)   (welded foot)
%     pelvis : fr_bodyT(kin,q,'pelvis') * trans(-treg)
%              with treg = hipOS_pelvis - hipR_pelvis (hip-aligned
%              pelvis-datum registration)
%   and (thB, thS, thF) interpolate the recorded four-bar knee table over
%   q.knee_angle_r.
%
% CLOUD BODIES (Spine/Pelvis_R/Sacrum/Talus/Calcaneus/Toes): distance =
% nearest cloud point in the owning OpenSim body frame (fr_bodyT per pose).
% A sparse cloud OVERestimates the true surface distance (conservative for
% clearance); the per-cloud NN-spacing uncertainty bound is in
% registration_report.json / model.clouds. Bodies with NO CAD mesh and NO
% cloud carry NO obstacle entry and must be stated as uncovered (none in
% this build: every obstacle body is mesh- or cloud-backed).
%
% USAGE RULE (registration report): attachment searches may use this model
% for clearance, but a large attachment excursion justified by uncertain
% registration is FORBIDDEN - the campaign patch bounds stay anatomically
% tight (surface patch + mount standoff in mm, fr_patch/fr_setup).
%
% SDF convention: meters; negative inside; clamped to +/-30 mm; queries
% outside a grid return +clamp (conservative). Trilinear interpolation.

% ---- mode 1: manifest -----------------------------------------------------
if nargin == 1
    d = fr_manifest(C);
    return
end

% ---- mode 2: distance query ------------------------------------------------
persistent cache
if isempty(cache) || ~isfield(cache, 'dir') || ...
        ~strcmp(cache.dir, fullfile(C.camp, 'packaging_v2', 'sdf'))
    cache = fr_loadAll(C);
end
pts = double(ptsWorld);
n = size(pts, 1);
bodies = cache.bodies;
nb = numel(bodies);
perBody = struct('body', {}, 'kind', {}, 'dist', {}, 'conservative', {});
D = zeros(n, nb);
for b = 1:nb
    B = bodies(b);
    switch B.kind
        case 'cad'
            pp = fr_toPieceFrame(B, pts, kin, q, cache);
            D(:, b) = fr_sdfSample(B.grid, pp, B.clamp_m);
        case 'cloud'
            T = fr_bodyT(kin, q, B.osimCounterpart);
            pp = (T(1:3, 1:3) * pts' + T(1:3, 4))';
            D(:, b) = min(frustumDist(pp, B.pts));
    end
    perBody(b).body = B.body;
    perBody(b).kind = B.kind;
    perBody(b).dist = D(:, b);
    perBody(b).conservative = strcmp(B.kind, 'cloud');
end
[d, bid] = min(D, [], 2);
d = d(:);
bid = bid(:);
end

% ==========================================================================
function model = fr_manifest(C)
persistent mcache
if isempty(mcache) || ~strcmp(mcache.dir, fullfile(C.camp, ...
        'packaging_v2', 'sdf'))
    mcache = fr_loadAll(C);
end
model = struct();
model.bodies = {mcache.bodies.body};
model.kinds = {mcache.bodies.kind};
model.osimCounterpart = {mcache.bodies.osimCounterpart};
model.sdfDir = mcache.dir;
model.kneeTable = mcache.kneeTable;
model.stats = {mcache.bodies.stats};
model.clouds = mcache.clouds;
model.chainNote = ['robot frames ISB at the pelvis datum; OpenSim leg ' ...
    'frames ground-aligned at neutral; femur_r origin = hip center; ' ...
    'full chains + mm registration numbers in ' ...
    'packaging_v2/registration_report.json'];
model.usageRule = ['clearance use only: attachment excursions justified ' ...
    'by registration uncertainty are FORBIDDEN; patch bounds stay ' ...
    'anatomically tight (fr_patch/fr_setup conventions)'];
end

% ==========================================================================
function cache = fr_loadAll(C)
sdfDir = fullfile(C.camp, 'packaging_v2', 'sdf');
cache.dir = sdfDir;
d = dir(fullfile(sdfDir, 'sdf_*.mat'));
clouds = dir(fullfile(sdfDir, 'cloud_*.mat'));
bodies = struct('body', {}, 'kind', {}, 'osimCounterpart', {}, ...
    'grid', {}, 'clamp_m', {}, 'pts', {}, 'stats', {}, 'chain', {});
for i = 1:numel(d)
    S = fr_loadsafe(fullfile(sdfDir, d(i).name));
    b = struct();
    b.body = S.body;
    b.kind = S.kind;
    b.osimCounterpart = S.osimCounterpart;
    b.grid.origin = S.origin(:)';
    b.grid.spacing = double(S.spacing);
    b.grid.dims = double(S.dims(:)');
    b.grid.values = double(S.sdf);
    b.clamp_m = double(S.clamp_m);
    b.chain = jsondecode(S.chain);
    b.stats = jsondecode(S.stats);
    b.pts = [];
    if isempty(bodies)
        bodies = b;
    else
        bodies(end + 1) = b; %#ok<AGROW>
    end
end
for i = 1:numel(clouds)
    S = fr_loadsafe(fullfile(sdfDir, clouds(i).name));
    b = struct();
    b.body = S.body;
    b.kind = S.kind;
    b.osimCounterpart = S.osimCounterpart;
    b.pts = double(S.pts);
    b.grid = [];
    b.clamp_m = 0;
    b.stats = jsondecode(S.uncertainty);
    b.chain = [];
    bodies(end + 1) = b; %#ok<AGROW> %#ok<AGROW>
end
cache.bodies = bodies;
kt = fr_loadsafe(fullfile(sdfDir, 'knee_table.mat'));
cache.kneeTable = struct('knee_angle_rad', kt.knee_angle_rad(:), ...
    'theta_shank', kt.theta_shank(:), 'theta_back', kt.theta_back(:), ...
    'theta_front', kt.theta_front(:));
cc = jsondecode(fileread(fullfile(C.camp, 'packaging_v2', ...
    'registration_report.json')));
cache.clouds = cc.cloud_uncertainty;
end

function S = fr_loadsafe(f)
% load() returns one field per saved variable here (multi-variable mats);
% return it directly. Read-only use; never saved afterwards.
S = load(f);
end

% ==========================================================================
function pp = fr_toPieceFrame(B, pts, kin, q, cache)
% ground-frame points -> the piece (robot body) frame at pose q
thS = interp1(cache.kneeTable.knee_angle_rad, cache.kneeTable.theta_shank, ...
    q.knee_angle_r, 'linear');
thBk = interp1(cache.kneeTable.knee_angle_rad, cache.kneeTable.theta_back, ...
    q.knee_angle_r, 'linear');
thF = interp1(cache.kneeTable.knee_angle_rad, cache.kneeTable.theta_front, ...
    q.knee_angle_r, 'linear');
ch = B.chain;
T = fr_bodyT(kin, q, 'femur_r');          % ground <- robot thigh (hip)
switch B.body
    case 'thigh_r'
        % identity piece (meshes are in the thigh frame)
    case 'back_r'
        T = T * fr_trans(ch.b1) * fr_rot(ch.axBack, thBk);
        T = inv(T);
    case 'front_r'
        T = T * fr_trans(ch.f1) * fr_rot(ch.axFront, thF);
        T = inv(T);
    case 'shank_r'
        T = T * fr_trans(ch.b1) * fr_rot(ch.axBack, thBk) ...
            * fr_trans(ch.dbl) * fr_rot(ch.axShank, thS);
        T = inv(T);
    case 'foot_r'
        T = T * fr_trans(ch.b1) * fr_rot(ch.axBack, thBk) ...
            * fr_trans(ch.dbl) * fr_rot(ch.axShank, thS) ...
            * fr_trans(ch.footOff);
        T = inv(T);
    case 'pelvis'
        Tp = fr_bodyT(kin, q, 'pelvis') * fr_trans(-ch.tregP);
        T = inv(Tp);
    otherwise
        error('fr_obstacles_v2:Body', 'unknown CAD body %s', B.body);
end
pp = (T(1:3, 1:3) * pts' + T(1:3, 4))';
end

function d = fr_sdfSample(g, pp, clamp)
% trilinear SDF sample; outside-grid queries return +clamp (conservative)
dims = g.dims;
n = size(pp, 1);
d = repmat(clamp, n, 1);
idx = (pp - repmat(g.origin, n, 1)) ./ repmat(g.spacing, n, 1);
in = all(idx >= 0 & idx <= repmat(dims - 1, n, 1), 2);
if ~any(in)
    return
end
i0 = floor(idx(in, :));
i1 = min(i0 + 1, dims - 1);
f = idx(in, :) - i0;
lin = @(a, b, t) a .* (1 - t) + b .* t;
c000 = fr_at(g.values, i0, dims);
c100 = fr_at(g.values, [i0(:, 1) + 1, i0(:, 2), i0(:, 3)], dims);
c010 = fr_at(g.values, [i0(:, 1), i0(:, 2) + 1, i0(:, 3)], dims);
c001 = fr_at(g.values, [i0(:, 1), i0(:, 2), i0(:, 3) + 1], dims);
c110 = fr_at(g.values, [i0(:, 1) + 1, i0(:, 2) + 1, i0(:, 3)], dims);
c101 = fr_at(g.values, [i0(:, 1) + 1, i0(:, 2), i0(:, 3) + 1], dims);
c011 = fr_at(g.values, [i0(:, 1), i0(:, 2) + 1, i0(:, 3) + 1], dims);
c111 = fr_at(g.values, i1, dims);
v = lin(lin(lin(c000, c100, f(:, 1)), lin(c010, c110, f(:, 1)), f(:, 2)), ...
    lin(lin(c001, c101, f(:, 1)), lin(c011, c111, f(:, 1)), f(:, 2)), ...
    f(:, 3));
d(in) = v;
end

function v = fr_at(V, idx, dims)
v = V(sub2ind(dims, idx(:, 1) + 1, idx(:, 2) + 1, idx(:, 3) + 1));
end

function dd = frustumDist(pp, P)
% brute-force nearest cloud point distance, chunked
n = size(pp, 1);
dd = zeros(n, 1);
chunk = 2000;
for s = 1:chunk:n
    e = min(n, s + chunk - 1);
    Q = pp(s:e, :);
    dd(s:e) = min(sqrt(max(0, ...
        sum(Q .^ 2, 2) + sum(P .^ 2, 2)' - 2 * (Q * P'))), [], 2);
end
end

function T = fr_trans(t)
T = eye(4);
T(1:3, 4) = t(:);
end

function R = fr_rot(axis, theta)
% Rodrigues rotation about an axis through the origin
a = axis(:) / norm(axis);
K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(4);
R(1:3, 1:3) = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end
