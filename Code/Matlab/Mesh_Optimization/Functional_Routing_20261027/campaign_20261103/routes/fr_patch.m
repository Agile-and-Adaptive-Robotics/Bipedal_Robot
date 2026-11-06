function P = fr_patch(C, boneName, seedBodyFrame, seedLoc, patchRadius_mm)
%FR_PATCH Bounded anatomical surface patch + surface projection.
%
% master directive item 3: attachment locations move within bounded,
% anatomically appropriate surface patches on their assigned bones, with
% surface projection so a Cartesian search step cannot carry a mount far
% beyond the OpenSim bone surface.
%
% Patch definition (documented per attachment):
%   - bone surface cloud: Bone_Mesh_Plots/Open_Sim_Bone_Geometry/
%     <bone>_Mesh_Points.txt (points in the bone's own body frame; the
%     layering CSV of geometry_20261103 measured body-frame attachment
%     distances against exactly these clouds)
%   - patch = the ball of radius patchRadius_mm (millimeters) around the
%     OpenSim seed, parameterized by (u, v) in [-1, 1]^2 along the patch
%     frame (e1, e2)
%   - candidate (u, v, s) maps to
%       c* = nearest cloud point to (seed + Rpatch*(u*e1 + v*e2))
%       p  = c* + (s/1000) * n_hat(c*)
%     with s the MOUNT STANDOFF in millimeters (bracket plate + mount,
%     bounded by C.standoff_mm). Every candidate is therefore on the bone
%     surface band (nearest-cloud granularity, documented per cloud:
%     femur 456 pts, tibia 144, pelvis_R 124, sacrum 191, calcaneus 999,
%     talus 98, toes 957, spine 5433 - geometry_20261103 README) plus an
%     explicit hardware standoff.
%
% P.project(u, v, s_mm) -> 1x3 mount position in seedBodyFrame.
% P.n: outward surface normal at the patch center.

cloud = fr_loadCloud(C, boneName);
seed = seedLoc(:)';
assert(all(isfinite(seed)),'fr_patch:InvalidSeed','Non-finite attachment seed on %s (%s)',boneName,seedBodyFrame);
r_local = min(0.08, 3 * patchRadius_mm / 1000);
d2 = sum((cloud - seed).^2, 2);
near = cloud(d2 <= r_local^2, :);
if size(near, 1) < 8
    d2s = sqrt(d2);
    [~, ord] = sort(d2s, 'ascend');
    near = cloud(ord(1:min(15, numel(ord))), :);
end
[~, ~, S] = svd(near - mean(near, 1), 0);
n = S(:, 3)';                        % smallest-variance direction = normal
if dot(n, seed - mean(near, 1)) < 0  % orient outward
    n = -n;
end
e1 = S(:, 1)';
if abs(dot(e1, n)) > 0.9
    e1 = [1 0 0];
    if abs(dot(e1, n)) > 0.9; e1 = [0 1 0]; end
    e1 = e1 - dot(e1, n) * n;
    e1 = e1 / norm(e1);
end
e2 = cross(n, e1);
e2 = e2 / max(norm(e2), eps);

e1c = e1; e2c = e2; nc = n; Rpatch = patchRadius_mm / 1000;
P = struct();
P.bone = boneName;
P.body = seedBodyFrame;
P.seed = seed;
P.Rpatch_m = Rpatch;
P.n = nc;
P.e1 = e1c;
P.e2 = e2c;
P.cloud = cloud;
P.project = @(u, v, s_mm) fr_project(cloud, seed, e1c, e2c, Rpatch, ...
    nc, s_mm, u, v);
end

function p = fr_project(cloud, seed, e1, e2, Rpatch, nAtSeed, s_mm, u, v)
% tangential-in-patch + normal-standoff placement: anchor at the nearest
% cloud point to the in-plane candidate, then slide WITHIN the patch along
% the local tangent plane (e1, e2) and stand off along the normal. The
% tangential slide keeps the mount on the bone-surface band (patch radius
% bounds it; the bone-exclusion constraint still guards true clearance)
% while staying a CONTINUOUS variable (nearest-cloud-only placement would
% freeze u,v on sparse clouds).
uc = max(-1, min(1, u));
vc = max(-1, min(1, v));
r=hypot(uc,vc);if r>1;uc=uc/r;vc=vc/r;end
c = seed + Rpatch * (uc * e1 + vc * e2);
d2 = sum((cloud - c).^2, 2);
[~, k] = min(d2);
% c already contains the tangential displacement. Adding the same u/v
% displacement again after nearest-point selection doubles the search.
delta=c-cloud(k,:);
tangent=dot(delta,e1)*e1+dot(delta,e2)*e2;
p = cloud(k, :) + tangent + (s_mm / 1000) * nAtSeed;
p = p(:)';
end

function cloud = fr_loadCloud(C, boneName)
f = fullfile(C.boneCloudDir, [boneName '_Mesh_Points.txt']);
if ~exist(f, 'file')
    error('fr_patch:NoCloud', 'bone cloud not found: %s', f);
end
cloud = load(f);
if isstruct(cloud)
    fn = fieldnames(cloud);
    cloud = cloud.(fn{1});
end
cloud = reshape(cloud, [], 3);
end
