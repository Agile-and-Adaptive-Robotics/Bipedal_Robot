function obs = fr_obstacles(C, caseDef, loc, rowBodies, kin, qPoseList, ...
    envelopes)
%FR_OBSTACLES Hard obstacle margins for one candidate route.
%
% Two obstacle families (both HARD in the constraint set):
%   1. BONE: the route's straight sleeve segments (row r -> row r+1 for
%      BPA-side rows) against the owning bones' point clouds, clearance =
%      inflated sleeve radius + C.clearanceBone_m. This is the
%      buildGeoExclusion.m / nonlconExclusion.m lineage carried over to
%      bone point clouds (same "distance >= physical radius + rigid
%      clearance" form; the SolidWorks analytic profiles of the knee
%      context are knee-specific and the campaign uses the bone clouds
%      per the geometry layering table's convention).
%   2. RETAINED LAYER ENVELOPES (Ben's layering ruling): capsules /
%      spheres from fr_envelope of already-placed groups, evaluated at the
%      CURRENT case's pose list (physically simultaneous configurations),
%      clearance = C.clearanceEnvelope_m. A CAD-level (SolidWorks)
%      hardware-obstacle registration layer is a documented LATER campaign
%      item: obs.cadHook documents the exact call site - it is NOT faked.
%
% obs.c(c): c <= 0 feasible (violation form, ready for nonlcon).
% obs.active(boneIdx/envelopeIdx): TRUE where the margin <= tolerance
% (constraint active - used by the layering proof in the gate).

nG = numel(qPoseList);
sleeveR = (fr_diameter(C, caseDef.diaStr, 0) / 2) / 1000;

% --- candidate route world segments (BPA-side rows only: 1..exitRow) ----
exitRow = caseDef.cross(1);
segs = cell(1, nG);
for i = 1:nG
    q = qPoseList{i};
    w = zeros(exitRow, 3);
    for r = 1:exitRow
        T = fr_bodyT(kin, q, rowBodies{r});
        w(r, :) = (T(1:3, 1:3) * loc(r, :)' + T(1:3, 4))';
    end
    segs{i} = w;
end

% --- bone obstacle --------------------------------------------------------
obs.boneClouds = {};
obs.boneMargin = inf;
obs.boneActive = false;
for b = 1:numel(caseDef.obstacleBones)
    cloud = fr_loadCloudCache(C, caseDef.obstacleBones{b});
    m = inf;
    for i = 1:nG
        w = segs{i};
        for r = 1:size(w, 1) - 1
            d = segPointDist(cloud, w(r, :), w(r + 1, :));
            m = min(m, d - sleeveR - C.clearanceBone_m);
        end
    end
    obs.boneClouds{end+1} = struct('bone', caseDef.obstacleBones{b}, ...
        'margin_m', m);
    obs.boneMargin = min(obs.boneMargin, m);
end
obs.boneActive = obs.boneMargin <= 1e-4;

% --- retained layer envelopes --------------------------------------------
obs.envMargin = inf;
obs.envActive = false;
obs.envWhich = '';
obs.envDeflection_m = NaN;   % set by the layering driver when known
for e = 1:numel(envelopes)
    envE = envelopes{e};
    if ~isstruct(envE); continue; end
    m = inf;
    for i = 1:nG
        w = segs{i};
        for r = 1:size(w, 1) - 1
            for k = 1:numel(envE.capsules{i})
                cap = envE.capsules{i}(k);
                d = segSegDist(w(r, :), w(r + 1, :), cap.p1, cap.p2) ...
                    - cap.r - sleeveR;
                m = min(m, d - C.clearanceEnvelope_m);
            end
            sph = envE.spheres{i};   % single struct; spheres = rows of .p
            for k = 1:size(sph.p, 1)
                d = segPointDist(sph.p(k, :), w(r, :), w(r + 1, :)) ...
                    - sph.r(k) - sleeveR;
                m = min(m, d - C.clearanceEnvelope_m);
            end
        end
    end
    assert(isscalar(m), 'fr_obstacles: non-scalar envelope margin');
    if m < obs.envMargin
        obs.envMargin = m;
        obs.envWhich = envE.caseName;
    end
end
obs.envActive = obs.envMargin <= 1e-4;

% --- documented CAD hook (NOT faked) --------------------------------------
obs.cadHook = ['CAD-level obstacle registration (SolidWorks) is a LATER ' ...
    'campaign item; the registration call site is fr_obstacles:cadHook - ' ...
    'when the CAD obstacles land, append their primitives to the ' ...
    'envelope loop below with provenance frames. Nothing is simulated ' ...
    'for them in stage 1.'];
% constraint values (violation form): [bone env]
obs.c = [-obs.boneMargin, -obs.envMargin];
end

function d = segPointDist(pts, a, b)
% min distance from segment a-b to a point set (vectorized)
ab = b - a;
L2 = dot(ab, ab);
if L2 < 1e-18
    d = min(vecnorm(pts - a, 2, 2));
    return
end
t = max(0, min(1, ((pts - a) * ab') / L2));
proj = a + t * ab;
d = min(vecnorm(pts - proj, 2, 2));
end

function d = segSegDist(p1, p2, q1, q2)
% closest distance between two segments (standard clamped algorithm)
u = p2 - p1; v = q2 - q1; w0 = p1 - q1;
a = dot(u, u); b = dot(u, v); c = dot(v, v);
d0 = dot(u, w0); e = dot(v, w0);
D = a * c - b * b;
if D < 1e-24
    % near-parallel: clamp on q's segment only
    sc = 0;
    tc = max(0, min(1, e / max(c, eps)));
else
    % standard two-pass clamp
    tc = max(0, min(1, (a * e - b * d0) / D));
    sc = max(0, min(1, (b * tc - d0) / a));
    tc = max(0, min(1, (b * sc - e) / c));
end
d = norm(p1 + sc * u - (q1 + tc * v));
end

function cloud = fr_loadCloudCache(C, boneName)
persistent cache cacheKeys
if isempty(cacheKeys)
    cache = {};
    cacheKeys = {};
end
k = [boneName, '@', C.boneCloudDir];
idx = find(strcmp(cacheKeys, k), 1);
if isempty(idx)
    f = fullfile(C.boneCloudDir, [boneName '_Mesh_Points.txt']);
    cl = load(f);
    if isstruct(cl)
        fnl = fieldnames(cl);
        cl = cl.(fnl{1});
    end
    idx = numel(cacheKeys) + 1;
    cacheKeys{idx} = k;
    cache{idx} = reshape(cl, [], 3);
end
cloud = cache{idx};
end
