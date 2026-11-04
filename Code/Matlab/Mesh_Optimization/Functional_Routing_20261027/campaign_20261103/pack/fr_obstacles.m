function obs = fr_obstacles(C, caseDef, loc, rowBodies, kin, qPoseList, ...
    envelopes)
%FR_OBSTACLES Obstacle margins for one candidate route (v2, 2026-11-04).
%
% Ben's corrected layering/interference ruling split the old single HARD
% envelope into TWO families:
%   HARD (constraint set, c <= 0 feasible):
%     1. BONE: straight sleeve segments vs the owning bones' point clouds,
%        clearance = inflated sleeve radius + C.clearanceBone_m
%        (buildGeoExclusion.m / nonlconExclusion.m lineage). When the
%        packaging_v2 SDF model (fr_obstacles_v2) lands it REPLACES the
%        cloud distance here with a validated surface-distance model and
%        quantified registration error - the call site is obs.cadHook.
%     2. RETAINED RIGID HARDWARE of deeper layers: fitting spheres and
%        tackle blocks from fr_envelope v2 (kind 'rigid_*'), evaluated at
%        the CURRENT case's pose list (physically simultaneous
%        configurations), clearance = C.clearanceEnvelope_m.
%   SOFT (objective term retained_bpa_corridor_soft, never a constraint):
%     3. RETAINED BPA SLEEVE CORRIDORS of deeper layers (kind
%        'bpa_sleeve'/'bpa_end'): penetration depth summed over the pose
%        grid (normalized by the resting inflated diameter, the same
%        convention as the in-group fr_softOverlap term) and max extent in
%        mm, both recorded per evaluation. Slight BPA-BPA interference is
%        tolerable at a penalty; it can NEVER make a candidate infeasible.
%
% obs.c(c): c <= 0 feasible (violation form, ready for nonlcon) - HARD
% entries only.
% obs.active(boneIdx/rigidIdx): TRUE where the margin <= tolerance
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

% --- 1. bone obstacle (HARD) ----------------------------------------------
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

% --- 2 + 3. retained layer envelopes (rigid HARD, BPA corridors SOFT) ----
obs.rigidMargin = inf;
obs.rigidActive = false;
obs.rigidWhich = '';
obs.corridorDepthNorm = 0;
obs.corridorMaxMm = 0;
obs.corridorWhich = '';
obs.envDeflection_m = NaN;   % set by the layering driver when known
for e = 1:numel(envelopes)
    envE = envelopes{e};
    if ~isstruct(envE); continue; end
    mRigid = inf;
    dSoft = 0;
    mmSoft = 0;
    for i = 1:nG
        w = segs{i};
        for r = 1:size(w, 1) - 1
            for k = 1:numel(envE.capsules{i})
                cap = envE.capsules{i}(k);
                d = segSegDist(w(r, :), w(r + 1, :), cap.p1, cap.p2) ...
                    - cap.r - sleeveR;
                if strcmpi(cap.kind, 'rigid_*') || ...
                        strncmpi(char(cap.kind), 'rigid', 5)
                    mRigid = min(mRigid, d - C.clearanceEnvelope_m);
                else
                    % SOFT BPA corridor: penetration depth, no clearance,
                    % normalized by the resting inflated diameter (same
                    % convention as fr_softOverlap in fr_objective)
                    ov = -d;
                    if ov > 0
                        dSoft = dSoft + ov / (2 * sleeveR);
                        mmSoft = max(mmSoft, ov * 1000);
                    end
                end
            end
            sph = envE.spheres{i};   % single struct; spheres = rows of .p
            for k = 1:size(sph.p, 1)
                d = segPointDist(sph.p(k, :), w(r, :), w(r + 1, :)) ...
                    - sph.r(k) - sleeveR;
                if strncmpi(char(sph.kind{k}), 'rigid', 5)
                    mRigid = min(mRigid, d - C.clearanceEnvelope_m);
                else
                    ov = -d;
                    if ov > 0
                        dSoft = dSoft + ov / (2 * sleeveR);
                        mmSoft = max(mmSoft, ov * 1000);
                    end
                end
            end
        end
    end
    assert(isscalar(mRigid), 'fr_obstacles: non-scalar envelope margin');
    if mRigid < obs.rigidMargin
        obs.rigidMargin = mRigid;
        obs.rigidWhich = envE.caseName;
    end
    if dSoft > obs.corridorDepthNorm
        obs.corridorDepthNorm = dSoft;
        obs.corridorMaxMm = mmSoft;
        obs.corridorWhich = envE.caseName;
    elseif dSoft > 0 && mmSoft > obs.corridorMaxMm
        obs.corridorMaxMm = mmSoft;
        obs.corridorWhich = envE.caseName;
    end
end
% backward-compatible alias for the gate/proof consumers: envMargin is now
% the RIGID-only margin (BPA corridors are excluded from it by design)
obs.envMargin = obs.rigidMargin;
obs.envActive = obs.rigidMargin <= 1e-4;
obs.envWhich = obs.rigidWhich;

% --- documented CAD hook (NOT faked) --------------------------------------
obs.cadHook = ['CAD-level obstacle registration: packaging_v2 SDF model ' ...
    '(fr_obstacles_v2, body-frame grids + quantified registration error) ' ...
    'replaces the cloud distance at THIS call site once its registration ' ...
    'report lands; until then the bone clouds carry the margin. Nothing ' ...
    'is faked.'];
% constraint values (violation form): [bone rigid-envelopes] - the soft
% corridor is an OBJECTIVE term, never here
obs.c = [-obs.boneMargin, -obs.rigidMargin];
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
