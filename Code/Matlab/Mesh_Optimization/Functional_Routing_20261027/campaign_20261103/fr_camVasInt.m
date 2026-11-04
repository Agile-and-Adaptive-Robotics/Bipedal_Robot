function varargout = fr_camVasInt(C, cs, xCam, kin, qList, opts, varargin)
%FR_CAMVASINT vas_int cam-assisted transmission candidate (task 3).
% (varargin carries the bounds-mode payload cs: fr_camVasInt('bounds', C, cs);)
%FR_CAMVASINT vas_int cam-assisted transmission candidate (task 3).
%
% Handoff item 3: vas_int replaces the combined vasti torque with several
% parallel 20 mm BPAs and MAY use a cam-assisted transmission with the
% terminal attachment at the tibial patellar-tendon attachment or inferior
% to it within bounded surface patches; cam profile, pivot, cable contact,
% bracket and insertion are OPTIMIZATION VARIABLES. This file parameterizes
% and evaluates that candidate in the Opt_run_ExtCam lineage style (flat
% design vector + lb/ub + closure vector), plugged into the fr_cases27 /
% fr_evalRoute / fr_mechBuild case structure.
%
% TOPOLOGY (femur frame, documented in packaging_v2/README_camVasInt.md):
%   - nBPA parallel 20 mm sleeves run straight from the femur origin
%     bracket (the case's movable-1 patch, seed of record) to the CAM
%     bracket; both are femur-fixed, so the sleeve length is pose constant.
%   - a cable leaves the cam bracket exit ON the cam profile, wraps the
%     taut path around the cam polygon, exits at a tangent vertex and runs
%     straight (3D) to the tibial terminal attachment.
%   - the tibial attachment sits on a 15 mm surface patch centered at the
%     seed patellar-tendon insertion (the case's moving-point row,
%     fixed-tracked), with an inferior-only guard: the attachment may not
%     move proximal of the seed + 2 mm in the tibia frame.
%   - the cam profile is POLAR: r(th) = r0 + a1*cos(th) + b1*sin(th) +
%     a2*cos(2 th) + b2*sin(2 th), sampled at 2 deg into the contact
%     polygon the cable wraps. The cam plane is the femur sagittal plane
%     z0 = the seed insertion's femur-frame z (constant: the OpenSim knee
%     axis is the femur-frame z, so the insertion's z does not drift --
%     checked and reported).
%
% DESIGN VECTOR (14x1, meters / radians / millimeters as noted):
%   xCam(1:2)  pivot (px, py) in the FEMUR frame at the cam plane z0 [m]
%   xCam(3)    r0   base radius [m]
%   xCam(4:7)  a1 b1 a2 b2  harmonics [m]
%   xCam(8)    thB  bracket exit angle on the cam [rad]
%   xCam(9)    sB   bracket mount standoff [mm] (C.standoff_mm)
%   xCam(10:11) u v  tibial attachment patch coords [-1, 1]
%   xCam(12)   s_tib attachment mount standoff [mm]
%   xCam(13)   L0   BPA rest length (bookkeeping) [m]
%   xCam(14)   tendon  cable bookkeeping length [m]
%
% CALL FORMS:
%   [lb, ub] = fr_camVasInt('bounds', C, cs)
%   out = fr_camVasInt(C, cs, xCam, kin, qList)
%   out = fr_camVasInt(C, cs, xCam, kin, qList, opts)
% opts.useObstacles (default true): the cable path samples are cleared
% against the registered CAD obstacle model (fr_obstacles_v2). This is
% clearance-only use: per the registration RULE, patch bounds stay
% anatomically tight regardless of what the obstacle model allows.
%
% OUTPUT (plug contract for a driver / fr_objective):
%   out.spanCam(i)     knee-crossing cable span per pose [m]
%   out.dSpan(i)       dSpan/dtheta on the case grid (central differences,
%                      the fr_objective chord convention)
%   out.tau_perT(i)    -dSpan/dtheta: per-unit-tension knee torque
%                      (tau = -T*dSpan/dtheta, fr_objective line 86)
%   out.momentArm(i)   |tau_perT| [m] (report value)
%   out.wrapAngle(i)   swept contact arc [rad] (feed fr_hardwareLoads)
%   out.pivotReaction_N_perT(i)  per-unit-tension pivot reaction magnitude
%   out.c, out.clabels closure vector c <= 0 feasible
%   out.geometry       pivot/polygon/bracket/exit/attachment per pose
%   out.L0, out.tendon, out.mlTotal(i)  bookkeeping (sleeve + span)
%   out.travelFeasible  dSpan <= kmaxF * L0 (G = 1 cam; tackle variants are
%                      the later driver's discrete choice)
%   out.zDriftMax_m    insertion cam-plane drift (plane-assumption check)
%
% The vas_int case of record: fr_cases27 group 11 (pbKneeMono; femur sleeve
% rows + moving-point tibial insertion; knee rides T1). The femur via rows
% (rows 2..cross-1) are NOT tension-path points in the cam topology (the
% straight sleeve + cam replaces them); their seed positions are reported
% in out.femurViaNote for the driver's topology-substitution record.

if ischar(C) && strcmp(C, 'bounds')
    % bounds-mode call fr_camVasInt('bounds', C, cs): the payload cs is
    % bound to the 3rd declared parameter (the xCam slot) - documented
    [lb, ub] = fr_bounds(xCam);
    varargout = {lb, ub};
    return
end
if nargin < 6
    opts = struct();
end
if ~isfield(opts, 'useObstacles')
    opts.useObstacles = true;
end
if ~isfield(opts, 'uvOrigin')
    opts.uvOrigin = [0 0 5];   % seed patch center, 5 mm standoff (family seed)
end

TH_POLY = 180;            % polygon samples (2 deg)
FIT = cs.fittingLength;

% ---- unpack the design ----------------------------------------------------
px = xCam(1); py = xCam(2);
r0 = xCam(3); a1 = xCam(4); b1 = xCam(5); a2 = xCam(6); b2 = xCam(7);
thB = xCam(8); sB_mm = xCam(9);
uAtt = xCam(10); vAtt = xCam(11); sAtt_mm = xCam(12);
L0 = xCam(13); tendon = xCam(14);

% ---- case geometry constants ----------------------------------------------
q0 = kin.qDefault;
T10 = cs.passes{1}.T1(kin, q0);          % tibia -> femur at neutral
insRow = cs.movable(end).row;
seedIns = cs.seedLoc(insRow, :);         % tibia frame (moving-point seed)
A0_fem = fr_applyT(T10, seedIns);
z0 = A0_fem(3);

% attachment patch: TIGHT 15 mm patellar-tendon patch on the Tibia cloud
patchR = min(cs.movable(end).radius_mm, 15);
insPatch = fr_patch(C, 'Tibia', 'tibia_r', seedIns, patchR);
orgPatch = cs.movable(1).patch;
O3 = orgPatch.project(opts.uvOrigin(1), opts.uvOrigin(2), opts.uvOrigin(3));

nG = numel(qList);
span = zeros(nG, 1);
wrapAng = zeros(nG, 1);
pivotR = zeros(nG, 1);
zA = zeros(nG, 1);
geom = repmat(struct('pivot', [0 0 0], 'polygon', zeros(TH_POLY, 3), ...
    'bracket', [0 0 0], 'exit', [0 0 0], 'attach', [0 0 0], ...
    'path', zeros(0, 3), 'wrapped', false), 1, nG);
minPathClear = inf;
obst = [];
for i = 1:nG
    q = qList{i};
    T1 = cs.passes{1}.T1(kin, q);
    A3 = fr_applyT(T1, insPatch.project(uAtt, vAtt, sAtt_mm));
    zA(i) = A3(3);
    P = [px, py, z0];
    th = 2 * pi * (0:TH_POLY - 1)' / TH_POLY;
    rv = r0 + a1 * cos(th) + b1 * sin(th) + a2 * cos(2 * th) + ...
        b2 * sin(2 * th);
    poly2 = [P(1) + rv .* cos(th), P(2) + rv .* sin(th)];
    poly3 = [poly2, repmat(z0, TH_POLY, 1)];
    iB = round(mod(thB, 2 * pi) / (2 * pi) * TH_POLY) + 1;
    B3 = poly3(iB, :);
    B3 = B3 + sB_mm / 1000 * ...
        [cos(th(iB)), sin(th(iB)), 0];   % bracket standoff, radial out
    A2 = A3(1:2);
    % taut path around the polygon
    [len2d, ~, path3, wrapped, sweepAng] = fr_tautPath(poly3, iB, ...
        A3, z0);
    span(i) = len2d;
    wrapAng(i) = sweepAng;
    geom(i).pivot = P;
    geom(i).polygon = poly3;
    geom(i).bracket = B3;
    geom(i).attach = A3;
    geom(i).wrapped = wrapped;
    if ~isempty(path3)
        geom(i).exit = path3(max(size(path3, 1) - 1, 1), :);
        geom(i).path = path3;
        % sampled path clearance vs the registered obstacles
        if opts.useObstacles
            % fr_obstacles_v2 caches its own grids internally
            clr = fr_pathClearance(C, path3, kin, q);
            minPathClear = min(minPathClear, clr);
        end
    else
        geom(i).exit = A3;
    end
    % pivot reaction direction: cable pulls from the exit toward the
    % attachment; the bracket side reacts the straight-line pull
    uOut = A3 - geom(i).exit;
    nu = norm(uOut);
    if nu > 1e-9
        uOut = uOut / nu;
    else
        uOut = [1 0 0];
    end
    uIn = poly3(mod(iB, TH_POLY) + 1, :) - B3;
    nu = norm(uIn);
    if nu > 1e-9
        uIn = uIn / nu;
    else
        uIn = uOut;
    end
    pivotR(i) = norm(uIn + uOut);
end

% dSpan/dtheta on the case grid (fr_objective chord convention)
dTheta = cs.passes{1}.th1(2) - cs.passes{1}.th1(1);
dSpan = zeros(nG, 1);
if nG >= 3
    dSpan(2:nG - 1) = (span(3:nG) - span(1:nG - 2)) / (2 * dTheta);
end
dSpan(1) = (span(2) - span(1)) / dTheta;
dSpan(nG) = (span(nG) - span(nG - 1)) / dTheta;

% ---- closure vector -------------------------------------------------------
rvAll = r0 + abs(a1) + abs(b1) + abs(a2) + abs(b2);   % conservative max
rMin = r0 - (abs(a1) + abs(b1) + abs(a2) + abs(b2));
kmaxF = fr_kmaxF(C, cs.diaStr);
spanRange = max(span) - min(span);
seedY = seedIns(2);   % tibia frame: y superior; proximal = +y
zDrift = max(abs(zA - z0));

c = zeros(8, 1);
clabels = cell(8, 1);
c(1) = 0.006 - rMin;                 clabels{1} = 'cam_min_radius';
c(2) = rvAll - 0.045;                clabels{2} = 'cam_max_radius';
c(3) = max(wrapAng) - deg2rad(240);  clabels{3} = 'cam_wrap_window';
c(4) = 0;                            clabels{4} = 'attachment_proximal_guard';
% the guard is enforced on the PROJECTED attachment in the tibia frame:
attTib = insPatch.project(uAtt, vAtt, sAtt_mm);
c(4) = attTib(2) - (seedY + 0.002);
c(5) = 0;                            clabels{5} = 'cam_plane_drift';
c(5) = zDrift - 0.005;
c(6) = spanRange - kmaxF * L0;       clabels{6} = 'travel_band_G1';
c(7) = 0;                            clabels{7} = 'cable_path_clearance';
if opts.useObstacles && isfinite(minPathClear)
    c(7) = -(minPathClear - 0.003);  % 3 mm rigid clearance (fr_setup)
end
c(8) = L0 - spanRange / kmaxF - 0.03;   % rest length stays sensible
clabels{8} = 'rest_vs_travel';

out = struct();
out.xCam = xCam(:)';
out.lb = [];
out.ub = [];
out.spanCam = span;
out.dSpan = dSpan;
out.tau_perT = -dSpan;
out.momentArm = abs(dSpan);
out.wrapAngle = wrapAng;
out.pivotReaction_N_perT = pivotR;
out.c = c;
out.clabels = clabels;
out.feasible = all(c <= 1e-12);
out.zDriftMax_m = zDrift;
out.camPlaneZ_m = z0;
out.L0 = L0;
out.tendon = tendon;
out.mlTotal = norm(O3 - geom(1).bracket) + span;   % sleeve + cable
out.travelFeasible = c(6) <= 0;
out.sleeveLength_m = norm(O3 - geom(1).bracket);
out.nBPA = cs.seed.nBPA;
out.diaStr = cs.diaStr;
out.minPathClearance_m = minPathClear;
out.femurViaNote = sprintf( ...
    ['cam topology: %d femur via row(s) (rows 2..%d) are not tension-' ...
    'path points; the straight sleeve + cam bracket replaces them. Seed ' ...
    'positions retained in cs.seedLoc for the substitution record.'], ...
    max(cs.cross(1) - 2, 0), cs.cross(1) - 1);
out.geometry = geom;
out.seedAttachFem_neutral = A0_fem;
out.summary = sprintf( ...
    ['FR_CAMVASINT %s: span %.1f..%.1f mm (range %.1f), arm %.1f..%.1f ' ...
    'mm, wrap max %.1f deg, cmax %.2e, feasible %d'], ...
    cs.name, 1e3 * min(span), 1e3 * max(span), 1e3 * spanRange, ...
    1e3 * min(out.momentArm), 1e3 * max(out.momentArm), ...
    rad2deg(max(wrapAng)), max(c), out.feasible);
varargout{1} = out;
if nargout >= 2
    [~, ub] = fr_bounds(cs);
    varargout{2} = ub;
end
end

% ==========================================================================
function [lb, ub] = fr_bounds(cs)
q0seed = cs.seedLoc(cs.movable(end).row, :);
% pivot box around the distal-anterior femur relative to the neutral
% attachment (anatomically tight mechanical-bracket zone; see README)
lb = zeros(14, 1);
ub = zeros(14, 1);
lb(1) = q0seed(1) * 0 + -0.03;   % px [m]: box around the distal femur
ub(1) = 0.12;
lb(2) = -0.50;                   % py [m]: distal femur band (femur origin
ub(2) = -0.33;                   % = hip; knee ~ -0.40)
lb(3) = 0.008;  ub(3) = 0.030;   % r0
lb(4) = -0.008; ub(4) = 0.008;   % a1
lb(5) = -0.008; ub(5) = 0.008;   % b1
lb(6) = -0.006; ub(6) = 0.006;   % a2
lb(7) = -0.006; ub(7) = 0.006;   % b2
lb(8) = 0;      ub(8) = 2 * pi;  % thB
lb(9) = 1;      ub(9) = 12;      % sB mm
lb(10) = -1;    ub(10) = 1;      % u
lb(11) = -1;    ub(11) = 1;      % v
lb(12) = 1;     ub(12) = 12;     % s_tib mm
lb(13) = 0.4 * cs.L0Seed;
ub(13) = 2.5 * cs.L0Seed;
lb(14) = 2 * cs.fittingLength + 0.002;
ub(14) = 2 * cs.fittingLength + 0.05;
end

% ==========================================================================
function [lenTotal, iExit, path3, wrapped, sweepAng] = fr_tautPath(poly3, ...
    iB, A3, z0)
% Taut cable from the bracket exit vertex iB (ON the polygon) around the
% closed polygon to the attachment A3. Two candidate directions (CCW / CW
% boundary chains to every vertex) plus the direct no-wrap chord. The
% final straight segment must clear the polygon (sampled legality test,
% tolerance ~ half the 3 mm sample spacing -- documented in the README).
N = size(poly3, 1);
A2 = A3(1:2);
poly2 = poly3(:, 1:2);
B2 = poly2(iB, :);
lenDirect = norm(A3 - poly3(iB, :));
if fr_segClear(poly3(iB, 1:2), A2, poly2)
    lenTotal = lenDirect;
    iExit = 0;
    path3 = [poly3(iB, :); A3];
    wrapped = false;
    sweepAng = 0;
    return
end
wrapped = true;
edges = sqrt(sum((circshift(poly2, -1) - poly2).^2, 2));
cumCw = cumsum(edges);                       % CCW cumulative from vertex 1
ringLen = cumCw(end);
best = [inf, -1, 0];
for dir = [1, -1]
    for j = 1:N
        if j == iB
            continue
        end
        if dir == 1
            if j > iB
                L = cumCw(j) - cumCw(iB);
            else
                L = (ringLen - cumCw(iB)) + cumCw(j);
            end
        else
            if j < iB
                L = cumCw(iB) - cumCw(j);
            else
                L = cumCw(iB) + (ringLen - cumCw(j));
            end
        end
        if ~fr_segClear(poly2(j, :), A2, poly2)
            continue
        end
        cost = L + norm(A3 - poly3(j, :));
        if cost < best(1)
            best = [cost, j, dir];
        end
    end
end
lenTotal = inf;
iExit = -1;
path3 = zeros(0, 3);
sweepAng = 0;
if best(2) > 0
    j = best(2);
    dir = best(3);
    lenTotal = best(1);
    iExit = j;
    sweepAng = 2 * pi * abs(j - iB) / N;   % swept arc estimate
    if dir == 1
        if j > iB
            idx = (iB:j)';
        else
            idx = [iB:N, 1:j]';
        end
    else
        if j < iB
            idx = (iB:-1:j)';
        else
            idx = [iB:-1:1, N:-1:j]';
        end
    end
    path3 = [poly3(idx, :); A3];
end
end

function ok = fr_segClear(p2, q2, poly2)
% sampled legality: the straight segment from the plane exit to the
% projected attachment must stay outside the cam polygon. 41 samples;
% even-odd inside test per sample (documented approximation, tolerance
% ~ half the 2 deg polygon spacing). The exit->attachment 3D clearance
% against registered obstacles is check 7 (fr_pathClearance).
n = 41;
t = linspace(0, 1, n);
pts2 = p2 + t' * (q2 - p2);
inside = fr_inpoly(pts2, poly2);
ok = ~any(inside);
end

function inside = fr_inpoly(pts, poly)
% even-odd ray cast to +x, vectorized points x edges
n = size(pts, 1);
m = size(poly, 1);
j = circshift(poly, -1);
inside = false(n, 1);
for i = 1:n
    px = pts(i, 1);
    py = pts(i, 2);
    cnt = 0;
    for k = 1:m
        y1 = poly(k, 2);
        y2 = j(k, 2);
        if (y1 > py) ~= (y2 > py)
            xint = poly(k, 1) + (py - y1) * (j(k, 1) - poly(k, 1)) / ...
                (y2 - y1);
            if xint > px
                cnt = cnt + 1;
            end
        end
    end
    inside(i) = mod(cnt, 2) == 1;
end
end

function clr = fr_pathClearance(C, path3, kin, q)
% min SDF distance of the path samples to the registered CAD obstacles
% (clearance-only use; fr_obstacles_v2).
ns = size(path3, 1);
% transform path samples (femur frame) to ground
Tf = fr_bodyT(kin, q, 'femur_r');
w = zeros(ns, 3);
for i = 1:ns
    w(i, :) = (Tf(1:3, 1:3) * path3(i, :)' + Tf(1:3, 4))';
end
[~, ~, perBody] = fr_obstacles_v2(C, w, kin, q);
clr = min([perBody.dist], [], 'all');
end

function k = fr_kmaxF(C, diaStr)
if strcmp(diaStr, '10')
    k = C.maxContractionPct.d10 / 100;
else
    k = C.maxContractionPct.d20 / 100;
end
end

function p = fr_applyT(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end
