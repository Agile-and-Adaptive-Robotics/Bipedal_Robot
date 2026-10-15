function out = routingChainLib(fn, varargin)
% Routing Chain Lib
% Campaign: routing27 (2026-10-03)
% Geometry helpers shared by the routing campaign: chain-edge hinge
% transforms (rolling-knee aware), route segment mapping into the
% proximal frame and the global pelvis frame, segment-segment distance
% (Ericson), bone-cloud clearance, and the lateral-footprint metric.
%
% Dispatch-style library so every helper lives in one place:
%   routingChainLib('edgeTrans', edge, thetaRad, kneePv)
%   routingChainLib('routeSegments', points, cross, T)
%   routingChainLib('segSegDist', p1, q1, p2, q2)
%   routingChainLib('cloudMinDist', samples, clouds)
%   routingChainLib('segListToPelvis', S, frame1Body, chain)
%
% CONVENTIONS (documented simplifications of the campaign):
%   - Each actuator sweeps only its PRIMARY DOF; every other joint sits
%     at its default (0) angle, matching the human-target generation.
%   - Cross-actuator spacing/footprint use the GLOBAL PELVIS frame with
%     each route swept over its own grid (an envelope approximation:
%     same-joint pairs are exact, cross-joint pairs are conservative,
%     the same style as the committed batch's route-clearance check).
%   - The lateral axis is z of the chain frames at default pose (the
%     gait2392 sagittal plane is x-y; +z is mediolateral).

switch fn
    case 'edgeTrans'
        out = edgeTrans(varargin{:});
    case 'routeSegments'
        out = routeSegments(varargin{:});
    case 'segSegDist'
        out = segSegDist(varargin{:});
    case 'cloudMinDist'
        out = cloudMinDist(varargin{:});
    case 'segListToPelvis'
        out = segListToPelvis(varargin{:});
    otherwise
        error('routingChainLib:unknown', 'Unknown helper %s', fn)
end
end

%% ---------------------------------------------------------------------
function T = edgeTrans(edge, thetaRad, kneePv)
% Hinge transform about an edge's axis/pivot. edge: struct(axis, pivot,
% name, transTheta, transX, transY). kneePv: optional override pivot for
% the rolling knee (x/y from the splines at thetaRad; z kept).
axis = edge.axis(:) / norm(edge.axis);
pv = edge.pivot(:);
hasRoll = isfield(edge, 'transThetaX') && ~isempty(edge.transThetaX);
hasRoll = hasRoll || (isfield(edge, 'transThetaY') && ~isempty(edge.transThetaY));
if hasRoll && nargin >= 3 && ~isempty(kneePv)
    px = 0;
    py = 0;
    if isfield(edge, 'transThetaX') && ~isempty(edge.transThetaX)
        px = interp1(edge.transThetaX, edge.transX, thetaRad, ...
            'linear', 'extrap');
    end
    if isfield(edge, 'transThetaY') && ~isempty(edge.transThetaY)
        py = interp1(edge.transThetaY, edge.transY, thetaRad, ...
            'linear', 'extrap');
    end
    pv = [px, py, pv(3)];
end
Km = [0, -axis(3), axis(2); axis(3), 0, -axis(1); -axis(2), axis(1), 0];
R = eye(3) + sin(thetaRad) * Km + (1 - cos(thetaRad)) * (Km * Km);
T = [R, pv(:); 0 0 0 1];
end

%% ---------------------------------------------------------------------
function S = routeSegments(points, cross, T)
% All route segments at every grid cell of the 4x4xN1xN2 transform stack,
% expressed in the PROXIMAL frame: rows [pA pB] (6 columns), ordered
% cell-major (ii fastest, then iii), the committed driver's layout.
[N1, N2] = size(T, [3, 4]);
nPts = size(points, 1);
S = zeros(N1 * N2 * (nPts - 1), 6);
row = 0;
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for iii = 1:N2
        T2 = T(:, :, iii, 2);
        for r = 1:nPts - 1
            row = row + 1;
            S(row, 1:3) = mapRow(points, cross, r, T1, T2, 1);
            S(row, 4:6) = mapRow(points, cross, r + 1, T1, T2, 1);
        end
    end
end
end

function p = mapRow(points, cross, r, T1, T2, frame)
% Twin of BiPam_pulley's frame chaining for one row (rigid geometry).
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
        p = RowVecTransLocal(T2, p);
        f = 2;
    else
        p = RowVecTransLocal(T1, p);
        f = 1;
    end
end
end

function v = RowVecTransLocal(T, v)
v = (T(1:3, 1:3) * v.' + T(1:3, 4)).';
end

%% ---------------------------------------------------------------------
function d = segSegDist(p1, q1, p2, q2)
% Minimum distance between two 3-D segments (Ericson, Real-Time
% Collision Detection), ported verbatim from Opt_run_BiPulley.m.
d1 = q1 - p1;
d2 = q2 - p2;
r = p1 - p2;
a = dot(d1, d1);
e = dot(d2, d2);
f = dot(d2, r);
c = dot(d1, r);
b = dot(d1, d2);
if a <= 1e-18 && e <= 1e-18
    d = norm(r);
    return
elseif a <= 1e-18
    s = 0;
    t = min(1, max(0, f / e));
elseif e <= 1e-18
    t = 0;
    s = min(1, max(0, -c / a));
else
    denom = a * e - b * b;
    if denom > 1e-18
        s = min(1, max(0, (b * f - c * e) / denom));
    else
        s = 0;
    end
    t = (b * s + f) / e;
    if t < 0
        t = 0;
        s = min(1, max(0, -c / a));
    elseif t > 1
        t = 1;
        s = min(1, max(0, (b - c) / a));
    end
end
d = norm((p1 + s * d1) - (p2 + t * d2));
end

%% ---------------------------------------------------------------------
function dmin = cloudMinDist(samples, clouds)
% Min distance from sample points (M x 3) to a union of point clouds
% (cell array of P x 3). Surface-sampled clouds give distance-to-surface
% up to the sample spacing (documented); the constraint margin covers it.
samples = samples(: , 1:3);
dmin = inf;
for k = 1:numel(clouds)
    C = clouds{k};
    if isempty(C)
        continue
    end
    % Squared distances without a toolbox: ||s||^2 + ||c||^2 - 2 s.c
    ss = sum(samples .^ 2, 2);
    cc = sum(C .^ 2, 2)';
    D2 = max(ss + cc - 2 * (samples * C.'), 0);
    dmin = min(dmin, sqrt(min(D2(:))));
    if dmin < 1e-4
        return   % deep penetration, no point continuing
    end
end
end

%% ---------------------------------------------------------------------
function S = segListToPelvis(S, frame1Body, chain)
% Map frame-1 segments into the GLOBAL PELVIS frame by chaining the
% DEFAULT-pose transforms of the proximal chain edges above frame1Body.
% chain: struct with .edges (containers.Map 'parent>child' -> edge) and
% the chain order below.
switch frame1Body
    case 'pelvis'
        return
    case 'torso'
        order = {'pelvis>torso'};
    case 'femur_r'
        order = {'pelvis>femur_r'};
    case 'tibia_r'
        order = {'pelvis>femur_r', 'femur_r>tibia_r'};
    case 'talus_r'
        order = {'pelvis>femur_r', 'femur_r>tibia_r', 'tibia_r>talus_r'};
    case 'calcn_r'
        order = {'pelvis>femur_r', 'femur_r>tibia_r', 'tibia_r>talus_r', ...
            'talus_r>calcn_r'};
    case 'toes_r'
        order = {'pelvis>femur_r', 'femur_r>tibia_r', 'tibia_r>talus_r', ...
            'talus_r>calcn_r', 'calcn_r>toes_r'};
    otherwise
        error('routingChainLib:badBody', 'Unhandled frame1 body %s', ...
            frame1Body)
end
T = eye(4);
for k = 1:numel(order)
    key = order{k};
    if ~isKey(chain.edges, key)
        error('routingChainLib:noEdge', 'Chain edge %s missing', key)
    end
    T = chain.edges(key).Tdefault * T;
end
if isempty(S)
    return
end
S(:, 1:3) = (T(1:3, 1:3) * S(:, 1:3).' + T(1:3, 4)).';
S(:, 4:6) = (T(1:3, 1:3) * S(:, 4:6).' + T(1:3, 4)).';
end
