function out = fr_evalRoute(route, kin, q)
%FR_EVALROUTE Evaluate a parsed route's declared behaviour at pose q.
%
% Implements, per the geometry extraction's protocol (geometry_20261103
% README): static points at their declared owning-body coordinates,
% ConditionalPathPoints by declared-range membership at the realized
% coordinate value, MovingPathPoints by their declared SimmSpline location
% functions - NOTHING frozen at the default pose.
%
% out.pts         struct array of active points: .name .body .loc (1x3,
%                 owning-body frame) .world (1x3 ground) .active
% out.activeMask  1xN logical over the route's declared points
% out.lengthM     polyline length over ACTIVE points (ground frame), m
% out.transitions conditional-point activation states (for the
%                 wrapping-transition checks; gait2392 has no wrap objects,
%                 the conditional points are the functional equivalent)

n = numel(route.points);
out.pts = repmat(struct('name', '', 'body', '', 'loc', [NaN NaN NaN], ...
    'world', [NaN NaN NaN], 'active', false), 1, n);
out.activeMask = false(1, n);
out.lengthM = NaN;
for i = 1:n
    P = route.points{i};
    pt = out.pts(i);
    pt.name = P.name;
    pt.body = P.body;
    act = true;
    loc = P.loc(:)';
    if isfield(P, 'range') && ~isempty(P.range)
        qv = q.(P.coord);
        act = (qv >= P.range(1) && qv <= P.range(2));
    end
    if isfield(P, 'fnXYZ')
        lx = fr_evalFn(P.fnXYZ.x, q.(P.coordXYZ{1}));
        ly = fr_evalFn(P.fnXYZ.y, q.(P.coordXYZ{2}));
        lz = fr_evalFn(P.fnXYZ.z, q.(P.coordXYZ{3}));
        loc = [lx, ly, lz];
    end
    pt.loc = loc;
    pt.active = act;
    if act
        T = fr_bodyT(kin, q, P.body);
        pt.world = (T(1:3, 1:3) * loc(:) + T(1:3, 4))';
        out.activeMask(i) = true;
    else
        pt.world = [NaN NaN NaN];
    end
    out.pts(i) = pt;
end
w = vertcat(out.pts(out.activeMask).world);
if size(w, 1) >= 2
    out.lengthM = sum(sqrt(sum(diff(w).^2, 2)));   % scalar polyline length
end
% conditional transitions record
out.transitions = repmat(struct('name', '', 'coord', '', 'range', ...
    [0 0], 'active', false, 'coordValue', 0), 0, 1);
for i = 1:n
    P = route.points{i};
    if isfield(P, 'range') && ~isempty(P.range)
        out.transitions(end+1) = struct('name', P.name, ...
            'coord', P.coord, 'range', P.range, ...
            'active', out.activeMask(i), 'coordValue', q.(P.coord));
    end
end
end

function v = fr_evalFn(F, xq)
switch F.type
    case 'linear'
        v = F.coeffs(1) * xq + F.coeffs(2);
    case 'constant'
        v = F.value;
    case 'spline'
        v = fr_splinemm(F.x, F.y, xq);
    otherwise
        error('fr_evalRoute:FnType', 'unknown function type %s', F.type);
end
end
