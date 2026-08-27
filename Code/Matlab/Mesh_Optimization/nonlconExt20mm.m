function [c, ceq] = nonlconExt20mm(x, ctx)
% nonlconExt20mm
%
% Nonlinear geometry constraints for the 20 mm extensor route.
%
% c <= 0 is feasible.
%
% x = [p1(1:3), pEnd(1:3), rest, tendon]

pEnd = x(4:6);
tendon = x(8);

geo = ctx.geo;

ceq = [];


%% Geometry-dependent tendon maximum

[~, tendonMaxGeom, ~, ~, cInside] = ...
    tendonLimit20mm(pEnd, geo);

[p8T1, okP8] = currentRouteP8T1(x, ctx);

if okP8
    cP8MinX = p8MinusXLimit20mm(p8T1, geo);
else
    cP8MinX = 1;
end


%% Constraints
%
% 1. Tendon cannot extend farther than the point where its straight
%    line reaches the BPA clearance envelope.
%
% 2. There must be enough geometrically available length for the
%    minimum 25 mm tendon.
%
% 3. If a circular tangent is being used, pEnd must permit a real
%    external tangent.

cTendonMax = tendon - tendonMaxGeom;

cMinimumPossible = geo.tendonMin - tendonMaxGeom;

c = [ ...
    cTendonMax;
    cMinimumPossible;
    cInside;
    cP8MinX];

end


function [p8T1, ok] = currentRouteP8T1(x, ctx)
% Use the route builder so this nonlinear constraint sees the same adjusted
% p8 that the prediction model uses for the current pEnd candidate.

p8T1 = [NaN NaN];
ok = false;

try
    [~, ~, routeInfo] = buildDistalRingLocation20mm( ...
        x(1:3), x(4:6), x(8), ctx);

    p8T1 = routeInfo.seed20(8,1:2);
    ok = all(isfinite(p8T1));
catch
    ok = false;
end

end


function c = p8MinusXLimit20mm(p8T1, geo)
% Positive c means p8 is too far in the -X direction.  The limit is a
% horizontal offset from the sloped CAD reference line in the t1 frame.

P1 = geo.p8MinXLineT1(1,:);
P2 = geo.p8MinXLineT1(2,:);

if abs(P2(2) - P1(2)) < 1e-12
    xLine = max(P1(1), P2(1));
else
    t = (p8T1(2) - P1(2)) / (P2(2) - P1(2));
    xLine = P1(1) + t*(P2(1) - P1(1));
end

xMin = xLine - geo.p8MinXMinusOffset;
c = xMin - p8T1(1);

end
