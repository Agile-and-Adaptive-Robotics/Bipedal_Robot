function [pTan, tendonMax, targetRadius, targetType, cInside] = tendonLimit20mm(pEnd, geo)
% tendonLimit20mm
%
% Geometry-dependent upper tendon length for the 20 mm extensor route.
%
% pEnd is in t1 frame.
%
% targetType:
%   1 = lower expanded cylinder
%   2 = vertical tangent line
%   3 = upper expanded cylinder
%
% cInside <= 0 is feasible.
%
% The tangent length is used instead of the radial distance because the
% tendon/BPA path approaches the clearance envelope tangentially.

P = pEnd(1:2);

yLower = geo.tibiaLowerCenter(2);
yUpper = geo.tibiaUpperCenter(2);

cInside = -1;


%% Below lower cylinder

if P(2) < yLower

    C = geo.tibiaLowerCenter;
    R = geo.tibiaLowerClearRadius;

    d = norm(P - C);

    cInside = R - d;

    [pTan, ok] = tangentPointPositiveX(P, C, R);

    if ~ok
        pTan = [geo.tibiaWallX, yLower];
        tendonMax = 0;
    else
        tendonMax = norm(pTan - P);
    end

    targetRadius = R;
    targetType = 1;


%% Between the two cylinders: vertical tangent region

elseif P(2) <= yUpper

    pTan = [geo.tibiaWallX, P(2)];

    tendonMax = abs(geo.tibiaWallX - P(1));

    targetRadius = 0;
    targetType = 2;


%% Above upper cylinder

else

    C = geo.tibiaUpperCenter;
    R = geo.tibiaUpperClearRadius;

    d = norm(P - C);

    cInside = R - d;

    [pTan, ok] = tangentPointPositiveX(P, C, R);

    if ~ok
        pTan = [geo.tibiaWallX, yUpper];
        tendonMax = 0;
    else
        tendonMax = norm(pTan - P);
    end

    targetRadius = R;
    targetType = 3;

end


end


function [q, ok] = tangentPointPositiveX(P, C, R)
% External tangent point from P to circle C,R.
% Select tangent point furthest in +x.

v = P - C;
d2 = dot(v,v);

if d2 <= R^2
    q = [NaN NaN];
    ok = false;
    return
end

base = C + (R^2/d2)*v;

h = R*sqrt(d2 - R^2)/d2;

vp = [-v(2), v(1)];

q1 = base + h*vp;
q2 = base - h*vp;

if q1(1) >= q2(1)
    q = q1;
else
    q = q2;
end

ok = true;

end