function [c, ceq] = nonlconExclusion(x, geo, idxP2)
%NONLCONEXCLUSION Nonlinear exclusion constraint for p2.
%
% c <= 0 is feasible.
% c > 0 means p2 is inside excluded object.

p2 = x(idxP2);
p2 = p2(:).';

px = p2(1);
py = p2(2);
pz = p2(3);

ceq = [];

% -------------------------------------------------------------------------
% Case 2:
% py < -8.87 mm
% Exclusion is the inner circle in X-Z.
% -------------------------------------------------------------------------
if py < geo.yCut
    dx = px - geo.circleCenter(1);
    dz = pz - geo.circleCenter(3);

    d = hypot(dx, dz);

    c = geo.circleRadius + geo.clearance - d;
    return
end

% -------------------------------------------------------------------------
% Case 1:
% py >= -8.87 mm
% Exclusion is the outer profile.
% Inside object means closer to origin than the profile boundary.
% -------------------------------------------------------------------------
theta = atan2(pz, px);
if theta < 0
    theta = theta + 2*pi;
end

rP = hypot(px, pz);

c = -1;

if theta >= min(geo.thetaVec) && theta <= max(geo.thetaVec)
    rBoundary = interp1(geo.thetaVec, geo.rVec, theta, 'pchip');

    c = rBoundary + geo.clearance - rP;
end

end