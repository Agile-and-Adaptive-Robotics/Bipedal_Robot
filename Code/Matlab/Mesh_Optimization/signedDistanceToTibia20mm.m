function [sdUnion, region] = signedDistanceToTibia20mm(points, geo)
%SIGNEDDISTANCETOTIBIA Signed distance to the union of:
%   1) the outer-profile prism for y >= yCut, and
%   2) the inner-circle cylinder for y <= yCut.

points = reshape(points, [], 3);
x = points(:,1);
y = points(:,2);
z = abs(points(:,3));   % Mirror tibial exclusion geometry about z = 0

% The original radial definition is exactly the polygon formed by the
% origin, the ordered outer profile, and the origin again.
polyX = [0; geo.xProf(:); 0];
polyZ = [0; geo.zProf(:); 0];
sdProfile2D = signedDistanceToPolygon(x, z, polyX, polyZ);

sdCircle2D = hypot(x - geo.circleCenter(1), ...
                   z - geo.circleCenter(3)) - geo.circleRadius;

% For an intersection of a 2-D region and a y half-space, combine the two
% signed distances using the standard extruded-solid SDF expression.
% Outer profile occupies y >= yCut; inner circle occupies y <= yCut.
sdOuter = intersectSignedDistances(sdProfile2D, geo.yCut - y);
sdInner = intersectSignedDistances(sdCircle2D, y - geo.yCut);

[sdUnion, whichRegion] = min([sdOuter, sdInner], [], 2);
region = strings(size(whichRegion));
region(whichRegion == 1) = "outer profile";
region(whichRegion == 2) = "inner circle";

end

function sd = intersectSignedDistances(sdA, sdB)
% Both component interiors use the convention sd <= 0.
outside = hypot(max(sdA, 0), max(sdB, 0));
inside = min(max(sdA, sdB), 0);
sd = outside + inside;
end

function sd = signedDistanceToPolygon(x, z, polyX, polyZ)
% Euclidean signed distance to a closed 2-D polygon.

inside = inpolygon(x, z, polyX, polyZ);
distance = inf(size(x));

for k = 1:(numel(polyX) - 1)
    ax = polyX(k);
    az = polyZ(k);
    bx = polyX(k+1);
    bz = polyZ(k+1);

    abx = bx - ax;
    abz = bz - az;
    denom = abx^2 + abz^2;

    if denom <= eps
        continue
    end

    t = ((x - ax)*abx + (z - az)*abz)/denom;
    t = min(1, max(0, t));

    qx = ax + t*abx;
    qz = az + t*abz;
    distance = min(distance, hypot(x - qx, z - qz));
end

sd = distance;
sd(inside) = -distance(inside);

end