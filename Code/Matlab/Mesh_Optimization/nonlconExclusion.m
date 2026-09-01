function [c, ceq, info] = nonlconExclusion(x, geo, ctx, idxP2)
%NONLCONEXCLUSION Attachment and inflated-BPA collision constraints.
%
% c(1) <= 0: p2 clears the simplified proximal-tibia solid.
% c(2) <= 0: the BPA/fittings and tendon clear that solid at +5 deg.
% c(3) <= 0: the BPA/fittings and tendon clear the femoral ellipse.
% c(4) <= 0: tendon + two fittings leave nonnegative current BPA length.
%% c(5) <= 0: the fixed pWrap y-z line intersects the BPA contact surface.

% The complete series path follows the same p1-wrap-p2 polyline used by the
% predictor.  At +5 deg the wrap is normally active.  The tendon occupies
% the final portion of that polyline and the BPA plus fittings occupy the
% proximal portion.

if nargin < 4 || isempty(idxP2)
    idxP2 = 4:6;
end

p1 = x(1:3);
p2 = x(idxP2);
p1 = p1(:).';
p2 = p2(:).';

rest = x(7);
tendon = x(8);

ceq = [];

% Use the sampled knee position nearest the requested +5 deg collision
% configuration.  p1 is stored in the femur frame; p2 is already in t1.
[~, idxCollision] = min(abs(ctx.phiD - geo.collisionAngleD));
[~, ~, routeInfo] = buildKneeFlexorRoute20mm(p1, p2, ctx);
cWrap = routeInfo.wrapContactConstraint;
p1T1 = routeInfo.p1T1(idxCollision,:);

if routeInfo.active(idxCollision)
    routeT1 = [p1T1; routeInfo.pWrapT1(idxCollision,:); p2];
else
    routeT1 = [p1T1; p2];
end

segmentLength = vecnorm(diff(routeT1,1,1),2,2);
pathLength = sum(segmentLength);

info = struct;
info.angleD = ctx.phiD(idxCollision);
info.p1T1 = p1T1;
info.p2T1 = p2;
info.pWrapT1 = routeInfo.pWrapT1(idxCollision,:);
info.wrapActive = routeInfo.active(idxCollision);
info.bpaRadius = geo.bpaRadius;
info.tendonRadius = geo.tendonRadius;

if ~isfinite(pathLength) || pathLength <= eps || ...
        ~isfinite(rest) || rest <= 0 || ...
        ~isfinite(tendon) || tendon < 0
    c = [1; 1; 1; 1; 1];
    info.minClearance = -Inf;
    info.constraintMargin = -Inf;
    info.failReason = "Invalid collision geometry";
    return
end

% The tendon always terminates at p2.  The junction moves toward p2 when
% the optimized tendon becomes shorter.  Everything proximal to the
% junction consists of the current BPA length plus both fittings and uses
% the inflated BPA radius for collision clearance.
tendonLengthChecked = min(tendon, pathLength);
bpaFittingsLength = pathLength - tendonLengthChecked;
currentMuscleLength = bpaFittingsLength - 2*ctx.fitting;

nPath = max(2, ceil(pathLength/geo.centerlineSampleStep) + 1);
pathPosition = linspace(0, pathLength, nPath).';
pathPosition = unique([pathPosition; bpaFittingsLength]);
centerlineT1 = pointsAlongPolyline(routeT1, pathPosition);
junctionT1 = pointsAlongPolyline(routeT1, bpaFittingsLength);

isBPA = pathPosition <= bpaFittingsLength + 10*eps;
pointRadius = geo.tendonRadius*ones(size(pathPosition));
pointRadius(isBPA) = geo.bpaRadius;
component = repmat("tendon", numel(pathPosition), 1);
component(isBPA) = "BPA + fittings";

% Signed distance is positive outside the simplified tibia, zero on its
% surface, and negative inside it.
[signedDistanceTibia, region] = signedDistanceToTibia20mm(centerlineT1, geo);

% Transform the complete wrapped centerline into the femur frame.
centerlineFemur = zeros(size(centerlineT1));
for i = 1:size(centerlineT1,1)
    qICR = RowVecTrans(ctx.T_ICR_t1(:,:,idxCollision), centerlineT1(i,:));
    centerlineFemur(i,:) = RowVecTrans(ctx.T_Pam(:,:,idxCollision), qICR);
end

signedDistanceFemur = signedDistanceToFemur(centerlineFemur, geo);

% A signed clearance of zero means the inflated BPA touches the bone.
% Half a sample interval is added to the required clearance because signed
% distance is 1-Lipschitz; this prevents a collision from hiding between
% adjacent centerline samples.
sampleAllowance = 0.5*geo.centerlineSampleStep;
surfaceClearanceTibia = signedDistanceTibia - pointRadius;
surfaceClearanceFemur = signedDistanceFemur - pointRadius;
[minClearanceTibia, idxWorstTibia] = min(surfaceClearanceTibia);
[minClearanceFemur, idxWorstFemur] = min(surfaceClearanceFemur);
cTibia = geo.clearance + sampleAllowance;
cFemur = geo.clearance + sampleAllowance - minClearanceFemur;
cSeries = -currentMuscleLength;

% Preserve an explicit attachment-point exclusion in addition to the BPA
% body constraint.  This uses the same 3-D union of the two simplified
% tibia cases rather than only the cross-section selected at p2.y.
signedDistanceP2 = signedDistanceToTibia20mm(p2, geo);
cP2 = geo.clearance - signedDistanceP2;

c = [cP2; cTibia; cFemur; cSeries; cWrap];

info.bpaStartT1 = p1T1;
info.bpaEndT1 = junctionT1;
info.junctionT1 = junctionT1;
info.bpaFittingsLengthChecked = bpaFittingsLength;
info.currentMuscleLength = currentMuscleLength;
info.tendon = tendon;
info.tendonLengthChecked = tendonLengthChecked;
info.minClearance = min(minClearanceTibia, minClearanceFemur);
info.minClearanceTibia = minClearanceTibia;
info.minClearanceFemur = minClearanceFemur;
info.requiredClearance = geo.clearance;
info.sampleAllowance = sampleAllowance;
info.constraintMargin = -max(cTibia, cFemur);
info.worstTibiaCenterT1 = centerlineT1(idxWorstTibia,:);
info.worstTibiaRegion = char(region(idxWorstTibia));
info.worstTibiaComponent = char(component(idxWorstTibia));
info.worstTibiaRadius = pointRadius(idxWorstTibia);
info.worstFemurCenterFemur = centerlineFemur(idxWorstFemur,:);
info.worstFemurComponent = char(component(idxWorstFemur));
info.worstFemurRadius = pointRadius(idxWorstFemur);
info.cP2 = cP2;
info.cTibia = cTibia;
info.cFemur = cFemur;
info.cSeries = cSeries;
info.failReason = "";
info.cWrap = cWrap;

end

function points = pointsAlongPolyline(nodes, distance)
%POINTSALONGPOLYLINE Evaluate one or more arclength positions.

nodes = reshape(nodes, [], 3);
distance = distance(:);
segmentLength = vecnorm(diff(nodes,1,1),2,2);
arc = [0; cumsum(segmentLength)];

if arc(end) <= eps
    points = repmat(nodes(1,:), numel(distance), 1);
    return
end

distance = min(arc(end), max(0, distance));
points = interp1(arc, nodes, distance, 'linear');

end

function sdFemur = signedDistanceToFemur(points, geo)
%SIGNEDDISTANCETOFEMUR Signed distance to the finite condyle ellipse.
% The physical ellipse occupies the femur X-Y plane and is extruded across
% the proximal-tibia Z width supplied by the user.

points = reshape(points, [], 3);
x = points(:,1);
y = points(:,2);
z = points(:,3);

sdEllipse2D = signedDistanceToPolygon( ...
    x, y, geo.femurEllipseX, geo.femurEllipseY);

% Exact signed distance to the Z slab: negative inside, positive outside.
sdZSlab = max(geo.femurZMin - z, z - geo.femurZMax);
sdFemur = intersectSignedDistances(sdEllipse2D, sdZSlab);

end

function [sdUnion, region] = signedDistanceToTibia(points, geo)
%SIGNEDDISTANCETOTIBIA Signed distance to the union of:
%   1) the outer-profile prism for y >= yCut, and
%   2) the inner-circle cylinder for y <= yCut.

points = reshape(points, [], 3);
x = points(:,1);
y = points(:,2);
z = points(:,3);

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
