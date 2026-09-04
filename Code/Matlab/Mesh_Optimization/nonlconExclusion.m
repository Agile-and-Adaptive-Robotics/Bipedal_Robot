function [c, ceq, info] = nonlconExclusion(x, geo, ctx, idxP2)
%NONLCONEXCLUSION Attachment and inflated-BPA collision constraints.
%
% c(1) <= 0: p2 clears the simplified proximal-tibia solid.
% c(2) <= 0: at +5 deg, the straight route stays outside the hard tibia by
%              bpaRs for BPA/fittings and tendonRadius for the tendon.
% c(3) <= 0: the BPA/fittings and tendon clear the femoral ellipse.
% c(4) <= 0: tendon + two fittings leave nonnegative current BPA length.
% c(5) <= 0: while wrap is active, the prescribed pWrap y-z line reaches
%              the nominal bpaRb standoff surface used to place pWrap.
% c(6) <= 0: tendon does not extend proximally beyond pWrap.
% c(7) <= 0: relatives strain between 0 and 1 for robot range of motion

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
[Location, ~, routeInfo] = buildKneeFlexorRoute20mm(p1, p2, tendon, ctx);
N = ctx.N;  %number of increments
M = size(Location,1)-1;         %Number of rows minus 1 since we're taking the difference
SL = zeros(N,M);

for ii = 1:N                    %Repeat for each orientation
    for i = 1:M               %Calculate all segments
        pointA = Location(i,:,ii);
        pointB = Location(i+1,:,ii);
        if i+1 == ctx.CrossPoint
            pointB = RowVecTrans(ctx.T_Pam(:,:,ii), pointB);
        end
        SL(ii,i) = norm(pointA - pointB);   %segment lengths
    end
end
Lmt = sum(SL,2);            % Length of musculotendon
Lm = Lmt - ctx.Xi0 - tendon - 2*ctx.fitting;    %BPA length
relstrain = ((rest-Lm)./rest)./ctx.KMAX;      %relative strain
idxRoM = ctx.phiD >= -120 & ctx.phiD <=5 ;   %RoM robot can actually perform
relRoM = relstrain(idxRoM);

cRelStrain = max([ ...
    max(relRoM - 1); ...   % violation of upper bound relstrain <= 1
    max(-relRoM)]);        % violation of lower bound relstrain >= 0

% Maximum tendon length that can fit entirely on the distal
% pWrap -> p2 segment.  Longer tendon would extend proximally
% past pWrap and require a different collision model.
tendonSL = SL((idxRoM & SL(:,2) > 0), 2);   %relevant nonzero tendon segment lengths

% <= 0 is feasible.
if isempty(tendonSL)
    cTendonWrap = 0;            %defensive against no wrapping point
    tendonMax = NaN;
else
    tendonMax = min(tendonSL) - ctx.fitting;  %The fitting cannot bend
    cTendonWrap = tendon - tendonMax;
end

cWrap = routeInfo.wrapContactConstraint;

p1T1 = routeInfo.p1T1(idxCollision,:);           %p1 in t1 at collision angle
pWrapT1 = routeInfo.pWrapT1(idxCollision,:);     %wrap point in t1 at collision angle

if SL(idxCollision,2) > 0
    routeT1 = [p1T1; pWrapT1; p2];
else
    routeT1 = [p1T1; p2];
end

pathLength = Lmt(idxCollision);
info = struct;
info.angleD = ctx.phiD(idxCollision);
info.p1T1 = p1T1;
info.p2T1 = p2;
info.pWrapT1 = pWrapT1;
info.wrapActive = routeInfo.active(idxCollision);
% Keep the three BPA radii explicit in diagnostics so their different roles
% cannot be confused later.
info.bpaRb = geo.bpaRb;        % nominal pWrap standoff from hard tibia
info.bpaRadius = geo.bpaRb;     % legacy diagnostic alias; nominal radius only
info.bpaRs = geo.bpaRs;        % minimum BPA centerline clearance to hard tibia
info.wRap = geo.wRap;          % Xi3 bend-length radius; not a collision radius
info.tendonRadius = geo.tendonRadius;

if ~isfinite(pathLength) || pathLength <= eps || ...
        ~isfinite(rest) || rest <= 0 || ...
        ~isfinite(tendon) || tendon < 0
    c = ones(7,1);
    info.minClearance = -Inf;
    info.constraintMargin = -Inf;
    info.failReason = "Invalid collision geometry";
    return
end

% The tendon always terminates at p2.  Everything proximal to the
% BPA/tendon junction is represented by the BPA + fittings portion.
% The actual optimized tendon length is used; it is not clipped here.

% Arc-length location, measured from p1, where the BPA/fittings end and the
% distal tendon begins along the straight-segment route.
bpaFittingsLength = pathLength - tendon;

% Modeled BPA length already calculated above from the same SL/Lmt route.
currentMuscleLength = Lm(idxCollision);

% Sample the straight p1 -> pWrap -> p2 approximation at about 1 mm spacing.
% This sampling is ONLY for collision checking; it does not change SL, Lmt,
% pWrap, the bend angle, or the Xi3 bend calculation.
nPath = max(2, ceil(pathLength/geo.centerlineSampleStep) + 1);
pathPosition = linspace(0, pathLength, nPath).';

% Force the BPA/tendon junction to be one of the sampled points.
pathPosition = unique([pathPosition; bpaFittingsLength]);
centerlineT1 = pointsAlongPolyline(routeT1, pathPosition);
junctionT1 = pointsAlongPolyline(routeT1, bpaFittingsLength);

% Identify which centerline samples belong to the BPA/fittings versus tendon.
isBPA = pathPosition <= bpaFittingsLength + 10*eps;
component = repmat("tendon", numel(pathPosition), 1);
component(isBPA) = "BPA + fittings";

% Tibia and femur intentionally use different BPA radii.  The tibial wrap
% point is nominally placed bpaRb = 20 mm from the hard surface, but the
% straight-line approximation is allowed to approach to bpaRs = 16 mm to
% represent changing BPA diameter/local squish.  The femur has no intentional
% wrap contact, so it retains the nominal bpaRb envelope.
tibiaRadius = geo.tendonRadius*ones(size(pathPosition));
tibiaRadius(isBPA) = geo.bpaRs;
femurRadius = geo.tendonRadius*ones(size(pathPosition));
femurRadius(isBPA) = geo.bpaRs;

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

% Convert centerline-to-hard-surface distances into component surface
% clearances.  For the tibia, zero means the BPA has reached the allowed
% compressed radius bpaRs (or the tendon has reached tendonRadius).
% For the femur, zero uses the nominal bpaRb envelope.
surfaceClearanceTibia = signedDistanceTibia - tibiaRadius;
surfaceClearanceFemur = signedDistanceFemur - femurRadius;
[minClearanceTibia, idxWorstTibia] = min(surfaceClearanceTibia);
[minClearanceFemur, idxWorstFemur] = min(surfaceClearanceFemur);

% Tibia: no additional stand-off is imposed beyond bpaRs/tendonRadius.
% Femur: retain the existing 1 mm rigid clearance plus half a sample spacing
% so a collision cannot hide between adjacent centerline samples.
sampleAllowance = 0.5*geo.centerlineSampleStep;
cTibia = -minClearanceTibia;
cFemur = geo.clearance + sampleAllowance - minClearanceFemur;
cSeries = -Lm(idxCollision);

% Preserve an explicit attachment-point exclusion in addition to the BPA
% body constraint.  This uses the same 3-D union of the two simplified
% tibia cases rather than only the cross-section selected at p2.y.
signedDistanceP2 = signedDistanceToTibia20mm(p2, geo);
cP2 = geo.clearance - signedDistanceP2;

c = [cP2; ...
     cTibia; ...
     cFemur; ...
     cSeries; ...
     cWrap; ...
     cTendonWrap; ...
     cRelStrain];

info.tendonMax = tendonMax;
info.cTendonWrap = cTendonWrap;
info.bpaStartT1 = p1T1;
info.bpaEndT1 = junctionT1;
info.junctionT1 = junctionT1;
info.bpaFittingsLengthChecked = bpaFittingsLength;
info.currentMuscleLength = currentMuscleLength;
info.tendon = tendon;
info.tendonLengthChecked = tendon;  % diagnostic only; tendon is not clipped
info.minClearance = min(minClearanceTibia, minClearanceFemur);
info.minClearanceTibia = minClearanceTibia;
info.minClearanceFemur = minClearanceFemur;
info.tibiaBPAMinRadius = geo.bpaRs;
info.femurBPANominalRadius = geo.bpaRb;
info.requiredFemurClearance = geo.clearance;
info.requiredClearance = geo.clearance;  % legacy diagnostic alias; femur/p2 only
info.sampleAllowance = sampleAllowance;
info.constraintMargin = -max(cTibia, cFemur);
info.worstTibiaCenterT1 = centerlineT1(idxWorstTibia,:);
info.worstTibiaRegion = char(region(idxWorstTibia));
info.worstTibiaComponent = char(component(idxWorstTibia));
info.worstTibiaRadius = tibiaRadius(idxWorstTibia);
info.worstFemurCenterFemur = centerlineFemur(idxWorstFemur,:);
info.worstFemurComponent = char(component(idxWorstFemur));
info.worstFemurRadius = femurRadius(idxWorstFemur);
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
z = abs(points(:,3));   % Mirror femoral exclusion geometry about z = 0

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
