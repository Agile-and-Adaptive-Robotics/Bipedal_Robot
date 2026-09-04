function [Location, bendMeasure, info] = buildKneeFlexorRoute20mm(p1, pEnd, tendon, ctx)
%BUILDKNEEFLEXORROUTE20MM Build the flexor route with one t1 wrap point.
%
% p1 is fixed in the femur frame. pEnd is fixed in the t1 frame. pWrap is
% an equivalent routing point used by the straight-segment model; it is not
% a solved tangent point. While contact is active, pWrap is placed at the
% nominal bpaRb standoff from the hard tibial surface. nonlconExclusion then
% permits the straight segments to approach only as close as bpaRs.
%
% The wrap is released during the extension-to-flexion sweep when the
% existing atan2 release rule is satisfied. An inactive wrap is represented
% by repeating pEnd so the second segment has zero length.

p1 = p1(:).';
pEnd = pEnd(:).';

if numel(p1) ~= 3 || numel(pEnd) ~= 3
    error('buildKneeFlexorRoute20mm:BadPointSize', ...
        'p1 and pEnd must each contain three coordinates.')
end

N = ctx.N;

p1T1 = zeros(N,3);
% aInD = zeros(N,1);
% aOutD = zeros(N,1);
turnWrapped = zeros(N,1);
aDirectD = zeros(N,1);
aWrapEndD = zeros(N,1);

% First express p1 in t1 at every knee angle. The combined femur-to-t1
% transform is precomputed once by buildKneeFlexorContext20mm.
for i = 1:N
    p1T1(i,:) = RowVecTrans(ctx.T_t1_f(:,:,i), p1);
end

%% Check whether the direct p1 -> pEnd route actually requires tibial wrap

directNeedsWrap = false(N,1);
directMinClearance = inf(N,1);

for i = 1:N

    directVector = pEnd - p1T1(i,:);
    directLength = norm(directVector);

    if directLength <= eps
        error('buildKneeFlexorRoute20mm:DegenerateDirectRoute', ...
            'Direct p1-pEnd route has zero length at sample %d.', i)
    end

    % Along the physical series path, the tendon occupies the distal
    % portion ending at pEnd. Everything proximal to it is BPA + fittings.
    bpaFittingsLengthDirect = max(0, directLength - tendon);

    nDirect = max(2, ...
        ceil(directLength/ctx.geo.centerlineSampleStep) + 1);

    sDirect = linspace(0, directLength, nDirect).';

    % Explicitly include the BPA/tendon junction.
    sDirect = unique([sDirect; bpaFittingsLengthDirect]);

    directPoints = p1T1(i,:) + ...
        (sDirect/directLength).*directVector;

    % Tendon is allowed to approach to its physical radius.
    directRadius = ...
        ctx.geo.tendonRadius*ones(size(sDirect));

    % BPA + fittings use the minimum allowable compressed BPA radius.
    isDirectBPA = ...
        sDirect <= bpaFittingsLengthDirect + 10*eps;

    directRadius(isDirectBPA) = ctx.geo.bpaRs;

    signedDistanceDirect = ...
        signedDistanceToTibia20mm(directPoints, ctx.geo);

    directClearance = signedDistanceDirect - directRadius;

    directMinClearance(i) = min(directClearance);

    % Zero means allowable contact. Negative means the direct path
    % penetrates beyond the permitted BPA/tendon envelope and needs wrap.
    directNeedsWrap(i) = directMinClearance(i) < 0;
end


% The wrap y seed is fixed in t1. At every pose, place pWrap z on the
% p1-pEnd line projected into the t1 y-z plane. Its x coordinate is then
% solved independently so pWrap lies on the nominal bpaRb standoff surface.
y1 = p1T1(:,2);
z1 = p1T1(:,3);
yEnd = pEnd(2);
zEnd = pEnd(3);
denYZ = yEnd-y1;

if any(abs(denYZ) <= 1e-12)
    error('buildKneeFlexorRoute20mm:DegenerateYZLine', ...
        'p1 and pEnd have the same t1 y coordinate at one or more poses.')
end

pWrapY = ctx.wrapPointT1XY(2);

wrapYZFraction = (pWrapY - y1)./(yEnd - y1);
pWrapZ = z1 + wrapYZFraction.*(zEnd - z1);

% bpaRb is the BIG/nominal BPA radius. It is used only to locate the
% equivalent pWrap point relative to the hard tibial geometry. It is NOT the
% minimum collision clearance for the straight p1-pWrap-pEnd segments;
% nonlconExclusion uses the smaller bpaRs for that purpose.
targetClearance = ctx.geo.bpaRb;

% Find the closest point on this fixed y-z line to the hard tibial geometry.
% If the line reaches the bpaRb standoff surface, fzero locates the outside
% intersection used as pWrap. No tangency is solved here.
xSearchLo = min(ctx.geo.xBack, ctx.geo.circleCenter(1)) ...
          - 4*targetClearance;

xSearchHi = max(ctx.geo.xAfter, ctx.geo.circleCenter(1)) ...
          + 4*targetClearance;

pWrapX = zeros(N,1);
xClosest = zeros(N,1);
minimumDistance = zeros(N,1);
contactConstraintFrame = zeros(N,1);

for i = 1:N
    distanceAtX = @(x) signedDistanceToTibia20mm( ...
        [x, pWrapY, pWrapZ(i)], ctx.geo);

    clearanceResidual = @(x) ...
        distanceAtX(x) - targetClearance;

    [xClosest(i), minimumDistance(i)] = fminbnd( ...
        distanceAtX, xSearchLo, xSearchHi);

    % <= 0 means this y-z line reaches or passes inside the nominal bpaRb
    % standoff surface, so a nominal pWrap intersection can be found.
    contactConstraintFrame(i) = ...
        minimumDistance(i) - targetClearance;

    if contactConstraintFrame(i) <= 0
        xOutside = xSearchLo;
        searchStep = max(4*targetClearance, 0.05);

        for k = 1:20
            if clearanceResidual(xOutside) > 0
                break
            end

            xOutside = xOutside - searchStep;
            searchStep = 2*searchStep;
        end

        if clearanceResidual(xOutside) > 0
            pWrapX(i) = fzero( ...
                clearanceResidual, [xOutside, xClosest(i)]);
        else
            pWrapX(i) = xClosest(i);
            contactConstraintFrame(i) = 1;
        end
    else
        % No nominal bpaRb intersection exists on this prescribed y-z line.
        % The active-wrap constraint will reject this pose if contact is
        % still required by the current release logic.
        pWrapX(i) = xClosest(i);
    end
end

pWrap = [ ...
    pWrapX, ...
    repmat(pWrapY,N,1), ...
    pWrapZ];

for i = 1:N

    % Full 3-D vectors are used for the bend angle. The x-y projections below
    % are used only by the existing wrap-release test.
    vIn = p1T1(i,:)-pWrap(i,:);
    vOut = pWrap(i,:)-pEnd;
    vWrap = vOut;
    vDirect = p1T1(i,:)-pEnd;

    if norm(vIn) <= eps || norm(vOut) <= eps
        error('buildKneeFlexorRoute20mm:DegenerateRoute', ...
            'A flexor route segment has zero length at sample %d.', i)
    end

    % Shift the release-test vectors by -90 degrees in the t1 x-y plane.
    R = [0 1; -1 0]; %-90 degree rotation matrix
    vWrap9 = R*vWrap(1:2).';
    vDir9 = R*vDirect(1:2).';

    % 3-D angle between the two straight routing segments, radians.
    turnWrapped(i) = acos(dot(vIn,vOut)/(norm(vIn)*norm(vOut)));

    % User-specified release comparison, with both lines originating at
    % pEnd. Release pWrap when aDirectD is greater than aWrapEndD.
    aDirectD(i) = atan2d(vDir9(2), vDir9(1));
    aWrapEndD(i) = atan2d(vWrap9(2), vWrap9(1));
end

% Sweep from full extension toward flexion. The first frame satisfying the
% direct atan2 comparison releases pWrap.
[~, order] = sort(ctx.phiD, 'descend');
aDirectSweepD = aDirectD(order);
aWrapEndSweepD = aWrapEndD(order);

% Test from the first sample, including full extension.
releasePosition = find( ...
    aDirectSweepD > aWrapEndSweepD + ctx.wrapAngleToleranceD, ...
    1, 'first');

active = true(1,N);
releaseIndex = [];
releaseAngleD = NaN;

if ~isempty(releasePosition)
    active(order(releasePosition:end)) = false;

    % An initially inactive wrap is not a later release event.
    if releasePosition > 1
        releaseIndex = order(releasePosition);
        releaseAngleD = ctx.phiD(releaseIndex);
    end
end

% A geometrically available wrap is used only when the direct physical
% BPA/tendon route would otherwise violate the tibial exclusion envelope.
active = active & directNeedsWrap.';

% Require a nominal bpaRb contact location only where the current release
% logic says the wrap is active. This is a contact-location requirement, not
% the bpaRs straight-segment collision constraint.
if any(active)
    wrapContactConstraint = ...
        max(contactConstraintFrame(active));
else
    wrapContactConstraint = -targetClearance;
end

Location = zeros(3,3,N);
bendMeasure = zeros(N,1);

for i = 1:N
    pWrapICR = RowVecTrans(ctx.T_ICR_t1(:,:,i), pWrap(i,:));
    pEndICR = RowVecTrans(ctx.T_ICR_t1(:,:,i), pEnd);

    if active(i)
        Location(:,:,i) = [p1; pWrapICR; pEndICR];

        % wRap is independent of the collision radii. Xi3 uses the modeled
        % bent length wRap*angle to estimate bend-related BPA length loss.
        bendMeasure(i) = ctx.geo.wRap*abs(turnWrapped(i));
    else
        Location(:,:,i) = [p1; pEndICR; pEndICR];
        bendMeasure(i) = 0;
    end
end

info = struct;

[~, idxFullExtension] = max(ctx.phiD);
info.directNeedsWrap = directNeedsWrap;
info.directMinClearance = directMinClearance;
% Preserve the scalar field consumed by existing printed diagnostics.
info.wrapYZFraction = wrapYZFraction(idxFullExtension);
info.wrapYZFractionByFrame = wrapYZFraction;
info.pWrapT1 = pWrap;
info.p1T1 = p1T1;
info.active = active;
info.releaseIndex = releaseIndex;
info.releaseAngleD = releaseAngleD;
info.releaseFound = ~isempty(releaseIndex);
info.activeAtExtension = active(order(1));
% info.aInD = aInD;
% info.aOutD = aOutD;
info.turnWrappedD = rad2deg(turnWrapped);
info.aDirectD = aDirectD;
info.aWrapEndD = aWrapEndD;
info.aDirectSweepD = aDirectSweepD;
info.aWrapEndSweepD = aWrapEndSweepD;
info.sweepOrder = order;
info.finalSignedTurnD = aDirectSweepD(end) - aWrapEndSweepD(end);
info.wrapClearanceTarget = targetClearance;
info.wrapMinimumDistance = minimumDistance;
info.wrapDistanceAtPoint = signedDistanceToTibia20mm(pWrap, ctx.geo);
info.contactConstraintFrame = contactConstraintFrame;
info.wrapContactConstraint = wrapContactConstraint;
info.wrapContactFound = wrapContactConstraint <= 1e-10;
info.bpaRb = ctx.geo.bpaRb;
info.bpaRs = ctx.geo.bpaRs;
info.wRap = ctx.geo.wRap;

end
