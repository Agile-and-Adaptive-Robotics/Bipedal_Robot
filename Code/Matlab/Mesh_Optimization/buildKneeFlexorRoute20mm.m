function [Location, bendMeasure, info] = buildKneeFlexorRoute20mm(p1, pEnd, ctx)
%BUILDKNEEFLEXORROUTE20MM Build the flexor route with one t1 wrap point.
%
% p1 is fixed in the femur frame.  pEnd and the intermediate wrap point
% are fixed in the t1 frame.  The wrap is active at full extension and is
% removed at the first extension-to-flexion sample for which aDirect is
% greater than aWrap.  An inactive wrap is represented by repeating pEnd.

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

% First express p1 in t1 at every knee angle.  The combined femur-to-t1
% transform is precomputed once by buildKneeFlexorContext20mm.
for i = 1:N
    p1T1(i,:) = RowVecTrans(ctx.T_t1_f(:,:,i), p1);
end

% The wrap x-y seed is fixed in t1.  At full extension, place its z
% coordinate on the p1-pEnd line projected into the t1 y-z plane.
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

% pWrap represents intentional BPA contact with the tibia. Therefore its
% centerline distance from the tibia equals the inflated BPA radius.
targetClearance = ctx.geo.bpaRadius;

% Find the closest point on this fixed y-z line to the tibial geometry.
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
        pWrapX(i) = xClosest(i);
    end
end

pWrap = [ ...
    pWrapX, ...
    repmat(pWrapY,N,1), ...
    pWrapZ];

for i = 1:N

    vIn = p1T1(i,:)-pWrap(i,:);
    vOut = pEnd-pWrap(i,:);
    vWrap = -vOut;
    vDirect = p1T1(i,:)-pEnd;

    if norm(vIn) <= eps || norm(vOut) <= eps
        error('buildKneeFlexorRoute20mm:DegenerateRoute', ...
            'A flexor route segment has zero length at sample %d.', i)
    end

    %shift vector by -90 degrees
    R = [0 1; -1 0]; %-90 degree rotation matrix
    vWrap9 = R*vWrap(1:2).';
    vDir9 = R*vDirect(1:2).';

    % Angle between the vectors
    turnWrapped(i) = acos(dot(vIn,vOut)/(norm(vIn)*norm(vOut)));

    % User-specified release comparison, with both lines originating at
    % pEnd.  Release pWrap when aDirectD is greater than aWrapEndD.
    aDirectD(i) = atan2d( vDir9(2), vDir9(1));
    aWrapEndD(i) = atan2d( vWrap9(2), vWrap9(1));
end

% Sweep from full extension toward flexion.  The first frame satisfying the
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

% Require a valid contact location only where contact is active.
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
        bendMeasure(i) = ctx.wrapBendRadius*abs(turnWrapped(i));
    else
        Location(:,:,i) = [p1; pEndICR; pEndICR];
        bendMeasure(i) = 0;
    end
end

info = struct;

[~, idxFullExtension] = max(ctx.phiD);

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

end
