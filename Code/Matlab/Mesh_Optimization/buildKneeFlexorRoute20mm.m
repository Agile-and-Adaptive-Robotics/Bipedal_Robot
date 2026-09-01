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
aInD = zeros(N,1);
aOutD = zeros(N,1);
turnWrappedD = zeros(N,1);
aDirectD = zeros(N,1);
aWrapEndD = zeros(N,1);

% First express p1 in t1 at every knee angle.  The combined femur-to-t1
% transform is precomputed once by buildKneeFlexorContext20mm.
for i = 1:N
    p1T1(i,:) = RowVecTrans(ctx.T_t1_f(:,:,i), p1);
end

% The wrap x-y seed is fixed in t1.  At full extension, place its z
% coordinate on the p1-pEnd line projected into the t1 y-z plane.
[~, idxExtension] = max(ctx.phiD);
y1 = p1T1(idxExtension,2);
z1 = p1T1(idxExtension,3);
yEnd = pEnd(2);
zEnd = pEnd(3);

if abs(yEnd - y1) <= eps
    error('buildKneeFlexorRoute20mm:DegenerateYZLine', ...
        'p1 and pEnd have the same t1 y coordinate at full extension.')
end

pWrapY = ctx.wrapPointT1XY(2);

wrapYZFraction = (pWrapY - y1)/(yEnd - y1);
pWrapZ = z1 + wrapYZFraction*(zEnd - z1);

% pWrap represents intentional BPA contact with the tibia. Therefore its
% centerline distance from the tibia equals the inflated BPA radius.
targetClearance = ctx.geo.bpaRadius;

distanceAtX = @(x) signedDistanceToTibia20mm( ...
    [x, pWrapY, pWrapZ], ctx.geo);

clearanceResidual = @(x) distanceAtX(x) - targetClearance;

% Find the closest point on this fixed y-z line to the tibial geometry.
xSearchLo = min(ctx.geo.xBack, ctx.geo.circleCenter(1)) ...
          - 4*targetClearance;

xSearchHi = max(ctx.geo.xAfter, ctx.geo.circleCenter(1)) ...
          + 4*targetClearance;

[xClosest, minimumDistance] = fminbnd( ...
    distanceAtX, xSearchLo, xSearchHi);

% This is a nonlinear feasibility constraint:
%   <= 0: the fixed y-z line reaches the clearance surface
%    > 0: the entire line is too far away to establish wrap contact
wrapContactConstraint = minimumDistance - targetClearance;

if wrapContactConstraint <= 1e-10
    % Find the negative-x clearance root.
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
        pWrapX = fzero(clearanceResidual, [xOutside, xClosest]);
    else
        % Do not crash an optimizer evaluation.
        pWrapX = xClosest;
        wrapContactConstraint = 1;
    end
else
    % No exact contact point exists for this candidate y-z line. Use its
    % closest point so the model remains finite, then reject the candidate
    % through wrapContactConstraint.
    pWrapX = xClosest;
end

pWrap = [pWrapX, pWrapY, pWrapZ];

for i = 1:N

    vIn = pWrap(1:2) - p1T1(i,1:2);
    vOut = pEnd(1:2) - pWrap(1:2);

    if norm(vIn) <= eps || norm(vOut) <= eps
        error('buildKneeFlexorRoute20mm:DegenerateRoute', ...
            'A flexor route segment has zero length at sample %d.', i)
    end

    % Preserve the current clockwise-positive convention for the bend
    % angle.  These values affect bendMeasure, not wrap activation.
    aInD(i) = atan2d(-vIn(2), vIn(1));
    aOutD(i) = atan2d(-vOut(2), vOut(1));

    % Signed difference between the two line angles.  No dot product,
    % cross product, or separate vector-angle helper is used.
    turnWrappedD(i) = atan2d( ...
        sind(aOutD(i) - aInD(i)), ...
        cosd(aOutD(i) - aInD(i)));

    % User-specified release comparison, with both lines originating at
    % pEnd.  Release pWrap when aDirectD is greater than aWrapEndD.
    aDirectD(i) = atan2d( ...
        p1T1(i,2) - pEnd(2), p1T1(i,1) - pEnd(1));
    aWrapEndD(i) = atan2d( ...
        pWrap(2) - pEnd(2), pWrap(1) - pEnd(1));
end

% Sweep from full extension toward flexion.  The first frame satisfying the
% direct atan2 comparison releases pWrap.
[~, order] = sort(ctx.phiD, 'descend');
aDirectSweepD = aDirectD(order);
aWrapEndSweepD = aWrapEndD(order);

% The fully extended sample is always active.  Begin applying the release
% inequality at the next sample in the extension-to-flexion sweep.
releaseOffset = find( ...
    aDirectSweepD(2:end) > ...
        aWrapEndSweepD(2:end) + ctx.wrapAngleToleranceD, ...
    1, 'first');

if isempty(releaseOffset)
    releasePosition = [];
else
    releasePosition = releaseOffset + 1;
end

active = true(1,N);
releaseIndex = [];
releaseAngleD = NaN;

if ~isempty(releasePosition)
    releaseIndex = order(releasePosition);
    active(order(releasePosition:end)) = false;

    releaseAngleD = ctx.phiD(releaseIndex);
end

Location = zeros(3,3,N);
bendMeasure = zeros(N,1);

for i = 1:N
    pWrapICR = RowVecTrans(ctx.T_ICR_t1(:,:,i), pWrap);
    pEndICR = RowVecTrans(ctx.T_ICR_t1(:,:,i), pEnd);

    if active(i)
        Location(:,:,i) = [p1; pWrapICR; pEndICR];
        bendMeasure(i) = ctx.wrapBendRadius*abs(deg2rad(turnWrappedD(i)));
    else
        Location(:,:,i) = [p1; pEndICR; pEndICR];
    end
end

info = struct;
info.pWrapT1 = pWrap;
info.wrapYZFraction = wrapYZFraction;
info.p1T1 = p1T1;
info.active = active;
info.releaseIndex = releaseIndex;
info.releaseAngleD = releaseAngleD;
info.releaseFound = ~isempty(releasePosition);
info.aInD = aInD;
info.aOutD = aOutD;
info.turnWrappedD = turnWrappedD;
info.aDirectD = aDirectD;
info.aWrapEndD = aWrapEndD;
info.aDirectSweepD = aDirectSweepD;
info.aWrapEndSweepD = aWrapEndSweepD;
info.sweepOrder = order;
info.finalSignedTurnD = aDirectSweepD(end) - aWrapEndSweepD(end);
info.wrapClearanceTarget = targetClearance;
info.wrapMinimumDistance = minimumDistance;
info.wrapDistanceAtPoint = distanceAtX(pWrapX);
info.wrapContactConstraint = wrapContactConstraint;
info.wrapContactFound = wrapContactConstraint <= 1e-10;

end
