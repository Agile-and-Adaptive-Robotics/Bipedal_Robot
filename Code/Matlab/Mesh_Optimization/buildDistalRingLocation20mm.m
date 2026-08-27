function [Location, bendMeasure, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx)
% buildDistalRingLocation20mm
%
% 20 mm extensor route built by initializing the fully-flexed route and
% sweeping toward extension while only eliminating contacts.
%
% Fixed-size route representation:
%   p1:p5  are femur-frame points.
%   p6:p9  are t1-frame points and are transformed to ICR before Location output.
%   CrossPoint = 6.
%
% Optimization vector:
%   x = [p1(1:3), pEnd(1:3), rest, tendon]
%
% p9 is the distal t1 attachment pEnd.
%
% Eliminated femur points repeat the previous unique femur point.
% Eliminated t1 points repeat the next unique t1 point.
%
% p1 and p9 are optimizer-supplied. p2:p8 come from the SolidWorks seed,
% except for local contact corrections needed after p1/p9 movement.
% p2 and p8 remain as the local bend regions while the broader p3-to-
% downstream condyle curvature model is deferred.

N = ctx.N;
phiD = ctx.phiD(:);

if N < 2
    error('buildDistalRingLocation20mm:NeedAtLeastTwoAngles', ...
        'ctx.phiD must contain at least two knee positions.')
end

dPhiD = abs(phiD(2) - phiD(1));
if ~isfinite(dPhiD) || dPhiD <= 0
    error('buildDistalRingLocation20mm:BadAngleSpacing', ...
        'ctx.phiD must have a finite nonzero increment.')
end

%% Flexion-to-extension sweep over the solver grid
sweepD = phiD;
Ns = numel(sweepD);

%% SolidWorks seed / branch reference
% The seed is defined centrally in buildKneeExtContext20mm.
seed = ctx.routeSeed;

if ~isequal(size(seed), [9 3])
    error('buildDistalRingLocation20mm:BadSeedSize', ...
        'ctx.routeSeed must be 9x3.')
end

% p1 and pEnd remain optimization variables.
seed(1,:) = p1;
seed(9,:) = pEnd;

% Preserve the SolidWorks x-y seed coordinates. Only the native-frame z
% coordinate follows the optimized endpoint on each bone.
seed20 = adjustedSeed20mm(seed, p1, pEnd, ctx);

%% Output allocation
Location = zeros(9,3,N);
bendRadius = zeros(9,N);

info.raw = zeros(9,3,N);
info.candidateRaw = zeros(9,3,N);
info.active = false(9,N);
info.tendonMax = zeros(N,1);
info.endcap_t1 = zeros(N,3);
info.targetType = zeros(N,1);
info.femurOverflow = false(N,1);
info.femurCrossingCollision = false(N,1);
info.angleCull.pointAngleD = nan(9,N);
info.angleCull.bypassAngleD = nan(9,N);
info.angleCull.marginD = nan(9,N);
info.angleCull.tested = false(9,N);
info.angleCull.angleGate = false(9,N);
info.angleCull.bypassClear = false(9,N);
info.angleCull.removed = false(9,N);
info.angleCull.candidateActive = false(9,N);
info.angleCull.retainedActive = false(9,N);

info.CrossPoint = 6;
info.routeRows = 9;
info.pointFrame = [ ...
    "femur";"femur";"femur";"femur";"femur"; ...
    "t1";"t1";"t1";"t1"];

info.addedAngleD = nan(9,1);
info.addedPoint = nan(9,3);
info.addedSweepIndex = nan(9,1);
info.eliminatedAngleD = nan(9,1);
info.eliminatedPoint = nan(9,3);
info.eliminatedSweepIndex = nan(9,1);

%% Route contact state
% The SolidWorks seed gives the fully-flexed candidate wrap chain. The route
% starts with p1:p9 unique at -120 deg, then carries state toward extension.
% Only optional p3:p7 contacts can be eliminated. p2 and p8 stay available
% as local bend/contact regions.

% Candidate point positions are fixed here in their native frames.
Pfixed = seed20;
Pfixed(1,:) = p1;
Pfixed(9,:) = pEnd;

activeState = true(9,1);

info.addedAngleD(:) = sweepD(1);
info.addedPoint(:,:) = Pfixed;
info.addedSweepIndex(:) = 1;

contactTol = 1e-8;
angleTolD = 0.05;

for s = 1:Ns

    thetaD = sweepD(s);
    solverIndex = s;

    T_Pam_i = ctx.T_Pam(:,:,solverIndex);
    T_Pam_inv_i = ctx.T_Pam_inv(:,:,solverIndex);
    T_ICR_t1_i = ctx.T_ICR_t1(:,:,solverIndex);
    T_t1_ICR_i = ctx.T_t1_ICR(:,:,solverIndex);

    [~, tendonMaxGeom, targetRadius, targetType, ~] = ...
        tendonLimit20mm(pEnd, ctx.geo);

    frameInfo.tendonMax = tendonMaxGeom;
    frameInfo.targetType = targetType;
    frameInfo.targetRadius = targetRadius;

    % Physical BPA endcap point along pEnd -> tangent direction.
    % Diagnostic only; it is not a route row.
    [pTanGeom, ~, ~, ~, ~] = tendonLimit20mm(pEnd, ctx.geo);

    if tendonMaxGeom > 1e-12
        frac = min(max(tendon/tendonMaxGeom,0),1);
        endcapXY = pEnd(1:2) + frac*(pTanGeom - pEnd(1:2));
    else
        endcapXY = pEnd(1:2);
    end

    frameInfo.endcap_t1 = [endcapXY, pEnd(3)];

    angleFrame = emptyAngleFrame();

    if s == 1
        removedNow = false(9,1);
    else
        [activeState, angleFrame, removedFemur] = eliminateFemurContacts( ...
            Pfixed, activeState, T_Pam_i, T_ICR_t1_i, ...
            ctx.geo, contactTol, angleTolD, angleFrame);

        [activeState, angleFrame, removedT1] = eliminateT1Contacts( ...
            Pfixed, activeState, T_Pam_inv_i, T_t1_ICR_i, ...
            ctx.geo, contactTol, angleTolD, angleFrame);

        removedNow = removedFemur | removedT1;
    end

    active = activeState;

    fixedBendRadius = routeBendRadii(active, seed20, ctx.geo);

    for jj = find(removedNow).'
        if isnan(info.eliminatedAngleD(jj))
            info.eliminatedAngleD(jj) = thetaD;
            info.eliminatedPoint(jj,:) = Pfixed(jj,:);
            info.eliminatedSweepIndex(jj) = s;
        end
    end

    raw = collapseInactive(Pfixed, active);

    % Interference checks are disabled while validating the absolute-angle
    % route state. The stored diagnostics stay false so downstream plots do
    % not mix clearance failures with topology transitions.
    % B_t1 = raw(6,:);
    % B_icr = RowVecTrans(T_ICR_t1_i, B_t1);
    % B_fem = RowVecTrans(T_Pam_i, B_icr);
    % crossingCollision = segmentPenetratesFemurOffset( ...
    %     raw(5,1:2), B_fem(1:2), ctx.geo, contactTol);
    % overflow = active(5) && crossingCollision;
    crossingCollision = false;
    overflow = false;

    rBend = fixedBendRadius;
    rBend(~active) = 0;

    % Repeated points are not physical bends.
    for j = 2:8
        if norm(raw(j,:) - raw(j-1,:)) < 1e-10 || ...
                norm(raw(j+1,:) - raw(j,:)) < 1e-10
            rBend(j) = 0;
        end
    end

    info.raw(:,:,solverIndex) = raw;
    info.candidateRaw(:,:,solverIndex) = Pfixed;
    info.active(:,solverIndex) = active;
    info.tendonMax(solverIndex) = frameInfo.tendonMax;
    info.endcap_t1(solverIndex,:) = frameInfo.endcap_t1;
    info.targetType(solverIndex) = frameInfo.targetType;
    info.femurOverflow(solverIndex) = overflow;
    info.femurCrossingCollision(solverIndex) = crossingCollision;
    info.angleCull.pointAngleD(:,solverIndex) = angleFrame.pointAngleD;
    info.angleCull.bypassAngleD(:,solverIndex) = angleFrame.bypassAngleD;
    info.angleCull.marginD(:,solverIndex) = angleFrame.marginD;
    info.angleCull.tested(:,solverIndex) = angleFrame.tested;
    info.angleCull.angleGate(:,solverIndex) = angleFrame.angleGate;
    info.angleCull.bypassClear(:,solverIndex) = angleFrame.bypassClear;
    info.angleCull.removed(:,solverIndex) = angleFrame.removed;
    info.angleCull.candidateActive(:,solverIndex) = active;
    info.angleCull.retainedActive(:,solverIndex) = active;
    bendRadius(:,solverIndex) = rBend;

    % p6:p9 are native t1 points. Convert them t1 -> ICR exactly as in
    % the existing model before MonoPamDataExplicit_balanceX3 applies
    % T_Pam at CrossPoint = 6.
    loc = raw;

    for j = 6:9
        loc(j,:) = RowVecTrans(T_ICR_t1_i, loc(j,:));
    end

    Location(:,:,solverIndex) = loc;
end

%% Xi3 geometric bend measure
[bendMeasure, wrapInfo] = geometricBendMeasure( ...
    Location, ctx.T_Pam, ctx.T_Pam_inv, ctx.T_t1_ICR, ctx.phiD, ...
    ctx.geo.alphaTol, ctx.geo, info.active, info.raw);

info.bendRadius = bendRadius;
info.wrap = wrapInfo;
info.wrapAngleRad = wrapInfo.angleRad;
info.wrapLength = wrapInfo.length;
info.sweepD = sweepD;
info.sweepStepD = dPhiD;
info.seedOriginal = ctx.routeSeed;
info.seed20 = seed20;

idx = find(isfinite(info.addedAngleD));

if isempty(idx)
    info.activationOrder = zeros(0,1);
else
    keyAngle = info.addedAngleD(idx);
    [~,ord] = sortrows([-keyAngle, idx],[1 2]);
    info.activationOrder = idx(ord);
end

idx = find(isfinite(info.eliminatedAngleD));

if isempty(idx)
    info.eliminationOrder = zeros(0,1);
else
    keyAngle = info.eliminatedAngleD(idx);
    tieRank = arrayfun(@eliminationTieRank, idx);
    [~,ord] = sortrows([keyAngle, tieRank(:)],[1 2]);
    info.eliminationOrder = idx(ord);
end

end


function rank = eliminationTieRank(row)

switch row
    case 5
        rank = 1;
    case 4
        rank = 2;
    case 3
        rank = 3;
    case 6
        rank = 4;
    case 7
        rank = 5;
    otherwise
        rank = row + 10;
end

end


%% =====================================================================
%% Seed geometry
%% =====================================================================

function seed20 = adjustedSeed20mm(seed, p1, pEnd, ctx)
% Preserve the SolidWorks route as the branch reference while allowing the
% local contacts that should respond to endpoint motion to move from their
% seed locations.
%
% p2: femur-cylinder contact.  This should only move slightly as p1 moves.
% p8: lower-tibia-cylinder contact.  This can move substantially when pEnd
%     moves, which is why it must not remain frozen at the original seed.
%
% p3:p7 remain the SolidWorks seed/reference coordinates and are retained or
% bypassed by the flexion-to-extension elimination logic.

geo = ctx.geo;

seed20 = seed;
seed20(1,:) = p1;
seed20(9,:) = pEnd;

seed20(2:5,3) = p1(3);
seed20(6:8,3) = pEnd(3);

[q2, ok2] = tangentPointCircleNearest( ...
    p1(1:2), geo.femurCylCenter, geo.femurCylClearRadius, seed(2,1:2));

if ok2
    seed20(2,1:2) = q2;
end

% This clearance-driven seed repair is disabled while validating the route
% state machine. It changes p3:p5 before any elimination tests run, which
% makes it harder to separate topology errors from seed-geometry changes.
% iFlex = ctx.idxMaxFlex;
% p6FlexFemur = t1PointToFemur( ...
%     seed20(6,:), ...
%     ctx.T_Pam(:,:,iFlex), ...
%     ctx.T_ICR_t1(:,:,iFlex));
%
% [qFemurChain, okFemurChain] = adjustedFemurSeedChainAtFlexion( ...
%     seed20(2,1:2), ...
%     seed20(3:5,1:2), ...
%     p6FlexFemur(1:2), ...
%     geo, ...
%     1e-8);
%
% if okFemurChain
%     seed20(3:5,1:2) = qFemurChain;
% end

[q8, ok8] = tangentPointCircleNearest( ...
    pEnd(1:2), geo.tibiaLowerCenter, geo.tibiaLowerClearRadius, seed(8,1:2));

if ok8
    seed20(8,1:2) = q8;
end

end


%% =====================================================================
%% Absolute-angle contact elimination
%% =====================================================================

function angleFrame = emptyAngleFrame()
angleFrame.pointAngleD = nan(9,1);
angleFrame.bypassAngleD = nan(9,1);
angleFrame.marginD = nan(9,1);
angleFrame.tested = false(9,1);
angleFrame.angleGate = false(9,1);
angleFrame.bypassClear = false(9,1);
angleFrame.removed = false(9,1);

end


function [activeState, angleFrame, removed] = eliminateFemurContacts( ...
    Pfixed, activeState, T_Pam_i, T_ICR_t1_i, geo, tol, angleTolD, angleFrame)
% Remove the last unique femur-side optional contact if the absolute-angle
% test passes. The bypass collision gate is retained below but disabled
% while the angle-only route state is being validated.

removed = false(9,1);

while true

    kOpt = find(activeState(3:5), 1, 'last');

    if isempty(kOpt)
        return
    end

    k = kOpt + 2;
    iPrev = find(activeState(1:k-1), 1, 'last');
    iNext = find(activeState(6:9), 1, 'first') + 5;

    if isempty(iPrev) || isempty(iNext)
        return
    end

    raw = collapseInactive(Pfixed, activeState);
    pNextFemur = t1PointToFemur(raw(iNext,:), T_Pam_i, T_ICR_t1_i);

    [angleGate, aBigD, aSmallD, marginD] = femurRemovalAngle( ...
        raw(iPrev,:), raw(k,:), pNextFemur, angleTolD);

    % The clearance polygon/tolerance gate can reject the intended topology
    % change before the absolute-angle state machine can be evaluated.
    % Re-enable this after the route order and p2/p8 wrap model are stable.
    % bypassClear = femurBypassCollisionFree( ...
    %     raw(iPrev,:), pNextFemur, geo, tol, k);
    bypassClear = true;

    angleFrame = recordEliminationTest( ...
        angleFrame, k, aSmallD, aBigD, marginD, angleGate, bypassClear);

    if angleGate
        activeState(k) = false;
        removed(k) = true;
        angleFrame.removed(k) = true;
    else
        return
    end
end

end


function [activeState, angleFrame, removed] = eliminateT1Contacts( ...
    Pfixed, activeState, T_Pam_inv_i, T_t1_ICR_i, geo, tol, angleTolD, angleFrame)
% Remove the first unique t1-side optional contact if the absolute-angle
% test passes. The bypass collision gate is retained below but disabled
% while the angle-only route state is being validated.

removed = false(9,1);

while true

    jOpt = find(activeState(6:7), 1, 'first');

    if isempty(jOpt)
        return
    end

    j = jOpt + 5;
    iPrev = find(activeState(1:5), 1, 'last');
    iNext = find(activeState(j+1:9), 1, 'first') + j;

    if isempty(iPrev) || isempty(iNext)
        return
    end

    raw = collapseInactive(Pfixed, activeState);
    pPrevT1 = femurPointToT1(raw(iPrev,:), T_Pam_inv_i, T_t1_ICR_i);

    [angleGate, aBigD, aSmallD, marginD] = t1RemovalAngle( ...
        pPrevT1, raw(j,:), raw(iNext,:), angleTolD);

    % The clearance polygon/tolerance gate can reject the intended topology
    % change before the absolute-angle state machine can be evaluated.
    % Re-enable this after the route order and p2/p8 wrap model are stable.
    % bypassClear = t1BypassCollisionFree( ...
    %     pPrevT1, raw(iNext,:), geo, tol, j);
    bypassClear = true;

    angleFrame = recordEliminationTest( ...
        angleFrame, j, aSmallD, aBigD, marginD, angleGate, bypassClear);

    if angleGate
        activeState(j) = false;
        removed(j) = true;
        angleFrame.removed(j) = true;
    else
        return
    end
end

end


function [removePoint, aBigD, aSmallD, marginD] = femurRemovalAngle( ...
    pPrev, pSmallPoint, pNextFemur, angleTolD)

pBig = pNextFemur(1:2) - pPrev(1:2);
pSmall = pNextFemur(1:2) - pSmallPoint(1:2);

% Femur-side removal is a clockwise comparison in the plotted/anatomical
% view. Flipping y before atan2 keeps the user's aBig > aSmall rule while
% avoiding the full-flexed seed being treated as already removable.
% The principal-angle values are compared directly.
aBigD = atan2d(-pBig(2), pBig(1));
aSmallD = atan2d(-pSmall(2), pSmall(1));
marginD = aBigD - aSmallD;

removePoint = marginD > angleTolD;

end


function [removePoint, aBigD, aSmallD, marginD] = t1RemovalAngle( ...
    pPrevT1, pSmallPoint, pNext, angleTolD)

pBig = pPrevT1(1:2) - pNext(1:2);
pSmall = pSmallPoint(1:2) - pNext(1:2);

aBigD = atan2d(pBig(2), pBig(1));
aSmallD = atan2d(pSmall(2), pSmall(1));
marginD = aSmallD - aBigD;

removePoint = marginD > angleTolD;

end


function angleFrame = recordEliminationTest( ...
    angleFrame, row, aSmallD, aBigD, marginD, angleGate, bypassClear)

angleFrame.pointAngleD(row) = aSmallD;
angleFrame.bypassAngleD(row) = aBigD;
angleFrame.marginD(row) = marginD;
angleFrame.tested(row) = true;
angleFrame.angleGate(row) = angleGate;
angleFrame.bypassClear(row) = bypassClear;

end


function qFemur = t1PointToFemur(qT1, T_Pam_i, T_ICR_t1_i)

qFemur = RowVecTrans(T_Pam_i*T_ICR_t1_i, qT1);

end


function qT1 = femurPointToT1(qFemur, T_Pam_inv_i, T_t1_ICR_i)

qICR = RowVecTrans(T_Pam_inv_i, qFemur);
qT1 = RowVecTrans(T_t1_ICR_i, qICR);

end


function clear = femurBypassCollisionFree(A, B, geo, tol, candidateRow)

A = A(1:2);
B = B(1:2);

switch candidateRow
    case 3
        clear = ...
            ~segmentPenetratesCircle(A, B, ...
                geo.femurCylCenter, geo.femurCylClearRadius, tol) && ...
            ~segmentPenetratesFemurOffset(A, B, geo, tol) && ...
            ~segmentIntersectsVerticalSpan(A, B, ...
                geo.femurLineX, geo.femurLineY, tol);

    case {4, 5}
        clear = ~segmentPenetratesFemurOffset(A, B, geo, tol);

    otherwise
        clear = true;
end

end


function clear = t1BypassCollisionFree(A, B, geo, tol, candidateRow)

A = A(1:2);
B = B(1:2);

clear = true;

switch candidateRow
    case 6
        clear = ~segmentPenetratesCircle(A, B, ...
            geo.tibiaUpperCenter, geo.tibiaUpperClearRadius, tol);

    case 7
        clear = ~segmentIntersectsVerticalSpan(A, B, ...
            geo.tibiaWallX, geo.tibiaWallY, tol);
end

end


function rBend = routeBendRadii(active, seed20, geo)

rBend = zeros(9,1);

if active(2)
    rBend(2) = geo.femurCylClearRadius;
end

if active(8)
    rBend(8) = geo.tibiaLowerClearRadius;
end

end


%% =====================================================================
%% Repeated-point fixed-size route
%% =====================================================================

function raw = collapseInactive(Pfixed, active)

raw = Pfixed;

% Femur side: inactive rows repeat the preceding active femur point.
raw(1,:) = Pfixed(1,:);

for j = 2:5
    if active(j)
        raw(j,:) = Pfixed(j,:);
    else
        raw(j,:) = raw(j-1,:);
    end
end

% Tibia side: inactive rows repeat the following active t1 point.
raw(9,:) = Pfixed(9,:);

for j = 8:-1:6
    if active(j)
        raw(j,:) = Pfixed(j,:);
    else
        raw(j,:) = raw(j+1,:);
    end
end

end


function [P, ok] = adjustedFemurSeedChainAtFlexion(A, Pseed, B, geo, tol)
% Move the SolidWorks femur-chain seed just far enough outside the corrected
% clearance envelope for the fully-flexed route to be collision-free.

P = Pseed;
ok = false;

Pbase = Pseed;

for k = 1:size(Pbase,1)
    [q, okNearest] = nearestPointOnPolyline(Pbase(k,:), geo.femurOffsetBoundary);

    if okNearest
        Pbase(k,:) = q;
    end
end

dir = Pbase - geo.femurProfileCenter;
dirNorm = vecnorm(dir, 2, 2);

bad = dirNorm < 1e-12;
if any(bad)
    dir(bad,:) = Pseed(bad,:) - geo.femurProfileCenter;
    dirNorm(bad) = vecnorm(dir(bad,:), 2, 2);
end

bad = dirNorm < 1e-12;
dir(bad,:) = repmat([1 0], nnz(bad), 1);
dirNorm(bad) = 1;

dir = dir ./ dirNorm;

for clearance = 0:0.00025:0.050
    Pcand = Pbase + clearance*dir;

    if femurSeedChainClear(A, Pcand, B, geo, tol)
        P = Pcand;
        ok = true;
        return
    end
end

end


function clear = femurSeedChainClear(A, P, B, geo, tol)

chain = [A; P; B];
clear = true;

for k = 1:size(chain,1)-1
    if segmentPenetratesFemurOffset(chain(k,:), chain(k+1,:), geo, tol)
        clear = false;
        return
    end
end

end


%% =====================================================================
%% Seeded contact helpers
%% =====================================================================

function [q, rhoEff, ok] = seededCondyleContactCandidate(A,B,geo,qRef,tol)
% Choose the tangent point from the DOWNSTREAM endpoint B to the corrected
% femoral clearance boundary that belongs to the SolidWorks seed branch.
%
% The candidate is accepted only if:
%   1. A -> q remains outside the corrected condyle clearance; and
%   2. q -> B remains outside the corrected condyle clearance.
%
% Because q is a true tangent point from B, q -> B is tangent by
% construction. The first test prevents us from choosing the wrong side of
% the ellipse simply because it is mathematically valid.

[Q, thetaRoots] = tangentPointsToFemurOffset(B,geo);

q = [NaN NaN];
rhoEff = 0;
ok = false;

if isempty(thetaRoots)
    return
end

bestScore = inf;

for k = 1:size(Q,1)

    qk = Q(k,:);

    if segmentPenetratesFemurOffset(A,qk,geo,tol)
        continue
    end

    if segmentPenetratesFemurOffset(qk,B,geo,tol)
        continue
    end

    % Stay on the SolidWorks branch. Distance to qRef dominates selection.
    score = norm(qk-qRef) + 0.05*(norm(qk-A)+norm(B-qk));

    if score < bestScore
        bestScore = score;
        q = qk;
        rhoEff = offsetEllipseRadiusOfCurvature(thetaRoots(k),geo);
        ok = true;
    end
end

end


function [q, ok] = seededCircleContactCandidate(A,B,C,R,qRef,tol)
% One seeded boundary contact for a circle. Try tangent points from the
% upstream endpoint A and retain the branch nearest the SolidWorks seed,
% provided the remaining q -> B segment does not enter the circle.

[Q, okTan] = tangentPointsCircle(A,C,R);

q = [NaN NaN];
ok = false;

if ~okTan
    return
end

bestScore = inf;

for k = 1:2

    qk = Q(k,:);

    if segmentPenetratesCircle(A,qk,C,R,tol)
        continue
    end

    if segmentPenetratesCircle(qk,B,C,R,tol)
        continue
    end

    score = norm(qk-qRef);

    if score < bestScore
        bestScore = score;
        q = qk;
        ok = true;
    end
end

end


function [q, ok] = tangentPointCircleNearest(P,C,R,qRef)
% External tangent point from P to circle C,R, choosing the branch nearest
% qRef. Used for the always-active p2 and p8 baseline contacts so those
% contacts respond to optimized p1/pEnd motion.

[Q, okTan] = tangentPointsCircle(P,C,R);

q = [NaN NaN];
ok = false;

if ~okTan
    return
end

[~,k] = min(vecnorm(Q-qRef,2,2));

q = Q(k,:);
ok = true;

end


function [q, ok] = nearestPointOnPolyline(P, poly)

q = [NaN NaN];
ok = false;

if isempty(poly)
    return
end

A = poly;
B = [poly(2:end,:); poly(1,:)];
D = B - A;
d2 = sum(D.^2, 2);
valid = d2 > 1e-18;

if ~any(valid)
    return
end

t = zeros(size(d2));
t(valid) = sum((P - A(valid,:)).*D(valid,:), 2)./d2(valid);
t = min(max(t, 0), 1);

Q = A + t.*D;
d = vecnorm(Q - P, 2, 2);
[~,idx] = min(d);

q = Q(idx,:);
ok = all(isfinite(q));

end


%% =====================================================================
%% Circle / wall contact helpers
%% =====================================================================

function hit = segmentPenetratesCircle(A,B,C,R,tol)

d = B-A;
d2 = dot(d,d);

if d2 < 1e-18
    hit = norm(A-C) < R-tol;
    return
end

t = dot(C-A,d)/d2;
t = min(max(t,0),1);

q = A + t*d;

hit = norm(q-C) < R-tol;

end


function [V, ok] = equivalentCircleAvoidanceVertex(A,B,C,R,qRef,tol)
% Represent the wrap around a circle by one equivalent corner formed by the
% intersection of an A-side tangent and a B-side tangent.

[QA, okA] = tangentPointsCircle(A,C,R);
[QB, okB] = tangentPointsCircle(B,C,R);

V = [NaN NaN];
ok = false;

if ~okA || ~okB
    return
end

bestScore = inf;

for ia = 1:2
    for ib = 1:2

        qA = QA(ia,:);
        qB = QB(ib,:);

        [v, okLine] = lineIntersection2D(A,qA,B,qB);

        if ~okLine || any(~isfinite(v))
            continue
        end

        if norm(v-A) > 1.0 || norm(v-B) > 1.0
            continue
        end

        tA = pointParameterOnLine(A,v,qA);
        tB = pointParameterOnLine(B,v,qB);

        if tA < -tol || tA > 1+tol || tB < -tol || tB > 1+tol
            continue
        end

        if segmentPenetratesCircle(A,v,C,R,tol) || ...
                segmentPenetratesCircle(v,B,C,R,tol)
            continue
        end

        score = norm(v-A) + norm(v-B) + 0.15*norm(v-qRef);

        if score < bestScore
            bestScore = score;
            V = v;
            ok = true;
        end
    end
end

end


function [Q, ok] = tangentPointsCircle(P,C,R)

v = P-C;
d2 = dot(v,v);

Q = nan(2,2);
ok = false;

if d2 <= R^2
    return
end

base = C + (R^2/d2)*v;
h = R*sqrt(d2-R^2)/d2;
vp = [-v(2), v(1)];

Q(1,:) = base + h*vp;
Q(2,:) = base - h*vp;

ok = true;

end


function [hit,q] = segmentCircleContactPoint(A,B,C,R,qRef,tol)
% Existing boundary-contact behavior retained for the tibia route.

d = B-A;
d2 = dot(d,d);

if d2 < 1e-18
    qLine = A;
else
    t = dot(C-A,d)/d2;
    t = min(max(t,0),1);
    qLine = A + t*d;
end

v = qLine-C;
dist = norm(v);

hit = dist <= R + tol;

if ~hit
    q = [NaN NaN];
    return
end

if dist > 1e-12

    q = C + R*(v/dist);

else

    vr = qRef-C;
    nr = norm(vr);

    if nr < 1e-12
        vr = [1 0];
        nr = 1;
    end

    q = C + R*(vr/nr);
end

end


function [hit,q] = segmentVerticalContactPoint(A,B,xValue,yRange,tol)

q = [NaN NaN];

dx = B(1)-A(1);

if abs(dx) < 1e-14
    hit = false;
    return
end

t = (xValue-A(1))/dx;

if t < -tol || t > 1+tol
    hit = false;
    return
end

y = A(2)+t*(B(2)-A(2));

if y < min(yRange)-tol || y > max(yRange)+tol
    hit = false;
    return
end

hit = true;
q = [xValue,y];

end


function hit = segmentIntersectsVerticalSpan(A,B,xValue,yRange,tol)

dx = B(1)-A(1);

if abs(dx) < 1e-14
    hit = false;
    return
end

t = (xValue-A(1))/dx;

if t <= tol || t >= 1-tol
    hit = false;
    return
end

y = A(2)+t*(B(2)-A(2));
hit = y >= min(yRange)-tol && y <= max(yRange)+tol;

end


%% =====================================================================
%% Corrected femoral-condyle clearance / tangent helpers
%% =====================================================================

function hit = segmentPenetratesFemurOffset(A,B,geo,tol)
% Strict collision test against the clipped normal-offset femur polygon.
% Tangency is allowed; activation occurs once the straight segment actually
% enters the clearance envelope.

poly = geo.femurOffsetBoundary;

[inA,onA] = inpolygon(A(1),A(2),poly(:,1),poly(:,2));
[inB,onB] = inpolygon(B(1),B(2),poly(:,1),poly(:,2));

if (inA && ~onA) || (inB && ~onB)
    hit = true;
    return
end

tCuts = segmentPolygonIntersectionParams(A,B,poly,tol);

if isempty(tCuts)
    hit = false;
    return
end

tCuts = sort(unique([0; tCuts(:); 1]));

hit = false;

for k = 1:numel(tCuts)-1

    ta = tCuts(k);
    tb = tCuts(k+1);

    if tb-ta < 1e-10
        continue
    end

    tm = 0.5*(ta+tb);
    Pm = A + tm*(B-A);

    [in,on] = inpolygon(Pm(1),Pm(2),poly(:,1),poly(:,2));

    if in && ~on
        hit = true;
        return
    end
end

end


function [V, rhoEff, ok] = equivalentCondyleAvoidanceVertex(A,B,geo,qRef,tol)
% Replace a many-point condyle wrap cluster by one equivalent fixed vertex.
%
% Tangents are calculated against the unclipped normal-offset curve:
%   physical ellipse point + BPA_radius * outward_unit_normal.
% Tangency points on the clipped-away side of the vertical wall are rejected.
%
% Candidate tangent-line pairs are intersected. The valid equivalent corner
% closest to the desired seed branch and with the shortest route is used.

[QA, thetaA] = tangentPointsToFemurOffset(A,geo);
[QB, thetaB] = tangentPointsToFemurOffset(B,geo);

V = [NaN NaN];
rhoEff = 0;
ok = false;

if isempty(thetaA) || isempty(thetaB)
    return
end

bestScore = inf;

for ia = 1:size(QA,1)
    for ib = 1:size(QB,1)

        qA = QA(ia,:);
        qB = QB(ib,:);

        [v, okLine] = lineIntersection2D(A,qA,B,qB);

        if ~okLine || any(~isfinite(v))
            continue
        end

        if norm(v-A) > 1.0 || norm(v-B) > 1.0
            continue
        end

        tA = pointParameterOnLine(A,v,qA);
        tB = pointParameterOnLine(B,v,qB);

        if tA < -tol || tA > 1+tol || tB < -tol || tB > 1+tol
            continue
        end

        % The two straight segments must remain outside the corrected
        % condyle clearance envelope.
        if segmentPenetratesFemurOffset(A,v,geo,tol) || ...
                segmentPenetratesFemurOffset(v,B,geo,tol)
            continue
        end

        routeLen = norm(v-A) + norm(v-B);

        % Branch selection is guided primarily by the expected location on
        % the condyle clearance curve, not by the globally shortest way
        % around the opposite side of the ellipse.
        branchPenalty = norm(qA-qRef) + norm(qB-qRef);

        score = routeLen + 6.0*branchPenalty;

        if score < bestScore

            bestScore = score;

            V = v;
            rhoA = offsetEllipseRadiusOfCurvature(thetaA(ia),geo);
            rhoB = offsetEllipseRadiusOfCurvature(thetaB(ib),geo);
            rhoEff = 0.5*(rhoA+rhoB);

            ok = true;
        end
    end
end

end


function [Q, thetaRoots] = tangentPointsToFemurOffset(P,geo)
% Find external tangencies from P to the unclipped normal-offset ellipse.
% Tangencies on the clipped-away side of the femur vertical wall are removed.
%
% For an offset point Q(theta), the tangent condition is
%   (P-Q) dot n_hat(theta) = 0.
%
% The offset curve has the same tangent direction / normal as the physical
% ellipse, so this condition is exact.

thetaGrid = geo.femurOffsetThetaSearch(:);

f = zeros(size(thetaGrid));

for k = 1:numel(thetaGrid)

    [Qk,nk] = femurOffsetPointNormal(thetaGrid(k),geo);

    f(k) = dot(P-Qk,nk);
end

rootsFound = [];

for k = 1:numel(thetaGrid)-1

    a = thetaGrid(k);
    b = thetaGrid(k+1);

    fa = f(k);
    fb = f(k+1);

    if abs(fa) < 1e-10
        rootsFound(end+1,1) = a; %#ok<AGROW>
    end

    if fa*fb < 0

        root = bisectTangentRoot(P,a,b,geo);
        rootsFound(end+1,1) = root; %#ok<AGROW>
    end
end

if abs(f(end)) < 1e-10
    rootsFound(end+1,1) = thetaGrid(end); %#ok<AGROW>
end

if isempty(rootsFound)

    Q = zeros(0,2);
    thetaRoots = zeros(0,1);
    return
end

% Normalize and remove duplicate roots at the +/-pi seam.
rootsFound = wrapToPiLocal(rootsFound);
rootsFound = sort(rootsFound);

keep = true(size(rootsFound));

for k = 2:numel(rootsFound)

    if abs(wrapToPiLocal(rootsFound(k)-rootsFound(k-1))) < 1e-6
        keep(k) = false;
    end
end

rootsFound = rootsFound(keep);

if numel(rootsFound) > 2

    % For a point outside a smooth convex offset ellipse there should be two
    % tangencies. If grid seam duplication survived, keep the two most
    % separated roots.
    bestPair = [1 2];
    bestSep = -inf;

    for i = 1:numel(rootsFound)-1
        for j = i+1:numel(rootsFound)

            sep = abs(wrapToPiLocal(rootsFound(j)-rootsFound(i)));

            if sep > bestSep
                bestSep = sep;
                bestPair = [i j];
            end
        end
    end

    rootsFound = rootsFound(bestPair);
end

thetaRoots = rootsFound(:);
Q = zeros(numel(thetaRoots),2);

for k = 1:numel(thetaRoots)
    Q(k,:) = femurOffsetPoint(thetaRoots(k),geo);
end

if isfield(geo, 'femurCondyleClipX')
    keep = Q(:,1) <= geo.femurCondyleClipX + 1e-10;
    Q = Q(keep,:);
    thetaRoots = thetaRoots(keep);
end

end


function root = bisectTangentRoot(P,a,b,geo)

[Qa,na] = femurOffsetPointNormal(a,geo);
fa = dot(P-Qa,na);

for iter = 1:40

    c = 0.5*(a+b);

    [Qc,nc] = femurOffsetPointNormal(c,geo);
    fc = dot(P-Qc,nc);

    if abs(fc) < 1e-12
        root = c;
        return
    end

    if fa*fc <= 0
        b = c;
    else
        a = c;
        fa = fc;
    end
end

root = 0.5*(a+b);

end


function Q = femurOffsetPoint(theta,geo)

[Q,~] = femurOffsetPointNormal(theta,geo);

end


function [Q,nWorld] = femurOffsetPointNormal(theta,geo)

a0 = geo.femurEllipsePhysicalA;
b0 = geo.femurEllipsePhysicalB;

ct = cos(theta);
st = sin(theta);

E = [a0*ct, b0*st];

nLocal = [ct/a0, st/b0];
nLocal = nLocal/norm(nLocal);

qLocal = E + geo.bpaRadius*nLocal;

th = geo.femurEllipseTheta;
Rell = [cos(th), -sin(th); sin(th), cos(th)];

Q = geo.femurProfileCenter + qLocal*Rell';
nWorld = nLocal*Rell';

end


function rho = offsetEllipseRadiusOfCurvature(theta,geo)
% Outward parallel curve radius of curvature = physical rho + offset.

a0 = geo.femurEllipsePhysicalA;
b0 = geo.femurEllipsePhysicalB;

ct = cos(theta);
st = sin(theta);

rhoPhysical = ((a0^2*st^2 + b0^2*ct^2)^(3/2))/(a0*b0);

rho = rhoPhysical + geo.bpaRadius;

end


function [q,theta] = projectSeedToFemurOffset(P,geo)
% Choose the exact normal-offset point whose physical-ellipse parameter is
% closest to the direction of P from the profile center.

th = geo.femurEllipseTheta;
Rell = [cos(th), -sin(th); sin(th), cos(th)];

local = (P-geo.femurProfileCenter)*Rell;

a0 = geo.femurEllipsePhysicalA;
b0 = geo.femurEllipsePhysicalB;

theta = atan2(local(2)/b0, local(1)/a0);
q = femurOffsetPoint(theta,geo);

end


function thetaMid = halfwayWrappedAngle(thetaA,thetaB)

d = wrapToPiLocal(thetaB-thetaA);
thetaMid = wrapToPiLocal(thetaA + 0.5*d);

end


function a = wrapToPiLocal(a)

a = mod(a+pi,2*pi)-pi;

end


%% =====================================================================
%% Small geometry utilities
%% =====================================================================

function tVals = segmentPolygonIntersectionParams(A,B,poly,tol)

tVals = zeros(0,1);

n = size(poly,1);

for e = 1:n

    C = poly(e,:);
    D = poly(mod(e,n)+1,:);

    t = segmentSegmentIntersectionParam(A,B,C,D,tol);

    if isfinite(t)
        tVals(end+1,1) = t; %#ok<AGROW>
    end
end

end


function t = segmentSegmentIntersectionParam(A,B,C,D,tol)

r = B-A;
s = D-C;

rxs = cross2(r,s);
qmp = C-A;

t = NaN;

if abs(rxs) <= tol
    return
end

tt = cross2(qmp,s)/rxs;
u = cross2(qmp,r)/rxs;

if tt >= -tol && tt <= 1+tol && u >= -tol && u <= 1+tol
    t = min(max(tt,0),1);
end

end


function [P,ok] = lineIntersection2D(A,B,C,D)

r = B-A;
s = D-C;

den = cross2(r,s);

if abs(den) < 1e-12
    P = [NaN NaN];
    ok = false;
    return
end

t = cross2(C-A,s)/den;

P = A+t*r;
ok = all(isfinite(P));

end


function t = pointParameterOnLine(A,B,Q)

d = B-A;
d2 = dot(d,d);

if d2 < 1e-18
    t = inf;
else
    t = dot(Q-A,d)/d2;
end

end


function z = cross2(a,b)

z = a(1)*b(2)-a(2)*b(1);

end


function q = projectSeedToCircle(P,C,R)

v = P-C;
n = norm(v);

if n < 1e-12
    v = [1 0];
    n = 1;
end

q = C + R*v/n;

end


function [bendMeasure, wrapInfo] = geometricBendMeasure( ...
    Location, T, T_Pam_inv, T_t1_ICR, phiD, alphaTol, geo, active, raw)
% bendMeasure is the geometric input to the X3 usable-length-loss model.
% Lmt stays as the piecewise-linear route length; the wrap terms below are
% only the extra sum(R*alpha) terms used by MonoPamDataExplicit_balanceX3.

N = size(Location,3);
M = size(Location,1);

labels = { ...
    'p2 femur cylinder'; ...
    'p6-p8 tibia curve'; ...
    'p3-p5 femur curve'; ...
    'p5-p6 ICR curve'};

nWrap = numel(labels);

bendMeasure = zeros(N,1);
wrapAngleRad = zeros(nWrap,N);
wrapLength = zeros(nWrap,N);
wrapRadius = zeros(nWrap,N);
pointAngleRad = zeros(M,N);
pointLength = zeros(M,N);

femurCylinderEntry = nan(N,2);
femurCylinderExit = nan(N,2);
tibiaCurveStart = nan(N,2);
tibiaCurveEnd = nan(N,2);
femurCondyleStart = nan(N,2);
femurCondyleEnd = nan(N,2);
femurCondyleRadiusA = nan(N,1);
femurCondyleRadiusB = nan(N,1);
crosspointStart = nan(N,2);
crosspointEnd = nan(N,2);

p2LineAngleRad = nan(N,1);
p2WrapFraction = nan(N,1);
tibiaCurveFraction = nan(N,1);
femurCurveFraction = nan(N,1);

Pall = zeros(M,3,N);

for i = 1:N
    Pall(:,:,i) = routeInFemurFrame(Location,T,i);
end

p2Ref = p2WrapReference(Pall, active, phiD, geo);
tibiaRef = tibiaCurveWrapReference(raw, active, phiD, geo);
femurRef = femurCurveWrapReference(raw, active, phiD, geo);
crossRef = crosspointWrapReference(raw, Location, active, phiD, T_Pam_inv);
for i = 1:N

    P = Pall(:,:,i);

    [alpha, R, qIn, qOut, ok, thetaLineRad, frac] = ...
        femurCylinderWrapAngle(P, active(:,i), p2Ref, phiD(i));

    if ok
        femurCylinderEntry(i,:) = qIn;
        femurCylinderExit(i,:) = qOut;
        p2LineAngleRad(i) = thetaLineRad;
        p2WrapFraction(i) = frac;
        wrapAngleRad(1,i) = alpha;
        wrapRadius(1,i) = R;
        pointAngleRad(2,i) = alpha;

        if alpha >= alphaTol
            wrapLength(1,i) = R*alpha;
            pointLength(2,i) = wrapLength(1,i);
        end
    end

    [alpha, R, qIn, qOut, ok, frac] = ...
        tibiaCurveWrapAngle(raw(:,:,i), active(:,i), tibiaRef, phiD(i));

    if ok
        tibiaCurveStart(i,:) = qIn;
        tibiaCurveEnd(i,:) = qOut;
        tibiaCurveFraction(i) = frac;
        wrapAngleRad(2,i) = alpha;
        wrapRadius(2,i) = R;
        pointAngleRad(7,i) = alpha;

        if alpha >= alphaTol
            wrapLength(2,i) = R*alpha;
            pointLength(7,i) = wrapLength(2,i);
        end
    end

    [alpha, R, qIn, qOut, ok, frac, Ra, Rb] = ...
        femurCurveWrapAngle(raw(:,:,i), active(:,i), femurRef, phiD(i), geo);

    if ok
        femurCondyleStart(i,:) = qIn;
        femurCondyleEnd(i,:) = qOut;
        femurCondyleRadiusA(i) = Ra;
        femurCondyleRadiusB(i) = Rb;
        femurCurveFraction(i) = frac;
        wrapAngleRad(3,i) = alpha;
        wrapRadius(3,i) = R;
        pointAngleRad(4,i) = alpha;

        if alpha >= alphaTol
            s = integratedVariableRadiusArc(femurRef, frac);
            wrapLength(3,i) = s;
            pointLength(4,i) = wrapLength(3,i);
        end
    end

    [alpha, R, qIn, qOut, ok] = ...
        crosspointICRWrapAngle(raw(:,:,i), Location(:,:,i), active(:,i), ...
        T_Pam_inv(:,:,i), crossRef, phiD(i));

    if ok
        crosspointStart(i,:) = qIn;
        crosspointEnd(i,:) = qOut;
        wrapAngleRad(4,i) = alpha;
        wrapRadius(4,i) = R;
        pointAngleRad(6,i) = alpha;

        if alpha >= alphaTol
            wrapLength(4,i) = R*alpha;
            pointLength(6,i) = wrapLength(4,i);
        end
    end

    bendMeasure(i) = sum(wrapLength(:,i));
end

wrapInfo.labels = labels;
wrapInfo.pointRows = [2;7;4;6];
wrapInfo.angleRad = wrapAngleRad;
wrapInfo.angleDeg = wrapAngleRad*180/pi;
wrapInfo.radius = wrapRadius;
wrapInfo.length = wrapLength;
wrapInfo.bendLossLength = wrapLength;
wrapInfo.pointAngleRad = pointAngleRad;
wrapInfo.pointAngleDeg = pointAngleRad*180/pi;
wrapInfo.pointLength = pointLength;

wrapInfo.femurCylinderEntry = femurCylinderEntry;
wrapInfo.femurCylinderExit = femurCylinderExit;
wrapInfo.tibiaCurveCenter = geo.tibiaCurveCenter;
wrapInfo.tibiaCurveStart = tibiaCurveStart;
wrapInfo.tibiaCurveEnd = tibiaCurveEnd;
wrapInfo.tibiaCurveFraction = tibiaCurveFraction;

wrapInfo.femurCondyleCenter = geo.femurProfileCenter;
wrapInfo.femurCondyleStart = femurCondyleStart;
wrapInfo.femurCondyleEnd = femurCondyleEnd;
wrapInfo.femurCondyleRadiusA = femurCondyleRadiusA;
wrapInfo.femurCondyleRadiusB = femurCondyleRadiusB;
wrapInfo.femurCurveFraction = femurCurveFraction;

wrapInfo.crosspointCenter = [0 0];
wrapInfo.crosspointStart = crosspointStart;
wrapInfo.crosspointEnd = crosspointEnd;

wrapInfo.p2LineAngleRad = p2LineAngleRad;
wrapInfo.p2LineAngleDeg = p2LineAngleRad*180/pi;
wrapInfo.p2WrapFraction = p2WrapFraction;
wrapInfo.p2ThetaBaseRad = p2Ref.thetaBaseRad;
wrapInfo.p2ThetaBaseDeg = p2Ref.thetaBaseRad*180/pi;
wrapInfo.p2ThetaMaxRad = p2Ref.thetaMaxRad;
wrapInfo.p2ThetaMaxDeg = p2Ref.thetaMaxRad*180/pi;
wrapInfo.p2AlphaBaseRad = p2Ref.alphaBaseRad;
wrapInfo.p2AlphaBaseDeg = p2Ref.alphaBaseRad*180/pi;
wrapInfo.p2AlphaMaxRad = p2Ref.alphaPlateauRad;
wrapInfo.p2AlphaMaxDeg = p2Ref.alphaPlateauRad*180/pi;

    function s = integratedVariableRadiusArc(ref, frac)

            frac = max(0, min(1, frac));
            n = 40;
            u = linspace(0, frac, n).';
            
            alphaSpan = ref.alphaSpanRad;
            r0 = ref.radiusBase;
            r1 = ref.radiusPlateau;
            drdu = r1 - r0;
            
            r = r0 + drdu*u;
            dtheta_du = alphaSpan;
            
            ds_du = sqrt((r.*dtheta_du).^2 + drdu.^2);
            s = trapz(u, ds_du);
    
        end

end


function P = routeInFemurFrame(Location,T,i)

M = size(Location,1);
P = zeros(M,3);

% Femur-side rows.
P(1:5,:) = Location(1:5,:,i);

% t1 rows are stored in Location in ICR coordinates. Convert ICR -> femur
% with T_Pam before calculating route angles.
for j = 6:9
    P(j,:) = RowVecTrans(T(:,:,i),Location(j,:,i));
end

end


function [alpha, R, qIn, qOut, ok, thetaLineRad, wrapFrac] = ...
    femurCylinderWrapAngle(P, activeFrame, ref, kneeAngleD)

alpha = 0;
R = 0;
qIn = [NaN NaN];
qOut = [NaN NaN];
ok = false;
thetaLineRad = NaN;
wrapFrac = NaN;

if ~activeFrame(2)
    return
end

if activeFrame(3)
    thetaLineRad = ref.thetaMaxRad;
    qIn = P(2,1:2);
    qOut = P(3,1:2);
    wrapFrac = 1;
    alpha = ref.alphaPlateauRad;
    R = ref.radiusPlateau;
    ok = ref.valid;
    return
end

% Before p3 is active, p2 wrap grows from the fully-extended p8*-p2
% baseline toward the p3 plateau. p8 is stored in t1, but P is already the
% common femur-frame route, so row 8 here is p8*.
[thetaLineRad, okTheta] = p2LineAngleToRow(P, 8);

if ~okTheta
    return
end

qIn = P(2,1:2);
qOut = P(8,1:2);

[alpha, R, wrapFrac] = interpolatedWABKA(kneeAngleD, ref);
ok = isfinite(alpha);

end


function ref = p2WrapReference(Pall, active, phiD, geo)

N = size(Pall,3);
idxExtension = N;

[thetaBaseRad, okBase] = p2LineAngleToRow(Pall(:,:,idxExtension), 8);

idxP3On = find(active(3,:), 1, 'last');

if isempty(idxP3On)
    idxP3On = 1;
end

[thetaMaxRad, okThetaMax] = p2LineAngleToRow(Pall(:,:,idxP3On), 8);
[alphaPlateauRad, okPlateau] = p2LineAngleToRow(Pall(:,:,idxP3On), 3);

ref = makeWABKARef( ...
    phiD(idxExtension), ...
    phiD(idxP3On), ...
    thetaBaseRad, ...
    alphaPlateauRad, ...
    geo.femurCylClearRadius, ...
    geo.femurCylClearRadius, ...
    okBase && okThetaMax && okPlateau);

ref.idxExtension = idxExtension;
ref.idxP3On = idxP3On;
ref.thetaBaseRad = thetaBaseRad;
ref.thetaMaxRad = thetaMaxRad;

end


function ref = tibiaCurveWrapReference(raw, active, phiD, geo)
% Temporary tibia wrap-loss model through p8, p7, and p6.  The knee-angle
% map starts at extension and reaches its plateau when p6 is introduced
% while reading the plot right-to-left.

N = size(raw,3);
idxExtension = N;

C = geo.tibiaCurveCenter;
seedPts = squeeze(raw([6 7 8],1:2,1));
R = mean(vecnorm(seedPts - C, 2, 2));

idxP6On = find(active(6,:), 1, 'last');

if isempty(idxP6On)
    idxP6On = 1;
end

[alphaBaseRad, okBase] = centralAngleUnsigned( ...
    C, raw(9,1:2,idxExtension), raw(8,1:2,idxExtension));
[alphaPlateauRad, okPlateau] = centralAngleUnsigned( ...
    C, raw(9,1:2,idxP6On), raw(6,1:2,idxP6On));

ref = makeWABKARef( ...
    phiD(idxExtension), ...
    phiD(idxP6On), ...
    alphaBaseRad, ...
    alphaPlateauRad, ...
    R, ...
    R, ...
    okBase && okPlateau);

ref.idxExtension = idxExtension;
ref.idxP6On = idxP6On;
ref.center = C;

end


function ref = femurCurveWrapReference(raw, active, phiD, geo)
% Temporary femur wrap-loss model through p3, p4, and p5.  It is zero until
% p3 is introduced and reaches its plateau when p5 is introduced.

idxP3On = find(active(3,:), 1, 'last');
idxP5On = find(active(5,:), 1, 'last');

if isempty(idxP3On)
    idxP3On = 1;
end

if isempty(idxP5On)
    idxP5On = idxP3On;
end

C = geo.femurProfileCenter;
R3 = norm(raw(3,1:2,idxP3On) - C);
R5 = norm(raw(5,1:2,idxP5On) - C);

[alphaPlateauRad, okPlateau] = centralAngleUnsigned( ...
    C, raw(3,1:2,idxP5On), raw(5,1:2,idxP5On));

ref = makeWABKARef( ...
    phiD(idxP3On), ...
    phiD(idxP5On), ...
    0, ...
    alphaPlateauRad, ...
    R3, ...
    R5, ...
    okPlateau);

ref.idxP3On = idxP3On;
ref.idxP5On = idxP5On;
ref.center = C;

end


function ref = makeWABKARef( ...
    phiBaseD, phiPlateauD, alphaBaseRad, alphaPlateauRad, ...
    radiusBase, radiusPlateau, ok)
% WABKA = wrap angle by knee angle.  The route state supplies the knee
% angles where a wrap term starts and plateaus; this struct stores the
% linear map between those knee-angle limits.

ref.valid = logical(ok) && ...
    all(isfinite([phiBaseD, phiPlateauD, alphaBaseRad, alphaPlateauRad, ...
                  radiusBase, radiusPlateau])) && ...
    radiusBase > 0 && radiusPlateau > 0 && ...
    abs(phiPlateauD - phiBaseD) > 1e-12;

ref.phiBaseD = phiBaseD;
ref.phiExtensionD = phiBaseD;
ref.phiPlateauD = phiPlateauD;
ref.phiSpanD = phiPlateauD - phiBaseD;
ref.alphaBaseRad = alphaBaseRad;
ref.alphaPlateauRad = alphaPlateauRad;
ref.alphaMaxRad = alphaPlateauRad;
ref.alphaSpanRad = alphaPlateauRad - alphaBaseRad;
ref.radiusBase = radiusBase;
ref.radiusPlateau = radiusPlateau;
ref.radiusSpan = radiusPlateau - radiusBase;

end


function [alpha, radius, frac] = interpolatedWABKA(kneeAngleD, ref)
% WABKA = wrap angle by knee angle.  A negative phiSpanD is expected here
% because the route is plotted from flexion on the left to extension on the
% right, while the wrap grows when reading the plot right-to-left.

alpha = 0;
radius = 0;
frac = 0;

if ~ref.valid
    return
end

frac = (kneeAngleD - ref.phiBaseD) / ref.phiSpanD;
frac = max(0, min(1, frac));

alpha = ref.alphaBaseRad + ref.alphaSpanRad*frac;
radius = ref.radiusBase + ref.radiusSpan*frac;

end


function [thetaRad, ok] = p2LineAngleToRow(P, row)

vIn = P(2,1:2) - P(1,1:2);
vOut = P(row,1:2) - P(2,1:2);

[thetaRad, ok] = vectorAngleUnsigned(vIn, vOut);

end


function [thetaRad, ok] = vectorAngleUnsigned(a, b)

na = norm(a);
nb = norm(b);

thetaRad = NaN;
ok = false;

if na < 1e-12 || nb < 1e-12
    return
end

ca = dot(a/na, b/nb);
ca = max(-1, min(1, ca));
thetaRad = acos(ca);
ok = true;

end


function [alpha, R, qIn, qOut, ok, frac] = ...
    tibiaCurveWrapAngle(rawFrame, activeFrame, ref, kneeAngleD)

alpha = 0;
R = 0;
qIn = [NaN NaN];
qOut = [NaN NaN];
ok = false;
frac = NaN;

if ~activeFrame(8)
    return
end

[alpha, R, frac] = interpolatedWABKA(kneeAngleD, ref);

if activeFrame(6)
    qIn = rawFrame(6,1:2);
    alpha = ref.alphaPlateauRad;
    R = ref.radiusPlateau;
    frac = 1;
elseif activeFrame(7)
    qIn = rawFrame(7,1:2);
else
    qIn = rawFrame(8,1:2);
end

qOut = rawFrame(9,1:2);
ok = ref.valid && isfinite(alpha) && isfinite(R);

end


function [alpha, R, qIn, qOut, ok, frac, Ra, Rb] = ...
    femurCurveWrapAngle(rawFrame, activeFrame, ref, kneeAngleD, geo)

alpha = 0;
R = 0;
qIn = [NaN NaN];
qOut = [NaN NaN];
ok = false;
frac = NaN;
Ra = 0;
Rb = 0;

if ~activeFrame(3)
    ok = ref.valid;
    return
end

[alpha, R, frac] = interpolatedWABKA(kneeAngleD, ref);

qIn = rawFrame(3,1:2);

if activeFrame(5)
    qOut = rawFrame(5,1:2);
    alpha = ref.alphaPlateauRad;
    R = ref.radiusPlateau;
    frac = 1;
elseif activeFrame(4)
    qOut = rawFrame(4,1:2);
else
    qOut = rawFrame(3,1:2);
end

Ra = norm(qIn - geo.femurProfileCenter);
Rb = norm(qOut - geo.femurProfileCenter);
ok = ref.valid && isfinite(alpha) && isfinite(R);

end

function ref = crosspointWrapReference(raw, Location, active, phiD, T_Pam_inv)

idxP5On = find(active(5,:), 1, 'last');
idxFlex = 1;

q5AtP5 = RowVecTrans(T_Pam_inv(:,:,idxP5On), raw(5,:,idxP5On));
q6AtP5 = Location(6,:,idxP5On);

q5Flex = RowVecTrans(T_Pam_inv(:,:,idxFlex), raw(5,:,idxFlex));
q6Flex = Location(6,:,idxFlex);

alphaAtP5 = 0;
okP5 = true;

[alphaFlex, okFlex] = centralAngleUnsigned([0 0], q5Flex(1:2), q6Flex(1:2));

RAtP5 = 0.5*(norm(q5AtP5(1:2)) + norm(q6AtP5(1:2)));
RFlex = 0.5*(norm(q5Flex(1:2)) + norm(q6Flex(1:2)));

ref = makeWABKARef(phiD(idxP5On), phiD(idxFlex), ...
    alphaAtP5, alphaFlex, RAtP5, RFlex, okP5 && okFlex);

end

function [alpha, R, qIn, qOut, ok] = ...
    crosspointICRWrapAngle(rawFrame, locFrame, activeFrame, T_Pam_inv_i, ref, kneeAngleD)

alpha = 0;
R = 0;
qIn = [NaN NaN];
qOut = [NaN NaN];
ok = false;

if ~(activeFrame(5) && activeFrame(6))
    return
end

q5ICR = RowVecTrans(T_Pam_inv_i, rawFrame(5,:));
q6ICR = locFrame(6,:);

[alpha, R] = interpolatedWABKA(kneeAngleD, ref);

qIn = q5ICR(1:2);
qOut = q6ICR(1:2);
ok = ref.valid && isfinite(alpha) && isfinite(R);

end


function [thetaRad, ok] = centralAngleUnsigned(C, A, B)

[thetaRad, ok] = vectorAngleUnsigned(A(1:2) - C, B(1:2) - C);

end


function alpha = circleMinorAngle(qA,qB,C)

a = qA-C;
b = qB-C;

na = norm(a);
nb = norm(b);

if na < 1e-12 || nb < 1e-12
    alpha = 0;
    return
end

ca = dot(a/na,b/nb);
ca = max(-1,min(1,ca));
alpha = acos(ca);

end
