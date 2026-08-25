function [Location, bendMeasure, info] = buildDistalRingLocation20mm(p1, pEnd, tendon, ctx)
% buildDistalRingLocation20mm
%
% 20 mm extensor route built by sweeping from extension toward flexion.
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
% Inactive points are represented by repeated neighboring active points.
% Once a contact activates it is frozen in its native bone frame.
%
% At the +10 deg solver-extension frame the seed route already contains:
%   p1, p2, p7, p8, p9.
%
% Expected geometric progression while flexing (NOT hard-coded angles):
%   p3 around -6 deg
%   p4 around -45 deg
%   p6 around -60 deg
%   p5 around -85 deg
%
% p3:p5 are successive femoral-body/condyle avoidance vertices.
% p6 is the later upper-tibia contact. p7:p8 are already present at extension.
%
% The hidden sweep uses the exact phiD increment. phi/phiD are not modified.

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

%% Grid-aligned hidden sweep from approximately +90 deg toward flexion
phiMax = max(phiD);
nAbove = max(0, ceil((90 - phiMax)/dPhiD));

if nAbove > 0
    preSweepD = phiMax + (nAbove:-1:1)'*dPhiD;
else
    preSweepD = zeros(0,1);
end

sweepD = [preSweepD; flipud(phiD)];
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
seed20 = adjustedSeed20mm(seed, p1, pEnd, ctx.geo);

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

info.CrossPoint = 6;
info.routeRows = 9;
info.pointFrame = [ ...
    "femur";"femur";"femur";"femur";"femur"; ...
    "t1";"t1";"t1";"t1"];

info.addedAngleD = nan(9,1);
info.addedPoint = nan(9,3);
info.addedSweepIndex = nan(9,1);

%% Progressive contact state
% The unsupported hidden pre-roll starts with only the two true endpoints.
% At the first measured solver frame (+10 deg), p2, p7 and p8 are installed
% from the SolidWorks seed, as specified by the physical design.
active = false(9,1);
active([1 9]) = true;
baseInitialized = false;

% Once a point activates, it is frozen here in its native frame.
Pfixed = seed20;
Pfixed(1,:) = p1;
Pfixed(9,:) = pEnd;

% Bend radius associated with each fixed physical contact.
fixedBendRadius = zeros(9,1);

info.addedAngleD(1) = sweepD(1);
info.addedAngleD(9) = sweepD(1);
info.addedPoint(1,:) = p1;
info.addedPoint(9,:) = pEnd;
info.addedSweepIndex(1) = 1;
info.addedSweepIndex(9) = 1;

contactTol = 1e-8;

for s = 1:Ns

    thetaD = sweepD(s);
    baseAddedThisSweep = false;

    if s <= nAbove
        [T_Pam_i, T_ICR_t1_i, T_t1_ICR_i] = ...
            extensionPrerollTransforms(thetaD, ctx);
        solverIndex = 0;
    else
        k = s - nAbove;
        solverIndex = N - k + 1;
        T_Pam_i = ctx.T_Pam(:,:,solverIndex);
        T_ICR_t1_i = ctx.T_ICR_t1(:,:,solverIndex);
        T_t1_ICR_i = ctx.T_t1_ICR(:,:,solverIndex);
    end

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

    %% Baseline contacts that physically exist at +10 deg
    if solverIndex > 0 && ~baseInitialized

        active([2 7 8]) = true;

        fixedBendRadius(2) = ctx.geo.femurCylClearRadius;
        fixedBendRadius(7) = 0; % straight tibia-wall contact
        fixedBendRadius(8) = ctx.geo.tibiaLowerClearRadius;

        for jj = [2 7 8]
            info.addedAngleD(jj) = thetaD;
            info.addedPoint(jj,:) = Pfixed(jj,:);
            info.addedSweepIndex(jj) = s;
        end

        baseInitialized = true;
        baseAddedThisSweep = true;
    end

    % Add the next physical contact, rebuild the repeated-point route, then
    % test again. Only one new point from the same curved obstacle is allowed
    % per sweep frame.
    if solverIndex > 0 && baseInitialized && ~baseAddedThisSweep

        activatedThisSweep = false(9,1);

        for pass = 1:6

            raw = collapseInactive(Pfixed, active);

            [newIdx, newPoint, newRadius] = nextContact( ...
                raw, active, activatedThisSweep, ...
                T_Pam_i, T_ICR_t1_i, T_t1_ICR_i, ...
                ctx.geo, seed20, contactTol);

            if newIdx == 0
                break
            end

            active(newIdx) = true;
            activatedThisSweep(newIdx) = true;
            Pfixed(newIdx,:) = newPoint;
            fixedBendRadius(newIdx) = newRadius;

            if isnan(info.addedAngleD(newIdx))
                info.addedAngleD(newIdx) = thetaD;
                info.addedPoint(newIdx,:) = newPoint;
                info.addedSweepIndex(newIdx) = s;
            end
        end
    end

    raw = collapseInactive(Pfixed, active);

    % Check the actual femur -> t1 crossing every frame. This catches a
    % visually obvious condyle penetration even when p5 has not activated.
    B_t1 = raw(6,:);
    B_icr = RowVecTrans(T_ICR_t1_i, B_t1);
    B_fem = RowVecTrans(T_Pam_i, B_icr);

    crossingCollision = segmentPenetratesFemurOffset( ...
        raw(5,1:2), B_fem(1:2), ctx.geo, contactTol);

    % If all available femur contacts have already been used and the
    % crossing still penetrates the corrected clearance, another femur row
    % would be required.
    overflow = active(5) && crossingCollision;

    rBend = fixedBendRadius;
    rBend(~active) = 0;

    % Repeated points are not physical bends.
    for j = 2:8
        if norm(raw(j,:) - raw(j-1,:)) < 1e-10 || ...
                norm(raw(j+1,:) - raw(j,:)) < 1e-10
            rBend(j) = 0;
        end
    end

    if solverIndex > 0

        info.raw(:,:,solverIndex) = raw;
        info.candidateRaw(:,:,solverIndex) = Pfixed;
        info.active(:,solverIndex) = active;
        info.tendonMax(solverIndex) = frameInfo.tendonMax;
        info.endcap_t1(solverIndex,:) = frameInfo.endcap_t1;
        info.targetType(solverIndex) = frameInfo.targetType;
        info.femurOverflow(solverIndex) = overflow;
        info.femurCrossingCollision(solverIndex) = crossingCollision;
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
end

%% Xi3 geometric bend measure
bendMeasure = geometricBendMeasure( ...
    Location, ctx.T_Pam, bendRadius, ctx.geo.alphaTol);

info.bendRadius = bendRadius;
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

end


%% =====================================================================
%% Seed geometry
%% =====================================================================

function seed20 = adjustedSeed20mm(seed, p1, pEnd, geo)
% Preserve the SolidWorks route as the branch reference while allowing the
% two baseline tangent contacts that should respond to endpoint motion to
% move with p1 and pEnd.
%
% p2: femur-cylinder contact.  This should only move slightly as p1 moves.
% p8: lower-tibia-cylinder contact.  This can move substantially when pEnd
%     moves, which is why it must not remain frozen at the original seed.
%
% p3:p7 remain the SolidWorks seed/reference coordinates until their own
% geometric activation logic creates/fixes them.

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

[q8, ok8] = tangentPointCircleNearest( ...
    pEnd(1:2), geo.tibiaLowerCenter, geo.tibiaLowerClearRadius, seed(8,1:2));

if ok8
    seed20(8,1:2) = q8;
end

end


%% =====================================================================
%% Progressive contact detection
%% =====================================================================

function [newIdx, newPoint, newRadius] = nextContact( raw, active, activatedThisSweep, ...
    T_Pam_i, T_ICR_t1_i, T_t1_ICR_i, geo, seed20, tol)
% Return at most one new contact per pass. The caller rebuilds the route
% before asking for another contact.
%
% p2, p7 and p8 already exist at +10 deg. The geometry then decides when
% p3, p4, p6 and p5 become necessary.

newIdx = 0;
newPoint = [NaN NaN NaN];
newRadius = 0;

% Current crossing segment endpoints in both native frames.
A_fem = raw(5,:);
B_t1 = raw(6,:);

B_icr = RowVecTrans(T_ICR_t1_i, B_t1);
B_fem = RowVecTrans(T_Pam_i, B_icr);

A_icr = RowVecTrans(T_Pam_i\eye(4), A_fem);
A_t1 = RowVecTrans(T_t1_ICR_i, A_icr);

%% p3: first femoral-body/condyle trigger
% Use the condyle-center -> p3 distance as the trigger radius.
% When the p2 -> tibia-side crossing enters this radius, activate the
% fixed CAD p3 seed directly.
if ~active(3)

    if segmentPenetratesCircle( ...
            raw(2,1:2), B_fem(1:2), ...
            geo.p3TriggerCenter, geo.p3TriggerRadius, tol)

        newIdx = 3;
        newPoint = [seed20(3,1:2), raw(2,3)];
        newRadius = geo.p3TriggerRadius;
        return
    end
end

%% p4: later femoral-condyle contact
% After p3 exists, check farther along the outgoing segment before deciding
% whether the true condyle clearance envelope needs the next seed point.
    if active(3) && ~active(4) && ~activatedThisSweep(3)
    
        A3 = raw(3,1:2);
        B3 = B_fem(1:2);
        Acheck = A3 + 0.25*(B3-A3);
    
        if segmentPenetratesFemurOffset(Acheck, B3, geo, tol)
    
            newIdx = 4;
            newPoint = [seed20(4,1:2), raw(3,3)];
            newRadius = norm(seed20(4,1:2) - geo.femurProfileCenter);
            return
        end
    end

%% p6: later tibia contact from the 3-point arc p6-p7-p8
    if active(4) && ~active(6)
    
        % With p5 inactive, raw(5) repeats p4. Express that femur-side crossing
        % point in t1 and test its line to the already-active p7 contact.
        A4_icr = RowVecTrans(T_Pam_i\eye(4), raw(5,:));
        A4_t1 = RowVecTrans(T_t1_ICR_i, A4_icr);
    
    hit = segmentPenetratesCircle(A4_t1(1:2), raw(7,1:2), ...
    geo.tibiaArcCenter, geo.tibiaArcRadius, tol);

    if hit
            newIdx = 6;
            newPoint = [seed20(6,1:2), raw(7,3)];
            newRadius = geo.tibiaArcRadius;
            return
        end
    end

%% p5: final high-flexion femoral-condyle / semi-minor clearance contact
    if active(4) && active(6) && ~active(5) && ~activatedThisSweep(4)
    
        B6_icr = RowVecTrans(T_ICR_t1_i, raw(6,:));
        B6_fem = RowVecTrans(T_Pam_i, B6_icr);
    
        A4 = raw(4,1:2);
        B4 = B6_fem(1:2);
        Acheck = A4 + 0.25*(B4-A4);
    
        if segmentPenetratesFemurOffset(Acheck, B4, geo, tol)
    
            newIdx = 5;
            newPoint = [seed20(5,1:2), raw(4,3)];
            newRadius = norm(seed20(5,1:2) - geo.femurProfileCenter);
            return
        end
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


%% =====================================================================
%% Corrected femoral-condyle clearance / tangent helpers
%% =====================================================================

function hit = segmentPenetratesFemurOffset(A,B,geo,tol)
% Strict collision test against the TRUE normal-offset ellipse polygon.
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
% Tangents are calculated against the TRUE normal-offset curve:
%   physical ellipse point + BPA_radius * outward_unit_normal.
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
% Find both external tangencies from P to the exact normal-offset ellipse.
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


%% =====================================================================
%% Extension pre-roll / Xi3 geometry
%% =====================================================================

function [T_Pam_i,T_ICR_t1_i,T_t1_ICR_i] = ...
    extensionPrerollTransforms(thetaD,ctx)
% Above the measured solver range, preserve the max-extension translations
% and ICR offset and continue only the knee rotation.

T_Pam_i = ctx.T_Pam(:,:,end);

ct = cosd(thetaD);
st = sind(thetaD);

T_Pam_i(1:3,1:3) = [ct,-st,0; st,ct,0; 0,0,1];

T_ICR_t1_i = ctx.T_ICR_t1(:,:,end);
T_t1_ICR_i = ctx.T_t1_ICR(:,:,end);

end


function bendMeasure = geometricBendMeasure( ...
    Location,T,bendRadius,alphaTol)

N = size(Location,3);
M = size(Location,1);

bendMeasure = zeros(N,1);

for i = 1:N

    P = zeros(M,3);

    % Femur-side rows.
    P(1:5,:) = Location(1:5,:,i);

    % t1 rows are stored in Location in ICR coordinates. Convert ICR ->
    % femur with T_Pam before calculating local bend angles.
    for j = 6:9
        P(j,:) = RowVecTrans(T(:,:,i),Location(j,:,i));
    end

    for j = 2:M-1

        Rj = bendRadius(j,i);

        if Rj <= 0
            continue
        end

        v1 = P(j,:)-P(j-1,:);
        v2 = P(j+1,:)-P(j,:);

        n1 = norm(v1);
        n2 = norm(v2);

        if n1 < 1e-9 || n2 < 1e-9
            continue
        end

        ca = dot(v1/n1,v2/n2);
        ca = max(-1,min(1,ca));

        alpha = acos(ca);

        if alpha < alphaTol
            continue
        end

        bendMeasure(i) = bendMeasure(i) + Rj*alpha;
    end
end

end
