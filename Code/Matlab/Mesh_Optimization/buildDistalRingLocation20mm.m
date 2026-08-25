function [Location, bendMeasure, info] = buildDistalRingLocation20mm(p1, p8, tendon, ctx)
% buildDistalRingLocation20mm
%
% 20 mm extensor route built by sweeping from extension toward flexion.
% The original Knee_Extensor_10mm p1:p8 route is used as the geometric seed.
% New contacts are created from the CURRENT straight route at first contact,
% then frozen in their native bone frame for all greater flexion.
%
% p1:p4 are femur-frame points.
% p5:p8 are t1-frame points and are transformed to ICR before Location output.
% Inactive points are represented by repeated neighboring active points.
%
% Contact order is constrained physically:
%   femur: p2 -> p3 -> p4
%   tibia: p7 -> p6 -> p5
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

%% Original route used as the seed / branch reference
% These are the original Knee_Extensor_10mm native-frame points.
seed = zeros(8,3);
seed(1,:) = [0.04000,  0.03500,  0.00000];
seed(2,:) = [0.07590, -0.27476,  0.00000];
seed(3,:) = [0.05982, -0.37427,  0.00000];
seed(4,:) = [0.03955, -0.42183,  0.00000];
seed(5,:) = [0.05871,  0.02500,  0.00000];
seed(6,:) = [0.05871,  0.01228,  0.00000];
seed(7,:) = [0.05400, -0.00612,  0.00000];
seed(8,:) = [0.03604, -0.02844,  0.00000];

% p1 and p8 are optimization variables, not fixed seed coordinates.
seed(1,:) = p1;
seed(8,:) = p8;

% 20 mm-adjusted seed references. These are only references/fallbacks;
% actual points are frozen from the first geometric contact.
seed20 = adjustedSeed20mm(seed, ctx.geo, p1, p8);

%% Output allocation
Location = zeros(8,3,N);
bendRadius = zeros(8,N);

info.raw = zeros(8,3,N);
info.candidateRaw = zeros(8,3,N);
info.active = false(8,N);
info.tendonMax = zeros(N,1);
info.endcap_t1 = zeros(N,3);
info.targetType = zeros(N,1);

info.pointFrame = ["femur";"femur";"femur";"femur"; ...
                   "t1";"t1";"t1";"t1"];
info.addedAngleD = nan(8,1);
info.addedPoint = nan(8,3);

%% Progressive contact state
active = false(8,1);
active([1 8]) = true;

% Once a point activates, it is frozen here in its native frame.
Pfixed = seed20;
Pfixed(1,:) = p1;
Pfixed(8,:) = p8;

info.addedAngleD(1) = sweepD(1);
info.addedAngleD(8) = sweepD(1);
info.addedPoint(1,:) = p1;
info.addedPoint(8,:) = p8;

contactTol = 1e-7;

for s = 1:Ns
    thetaD = sweepD(s);

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
        tendonLimit20mm(p8, ctx.geo);

    frameInfo.tendonMax = tendonMaxGeom;
    frameInfo.targetType = targetType;
    frameInfo.targetRadius = targetRadius;

    % Physical BPA endcap point along p8 -> tangent direction. Diagnostic only.
    [pTanGeom, ~, ~, ~, ~] = tendonLimit20mm(p8, ctx.geo);
    if tendonMaxGeom > 1e-12
        frac = min(max(tendon/tendonMaxGeom,0),1);
        endcapXY = p8(1:2) + frac*(pTanGeom - p8(1:2));
    else
        endcapXY = p8(1:2);
    end
    frameInfo.endcap_t1 = [endcapXY, p8(3)];

    % At each sweep angle, repeatedly add the NEXT allowable contact if the
    % current straight route reaches its clearance geometry. Because the new
    % point is taken from the current segment at contact, activating it does
    % not create a finite path-length jump.
    if s > 1
        for pass = 1:6
            raw = collapseInactive(Pfixed, active);

            [newIdx, newPoint] = nextContact( ...
                raw, active, T_Pam_i, T_ICR_t1_i, T_t1_ICR_i, ...
                ctx.geo, frameInfo, seed20, contactTol);

            if newIdx == 0
                break
            end

            active(newIdx) = true;
            Pfixed(newIdx,:) = newPoint;

            if isnan(info.addedAngleD(newIdx))
                info.addedAngleD(newIdx) = thetaD;
                info.addedPoint(newIdx,:) = newPoint;
            end
        end
    end

    raw = collapseInactive(Pfixed, active);

    % Active bend radii only.
    rBend = activeBendRadii(Pfixed, active, frameInfo, ctx.geo);

    % Repeated points are not physical bends.
    for j = 2:7
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
        bendRadius(:,solverIndex) = rBend;

        % Match Knee_Extensor_10mm: p5:p8 t1 -> ICR before Location output.
        loc = raw;
        for j = 5:8
            loc(j,:) = RowVecTrans(T_ICR_t1_i, loc(j,:));
        end
        Location(:,:,solverIndex) = loc;
    end
end

%% Xi3 geometric bend measure
bendMeasure = geometricBendMeasure(Location, ctx.T_Pam, bendRadius, ctx.geo.alphaTol);

info.bendRadius = bendRadius;
info.sweepD = sweepD;
info.sweepStepD = dPhiD;
info.seedOriginal = seed;
info.seed20 = seed20;

idx = (1:8)';
keyAngle = info.addedAngleD;
keyAngle(isnan(keyAngle)) = -Inf;
[~,ord] = sortrows([-keyAngle, idx],[1 2]);
info.activationOrder = idx(ord);

end


function seed20 = adjustedSeed20mm(seed, geo, p1, p8)
% Adjust the original route approximately onto the 20 mm BPA clearance
% geometry while preserving the original point directions/locations as much
% as possible. These values are branch references and safe fallbacks.

seed20 = seed;
seed20(1,:) = p1;
seed20(8,:) = p8;

% p2: same radial direction as original p2, but on expanded femur cylinder.
seed20(2,1:2) = projectSeedToCircle( ...
    seed(2,1:2), geo.femurCylCenter, geo.femurCylClearRadius);
seed20(2,3) = p1(3);

% p3: preserve original y, move x to expanded vertical clearance line.
seed20(3,1) = geo.femurLineX;
seed20(3,2) = min(max(seed(3,2), min(geo.femurLineY)), max(geo.femurLineY));
seed20(3,3) = p1(3);

% p4: preserve original radial direction from profile center, move to the
% expanded rotated ellipse.
seed20(4,1:2) = projectSeedToEllipse(seed(4,1:2), geo);
seed20(4,3) = p1(3);

% p5: preserve original direction from upper tibia cylinder center.
seed20(5,1:2) = projectSeedToCircle( ...
    seed(5,1:2), geo.tibiaUpperCenter, geo.tibiaUpperClearRadius);
seed20(5,3) = p8(3);

% p6: preserve original y, move x to expanded common tangent wall.
seed20(6,1) = geo.tibiaWallX;
seed20(6,2) = min(max(seed(6,2), min(geo.tibiaWallY)), max(geo.tibiaWallY));
seed20(6,3) = p8(3);

% p7 reference: preserve original direction relative to lower tibia cylinder.
seed20(7,1:2) = projectSeedToCircle( ...
    seed(7,1:2), geo.tibiaLowerCenter, geo.tibiaLowerClearRadius);
seed20(7,3) = p8(3);

end


function [newIdx, newPoint] = nextContact( ...
    raw, active, T_Pam_i, T_ICR_t1_i, T_t1_ICR_i, ...
    geo, frameInfo, seed20, tol)
% Return at most one new contact per pass. The caller rebuilds the route and
% calls again, preserving physical contact order.

newIdx = 0;
newPoint = [NaN NaN NaN];

% Current crossing segment endpoints in both native frames.
A_fem = raw(4,:);
B_t1 = raw(5,:);
B_icr = RowVecTrans(T_ICR_t1_i, B_t1);
B_fem = RowVecTrans(T_Pam_i, B_icr);

A_icr = RowVecTrans(T_Pam_i\eye(4), A_fem);
A_t1 = RowVecTrans(T_t1_ICR_i, A_icr);

%% Femur sequence p2 -> p3 -> p4
if ~active(2)
    [hit,q] = segmentCircleContactPoint( ...
        raw(1,1:2), B_fem(1:2), ...
        geo.femurCylCenter, geo.femurCylClearRadius, seed20(2,1:2), tol);
    if hit
        newIdx = 2;
        newPoint = [q, raw(1,3)];
        return
    end

elseif ~active(3)
    [hit,q] = segmentVerticalContactPoint( ...
        raw(2,1:2), B_fem(1:2), geo.femurLineX, geo.femurLineY, tol);
    if hit
        newIdx = 3;
        newPoint = [q, raw(2,3)];
        return
    end

elseif ~active(4)
    [hit,q] = segmentEllipseContactPoint( ...
        raw(3,1:2), B_fem(1:2), geo, seed20(4,1:2), tol);
    if hit
        newIdx = 4;
        newPoint = [q, raw(3,3)];
        return
    end
end

%% Tibia sequence p7 -> p6 -> p5
if ~active(7)
    switch frameInfo.targetType
        case 1
            [hit,q] = segmentCircleContactPoint( ...
                A_t1(1:2), raw(8,1:2), ...
                geo.tibiaLowerCenter, geo.tibiaLowerClearRadius, ...
                seed20(7,1:2), tol);
        case 2
            [hit,q] = segmentVerticalContactPoint( ...
                A_t1(1:2), raw(8,1:2), ...
                geo.tibiaWallX, geo.tibiaWallY, tol);
        case 3
            [hit,q] = segmentCircleContactPoint( ...
                A_t1(1:2), raw(8,1:2), ...
                geo.tibiaUpperCenter, geo.tibiaUpperClearRadius, ...
                seed20(7,1:2), tol);
        otherwise
            hit = false;
            q = [NaN NaN];
    end

    if hit
        newIdx = 7;
        newPoint = [q, raw(8,3)];
        return
    end

elseif ~active(6)
    [hit,q] = segmentVerticalContactPoint( ...
        A_t1(1:2), raw(7,1:2), geo.tibiaWallX, geo.tibiaWallY, tol);
    if hit
        newIdx = 6;
        newPoint = [q, raw(7,3)];
        return
    end

elseif ~active(5)
    [hit,q] = segmentCircleContactPoint( ...
        A_t1(1:2), raw(6,1:2), ...
        geo.tibiaUpperCenter, geo.tibiaUpperClearRadius, ...
        seed20(5,1:2), tol);
    if hit
        newIdx = 5;
        newPoint = [q, raw(6,3)];
        return
    end
end

end


function raw = collapseInactive(Pfixed, active)
raw = Pfixed;

raw(1,:) = Pfixed(1,:);
for j = 2:4
    if active(j)
        raw(j,:) = Pfixed(j,:);
    else
        raw(j,:) = raw(j-1,:);
    end
end

raw(8,:) = Pfixed(8,:);
for j = 7:-1:5
    if active(j)
        raw(j,:) = Pfixed(j,:);
    else
        raw(j,:) = raw(j+1,:);
    end
end
end


function r = activeBendRadii(Pfixed, active, frameInfo, geo)
r = zeros(8,1);

if active(2)
    r(2) = geo.femurCylClearRadius;
end
if active(3)
    r(3) = geo.femurFilletClearRadius;
end
if active(4)
    r(4) = ellipseRadiusOfCurvature(Pfixed(4,1:2), geo);
end
if active(5)
    r(5) = geo.tibiaUpperClearRadius;
end
if active(6) && frameInfo.targetType ~= 2
    r(6) = frameInfo.targetRadius;
end
end


function [hit,q] = segmentCircleContactPoint(A,B,C,R,qRef,tol)
% Detect contact/penetration and return a surface point on the current
% straight segment direction. At first contact this is the tangent point.

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


function [hit,q] = segmentEllipseContactPoint(A,B,geo,qRef,tol)
% Segment / expanded rotated ellipse intersection. When two intersections
% exist, choose the one closest to the original p4 seed branch.

th = geo.femurEllipseTheta;
Rell = [cos(th), -sin(th); sin(th), cos(th)];
a = geo.femurEllipseA;
b = geo.femurEllipseB;

A0 = (A-geo.femurProfileCenter)*Rell;
B0 = (B-geo.femurProfileCenter)*Rell;
d = B0-A0;

qa = (d(1)/a)^2 + (d(2)/b)^2;
qb = 2*(A0(1)*d(1)/a^2 + A0(2)*d(2)/b^2);
qc = (A0(1)/a)^2 + (A0(2)/b)^2 - 1;

if qa < 1e-18
    hit = false;
    q = [NaN NaN];
    return
end

disc = qb^2 - 4*qa*qc;
if disc < -tol
    hit = false;
    q = [NaN NaN];
    return
end

disc = max(disc,0);
rdisc = sqrt(disc);
tVals = [(-qb-rdisc)/(2*qa), (-qb+rdisc)/(2*qa)];
tVals = tVals(tVals >= -tol & tVals <= 1+tol);

if isempty(tVals)
    hit = false;
    q = [NaN NaN];
    return
end

Q = zeros(numel(tVals),2);
for k = 1:numel(tVals)
    qLocal = A0 + tVals(k)*d;
    Q(k,:) = geo.femurProfileCenter + qLocal*Rell';
end

[~,kBest] = min(vecnorm(Q-qRef,2,2));
q = Q(kBest,:);
hit = true;
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


function q = projectSeedToEllipse(P,geo)
th = geo.femurEllipseTheta;
Rell = [cos(th), -sin(th); sin(th), cos(th)];
local = (P-geo.femurProfileCenter)*Rell;
a = geo.femurEllipseA;
b = geo.femurEllipseB;
val = (local(1)/a)^2 + (local(2)/b)^2;
if val < 1e-18
    local = [0,-b];
else
    local = local/sqrt(val);
end
q = geo.femurProfileCenter + local*Rell';
end


function rho = ellipseRadiusOfCurvature(P,geo)
th = geo.femurEllipseTheta;
Rell = [cos(th), -sin(th); sin(th), cos(th)];
local = (P-geo.femurProfileCenter)*Rell;
a = geo.femurEllipseA;
b = geo.femurEllipseB;

ct = local(1)/a;
st = local(2)/b;
n = hypot(ct,st);
if n > 1e-12
    ct = ct/n;
    st = st/n;
end
rho = ((a^2*st^2 + b^2*ct^2)^(3/2))/(a*b);
end


function [T_Pam_i,T_ICR_t1_i,T_t1_ICR_i] = extensionPrerollTransforms(thetaD,ctx)
% Above the measured solver range, preserve the max-extension translations
% and ICR offset and continue only the knee rotation.
T_Pam_i = ctx.T_Pam(:,:,end);
ct = cosd(thetaD);
st = sind(thetaD);
T_Pam_i(1:3,1:3) = [ct,-st,0; st,ct,0; 0,0,1];
T_ICR_t1_i = ctx.T_ICR_t1(:,:,end);
T_t1_ICR_i = ctx.T_t1_ICR(:,:,end);
end


function bendMeasure = geometricBendMeasure(Location,T,bendRadius,alphaTol)
N = size(Location,3);
bendMeasure = zeros(N,1);

for i = 1:N
    P = zeros(8,3);
    P(1:4,:) = Location(1:4,:,i);

    for j = 5:8
        P(j,:) = RowVecTrans(T(:,:,i),Location(j,:,i));
    end

    for j = 2:7
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
