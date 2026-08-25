%% TestDynamicRoutePoints20mm
% Standalone diagnostic for the 20 mm knee-extensor route.
%
% PURPOSE
%   Test the idea that the geometry, rather than a fixed p1:p8 state machine,
%   decides how many route/contact points are needed in the femur frame and
%   in the t1 frame at each knee angle.
%
% IMPORTANT
%   - This script DOES NOT replace buildDistalRingLocation20mm.m.
%   - It does not change the optimizer.
%   - It independently constructs a shortest collision-free polyline through
%     sampled clearance-envelope boundaries at each knee angle.
%   - The number of selected boundary points therefore emerges from the
%     geometry. A curved wrap can legitimately use several points.
%   - All collision testing is done in the femur frame. t1 obstacles are
%     transformed t1 -> ICR -> femur using the existing model transforms.
%
% OUTPUTS LEFT IN WORKSPACE
%   routeDynamic        100x1 struct, variable-length route for every angle
%   routeDynamicSummary table with point counts/path length
%
% If xBest exists in the workspace, the script tests xBest. Otherwise it
% tests ctx.x0.

clc
close all

%% Existing model context / design
if ~exist('ctx','var') || ~isstruct(ctx) || ~isfield(ctx,'geo')
    ctx = buildKneeExtContext20mm();
end

if exist('xBest','var') && isnumeric(xBest) && numel(xBest) >= 8
    xRouteTest = xBest(:).';
    designName = 'xBest';
else
    xRouteTest = ctx.x0(:).';
    designName = 'ctx.x0';
end

p1 = xRouteTest(1:3);
p8 = xRouteTest(4:6);

fprintf('\n============================================================\n')
fprintf('DYNAMIC ROUTE-POINT TEST\n')
fprintf('Design: %s\n', designName)
fprintf('p1 = [%.6f %.6f %.6f] m, femur frame\n', p1)
fprintf('p8 = [%.6f %.6f %.6f] m, t1 frame\n', p8)
fprintf('============================================================\n')

%% Diagnostic resolution
% More samples give a closer curved-envelope approximation, but create more
% candidate graph nodes and take longer. These values are a reasonable first
% test for all 100 knee positions.
nCircle  = 36;   % boundary vertices per circle
nEllipse = 48;   % boundary vertices on femur ellipse
nWall    = 12;   % candidate points along each exposed vertical wall

% Number of interior samples used to decide whether a straight graph edge
% passes through an obstacle.
nEdgeProbe = 9;

% Small geometry tolerance used only for inside/on-boundary classification.
insideTol = 1e-10;

% Plot a manageable subset of the 100 angles.
plotAnglesD = [10 0 -10 -30 -60 -90 -120];
overlayCurrentBuilder = true;

%% Optional current-builder prediction for visual comparison only
predCurrent = [];
if overlayCurrentBuilder
    try
        predCurrent = predictKneeExt20mm(xRouteTest, ctx);
        if ~predCurrent.ok
            predCurrent = [];
        end
    catch
        predCurrent = [];
    end
end

%% Solve a variable-length route independently at each solver angle
N = ctx.N;
routeDynamic = repmat(struct( ...
    'phiD',NaN, ...
    'ok',false, ...
    'pointsFemur',[], ...
    'pointsNative',[], ...
    'pointFrame',strings(0,1), ...
    'pointObstacle',strings(0,1), ...
    'pathLength',NaN, ...
    'nFemurContacts',NaN, ...
    'nTibiaContacts',NaN, ...
    'nGeneratedContacts',NaN, ...
    'frameSwitches',NaN, ...
    'obstacles',[]), N, 1);

nFemur = nan(N,1);
nTibia = nan(N,1);
nGenerated = nan(N,1);
pathLength = nan(N,1);
frameSwitches = nan(N,1);
okRoute = false(N,1);

for i = 1:N
    obs = buildDynamicObstacles( ...
        ctx.geo, p1, p8, ...
        ctx.T_Pam(:,:,i), ...
        ctx.T_ICR_t1(:,:,i), ...
        nCircle, nEllipse, nWall);

    % p8 expressed in femur frame for the global 2-D visibility problem.
    p8_icr = RowVecTrans(ctx.T_ICR_t1(:,:,i), p8);
    p8_fem = RowVecTrans(ctx.T_Pam(:,:,i), p8_icr);

    [nodes, G] = buildVisibilityGraph( ...
        p1, p8, p8_fem, obs, nEdgeProbe, insideTol);

    [nodePath, distPath] = shortestpath(G, 1, 2);

    routeDynamic(i).phiD = ctx.phiD(i);
    routeDynamic(i).obstacles = obs;

    if isempty(nodePath) || ~isfinite(distPath)
        continue
    end

    % Remove any graph nodes that can be skipped without cutting through an
    % obstacle. Curved wraps will retain as many samples as they actually need.
    nodePath = simplifyDynamicPath( ...
        nodePath, nodes, obs, nEdgeProbe, insideTol);

    P = nodes.xy(nodePath,:);
    Pnative = nodes.nativeXYZ(nodePath,:);
    frames = nodes.frame(nodePath);
    names = nodes.obstacle(nodePath);

    % 3-D path length: femur-frame XY plus the native-frame z assignment.
    % All candidate obstacle points inherit p1.z or p8.z depending on frame;
    % after transformation, p8-side points carry the correct transformed z.
    P3 = zeros(size(P,1),3);
    P3(:,1:2) = P;
    for k = 1:size(P,1)
        if frames(k) == "femur"
            P3(k,3) = p1(3);
        else
            q_icr = RowVecTrans(ctx.T_ICR_t1(:,:,i), Pnative(k,:));
            q_fem = RowVecTrans(ctx.T_Pam(:,:,i), q_icr);
            P3(k,3) = q_fem(3);
        end
    end

    L = sum(vecnorm(diff(P3,1,1),2,2));

    % Generated contacts exclude the permanent endpoints p1 and p8.
    isEndpoint = names == "p1" | names == "p8";
    nf = sum(frames == "femur" & ~isEndpoint);
    nt = sum(frames == "t1" & ~isEndpoint);

    fs = sum(frames(2:end) ~= frames(1:end-1));

    routeDynamic(i).ok = true;
    routeDynamic(i).pointsFemur = P3;
    routeDynamic(i).pointsNative = Pnative;
    routeDynamic(i).pointFrame = frames;
    routeDynamic(i).pointObstacle = names;
    routeDynamic(i).pathLength = L;
    routeDynamic(i).nFemurContacts = nf;
    routeDynamic(i).nTibiaContacts = nt;
    routeDynamic(i).nGeneratedContacts = nf + nt;
    routeDynamic(i).frameSwitches = fs;

    okRoute(i) = true;
    nFemur(i) = nf;
    nTibia(i) = nt;
    nGenerated(i) = nf + nt;
    pathLength(i) = L;
    frameSwitches(i) = fs;
end

routeDynamicSummary = table( ...
    ctx.phiD(:), okRoute, nFemur, nTibia, nGenerated, pathLength, frameSwitches, ...
    'VariableNames', { ...
    'phiD','ok','FemurContacts','TibiaContacts','GeneratedContacts', ...
    'PathLength_m','FrameSwitches'});

%% Print only rows where the point-count state changes
fprintf('\n========== DYNAMIC CONTACT-COUNT CHANGES ==========\n')

lastState = [NaN NaN NaN];
for i = N:-1:1   % extension -> flexion in print order
    if ~okRoute(i)
        fprintf('phi = %9.4f deg : NO COLLISION-FREE GRAPH PATH FOUND\n', ctx.phiD(i))
        continue
    end

    state = [nFemur(i), nTibia(i), frameSwitches(i)];
    if any(state ~= lastState)
        fprintf(['phi = %9.4f deg : femur contacts = %2d, ' ...
                 't1 contacts = %2d, total generated = %2d, ' ...
                 'frame switches = %d\n'], ...
            ctx.phiD(i), nFemur(i), nTibia(i), nGenerated(i), frameSwitches(i));
        lastState = state;
    end
end
fprintf('===================================================\n')

%% Warn if the shortest path changes native frame more than once
badSwitch = find(okRoute & frameSwitches > 1);
if ~isempty(badSwitch)
    fprintf('\nNOTE: %d of %d routes switch femur/t1 frame more than once.\n', ...
        numel(badSwitch), N)
    fprintf(['That does not automatically mean the geometry is wrong, but it is a ' ...
             'useful flag before converting this diagnostic into the solver route.\n'])
end

%% Plot number of automatically generated contacts vs knee angle
figure('Name','Dynamic route point count','Color','w')
hold on
plot(ctx.phiD, nFemur, 'o-', 'LineWidth', 1.5)
plot(ctx.phiD, nTibia, 's-', 'LineWidth', 1.5)
plot(ctx.phiD, nGenerated, '^-', 'LineWidth', 1.5)
grid on
xlabel('Knee angle, deg')
ylabel('Automatically generated contact points')
legend('Femur-frame contacts','t1-frame contacts','Total generated contacts', ...
    'Location','best')
title('Geometry-selected route complexity')

%% Plot dynamic route snapshots
figure('Name','Dynamic route snapshots','Color','w')
tiledlayout(2,4,'TileSpacing','compact','Padding','compact')

for q = 1:numel(plotAnglesD)
    [~,i] = min(abs(ctx.phiD - plotAnglesD(q)));
    nexttile
    hold on

    obs = routeDynamic(i).obstacles;
    for m = 1:numel(obs)
        xy = obs(m).polyXY;
        plot([xy(:,1);xy(1,1)], [xy(:,2);xy(1,2)], '-', 'LineWidth', 1)
    end

    if routeDynamic(i).ok
        P = routeDynamic(i).pointsFemur;
        plot(P(:,1), P(:,2), 'o-', 'LineWidth', 2, 'MarkerSize', 4)

        for k = 1:size(P,1)
            fr = routeDynamic(i).pointFrame(k);
            ob = routeDynamic(i).pointObstacle(k);
            text(P(k,1), P(k,2), sprintf(' %d:%s/%s', k, fr, ob), ...
                'FontSize', 7, 'Interpreter','none')
        end
    end

    % Overlay the existing fixed-state builder for comparison.
    if ~isempty(predCurrent)
        Pc = predCurrent.Location(:,:,i);
        for j = 5:8
            Pc(j,:) = RowVecTrans(ctx.T_Pam(:,:,i), Pc(j,:));
        end
        plot(Pc(:,1), Pc(:,2), '--', 'LineWidth', 1.2)
    end

    axis equal
    grid on
    xlabel('Femur-frame x, m')
    ylabel('Femur-frame y, m')
    title(sprintf('\\phi = %.2f deg | F=%g, T=%g', ...
        ctx.phiD(i), nFemur(i), nTibia(i)))
end

nexttile
axis off
text(0,1, { ...
    'Solid line/markers: dynamic visibility route', ...
    'Dashed line: current buildDistalRingLocation20mm route', ...
    'F = generated femur-frame contacts', ...
    'T = generated t1-frame contacts', ...
    '', ...
    'routeDynamic(i).pointsNative contains the', ...
    'native-frame coordinates of every selected point.'}, ...
    'VerticalAlignment','top','Interpreter','none')

fprintf('\nFull results are in routeDynamic and routeDynamicSummary.\n')
fprintf('The script has NOT changed buildDistalRingLocation20mm.m or the optimizer.\n')

%% =====================================================================
%% Local functions
%% =====================================================================

function obs = buildDynamicObstacles( ...
    geo, p1, p8, T_Pam_i, T_ICR_t1_i, nCircle, nEllipse, nWall)
% Build all clearance obstacles as polygons in the femur frame, while also
% retaining the candidate boundary points in each obstacle's native frame.

obs = struct('name',{},'frame',{},'polyXY',{},'candXY',{},'candNativeXYZ',{});

%% Femur cylinder
[poly,cand] = circleBoundary(geo.femurCylCenter, geo.femurCylClearRadius, nCircle);
obs(end+1) = makeFemurObstacle("femurCylinder", poly, cand, p1(3)); %#ok<AGROW>

%% Femur vertical clearance body
% The physical body is on the -x side of the exposed clearance line.
xLeftFemur = min([p1(1), geo.femurProfileCenter(1), -0.10]) - 0.15;
[poly,cand] = verticalWallBoundary( ...
    xLeftFemur, geo.femurLineX, geo.femurLineY, nWall);
obs(end+1) = makeFemurObstacle("femurWall", poly, cand, p1(3)); %#ok<AGROW>

%% Femur expanded condyle ellipse
[poly,cand] = ellipseBoundary(geo, nEllipse);
obs(end+1) = makeFemurObstacle("femurEllipse", poly, cand, p1(3)); %#ok<AGROW>

%% Tibia lower expanded circle, native t1
[polyT,candT] = circleBoundary( ...
    geo.tibiaLowerCenter, geo.tibiaLowerClearRadius, nCircle);
obs(end+1) = makeTibiaObstacle("tibiaLower", polyT, candT, p8(3), ...
    T_Pam_i, T_ICR_t1_i); %#ok<AGROW>

%% Tibia vertical clearance body, native t1
xLeftTibia = min([p8(1), geo.tibiaLowerCenter(1), -0.10]) - 0.15;
[polyT,candT] = verticalWallBoundary( ...
    xLeftTibia, geo.tibiaWallX, geo.tibiaWallY, nWall);
obs(end+1) = makeTibiaObstacle("tibiaWall", polyT, candT, p8(3), ...
    T_Pam_i, T_ICR_t1_i); %#ok<AGROW>

%% Tibia upper expanded circle, native t1
[polyT,candT] = circleBoundary( ...
    geo.tibiaUpperCenter, geo.tibiaUpperClearRadius, nCircle);
obs(end+1) = makeTibiaObstacle("tibiaUpper", polyT, candT, p8(3), ...
    T_Pam_i, T_ICR_t1_i); %#ok<AGROW>

end


function out = makeFemurObstacle(name, polyXY, candXY, z)
out.name = string(name);
out.frame = "femur";
out.polyXY = polyXY;
out.candXY = candXY;
out.candNativeXYZ = [candXY, repmat(z,size(candXY,1),1)];
end


function out = makeTibiaObstacle(name, polyNative, candNative, z, T_Pam_i, T_ICR_t1_i)
out.name = string(name);
out.frame = "t1";
out.polyXY = transformT1XYtoFemur(polyNative, z, T_Pam_i, T_ICR_t1_i);
out.candXY = transformT1XYtoFemur(candNative, z, T_Pam_i, T_ICR_t1_i);
out.candNativeXYZ = [candNative, repmat(z,size(candNative,1),1)];
end


function xyFemur = transformT1XYtoFemur(xyT1, z, T_Pam_i, T_ICR_t1_i)
xyFemur = zeros(size(xyT1));
for k = 1:size(xyT1,1)
    qT1 = [xyT1(k,:), z];
    qICR = RowVecTrans(T_ICR_t1_i, qT1);
    qFemur = RowVecTrans(T_Pam_i, qICR);
    xyFemur(k,:) = qFemur(1:2);
end
end


function [nodes,G] = buildVisibilityGraph( ...
    p1, p8, p8_fem, obs, nEdgeProbe, insideTol)
% Node 1 is p1 and node 2 is p8 so shortestpath(G,1,2) is always valid.

xy = [p1(1:2); p8_fem(1:2)];
frame = ["femur"; "t1"];
obstacle = ["p1"; "p8"];
nativeXYZ = [p1; p8];

for m = 1:numel(obs)
    C = obs(m).candXY;
    Cnative = obs(m).candNativeXYZ;

    for k = 1:size(C,1)
        if candidateInsideOtherObstacle(C(k,:), m, obs, insideTol)
            continue
        end

        xy(end+1,:) = C(k,:); %#ok<AGROW>
        frame(end+1,1) = obs(m).frame; %#ok<AGROW>
        obstacle(end+1,1) = obs(m).name; %#ok<AGROW>
        nativeXYZ(end+1,:) = Cnative(k,:); %#ok<AGROW>
    end
end

% Remove exact/near duplicate global XY candidates while preserving p1/p8.
keep = true(size(xy,1),1);
for a = 3:size(xy,1)
    if ~keep(a)
        continue
    end
    d = vecnorm(xy(1:a-1,:) - xy(a,:),2,2);
    if any(d < 1e-9)
        keep(a) = false;
    end
end

xy = xy(keep,:);
frame = frame(keep);
obstacle = obstacle(keep);
nativeXYZ = nativeXYZ(keep,:);

nodes.xy = xy;
nodes.frame = frame;
nodes.obstacle = obstacle;
nodes.nativeXYZ = nativeXYZ;

n = size(xy,1);
I = zeros(n*(n-1)/2,1);
J = zeros(n*(n-1)/2,1);
W = zeros(n*(n-1)/2,1);
ne = 0;

for a = 1:n-1
    A = xy(a,:);
    for b = a+1:n
        B = xy(b,:);

        if norm(B-A) < 1e-12
            continue
        end

        if segmentClear(A,B,obs,nEdgeProbe,insideTol)
            ne = ne + 1;
            I(ne) = a;
            J(ne) = b;
            W(ne) = norm(B-A);
        end
    end
end

I = I(1:ne);
J = J(1:ne);
W = W(1:ne);

G = graph(I,J,W,n);
end


function tf = candidateInsideOtherObstacle(P, ownIndex, obs, insideTol)
tf = false;
for m = 1:numel(obs)
    if m == ownIndex
        continue
    end

    [in,on] = inpolygon(P(1),P(2),obs(m).polyXY(:,1),obs(m).polyXY(:,2));
    if in && ~on
        tf = true;
        return
    end

    %#ok<NASGU> insideTol is intentionally retained as a named diagnostic
    % tolerance placeholder if we later replace polygon membership with an
    % analytic signed-distance check.
end
end


function clearFlag = segmentClear(A,B,obs,nProbe,insideTol)
% A graph edge is allowed when its interior does not enter any obstacle.
% Endpoints may lie on obstacle boundaries.

t = linspace(0,1,nProbe+2).';
t = t(2:end-1);
P = A + t.*(B-A);

clearFlag = true;
for m = 1:numel(obs)
    [in,on] = inpolygon(P(:,1),P(:,2),obs(m).polyXY(:,1),obs(m).polyXY(:,2));

    % inpolygon treats boundary points as in=true,on=true. Reject only
    % points strictly inside the polygon.
    if any(in & ~on)
        clearFlag = false;
        return
    end
end

% Retain the argument explicitly for easy future tightening.
if insideTol < 0 %#ok<UNRCH>
    clearFlag = false;
end
end


function nodePath = simplifyDynamicPath(nodePath,nodes,obs,nEdgeProbe,insideTol)
% Remove intermediate graph vertices whenever the surrounding two points
% have direct line of sight. Repeat until no more points can be removed.

changed = true;
while changed && numel(nodePath) > 2
    changed = false;
    k = 2;

    while k <= numel(nodePath)-1
        A = nodes.xy(nodePath(k-1),:);
        B = nodes.xy(nodePath(k+1),:);

        if segmentClear(A,B,obs,nEdgeProbe,insideTol)
            nodePath(k) = [];
            changed = true;
        else
            k = k + 1;
        end
    end
end
end


function [poly,cand] = circleBoundary(C,R,n)
th = linspace(0,2*pi,n+1).';
th(end) = [];
poly = C + R*[cos(th), sin(th)];
cand = poly;
end


function [poly,cand] = ellipseBoundary(geo,n)
th = linspace(0,2*pi,n+1).';
th(end) = [];
local = [geo.femurEllipseA*cos(th), geo.femurEllipseB*sin(th)];
Rell = [cos(geo.femurEllipseTheta), -sin(geo.femurEllipseTheta); ...
        sin(geo.femurEllipseTheta),  cos(geo.femurEllipseTheta)];
poly = geo.femurProfileCenter + local*Rell';
cand = poly;
end


function [poly,cand] = verticalWallBoundary(xLeft,xRight,yRange,nWall)
y0 = min(yRange);
y1 = max(yRange);

% Closed rectangle used as the blocked body. The exposed clearance surface
% is the x=xRight side.
poly = [ ...
    xRight, y0; ...
    xRight, y1; ...
    xLeft,  y1; ...
    xLeft,  y0];

y = linspace(y0,y1,nWall).';
cand = [repmat(xRight,nWall,1), y];
end
