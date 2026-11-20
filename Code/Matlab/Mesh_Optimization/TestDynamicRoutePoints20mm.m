%% TestDynamicRoutePoints20mm
% Standalone diagnostic for the 20 mm knee-extensor route.
%
% PURPOSE
%   Let the actual clearance geometry decide how many physical contact
%   REGIONS are needed in the femur frame and in the t1 frame.
%
% IMPORTANT
%   - Fully inflated 20 mm BPA diameter = 38.5 mm.
%   - BPA clearance radius = 19.25 mm.
%   - Collision tests are analytic for circles and the rotated ellipse.
%   - The straight wall regions use exact segment/polygon intersection.
%   - Curved obstacles are followed with boundary-arc graph edges.
%   - Consecutive graph vertices on the SAME obstacle are clustered into
%     ONE representative physical contact region.
%   - This script DOES NOT modify buildDistalRingLocation20mm.m, the
%     context builder, predictor, or optimizer.
%
% OUTPUTS LEFT IN WORKSPACE
%   routeDynamic
%   routeDynamicSummary
%
% If xBest exists in the workspace, the script tests xBest.
% Otherwise it tests ctx.x0.

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

%% Corrected 20 mm BPA clearance geometry
bpaInflatedDiameter = 0.0385;
bpaClearRadius = bpaInflatedDiameter/2;       % 0.01925 m

% Local copy only. Do not alter ctx.geo.
geo = ctx.geo;
geo.bpaRadius = bpaClearRadius;

% Rebuild expanded clearances from the physical obstacle dimensions.
geo.femurCylClearRadius     = 0.0100 + bpaClearRadius;
geo.femurFilletClearRadius  = 0.0200 + bpaClearRadius;

geo.femurEllipseA = 0.06640/2 + bpaClearRadius;
geo.femurEllipseB = 0.05439/2 + bpaClearRadius;

geo.tibiaLowerClearRadius = 0.0050 + bpaClearRadius;
geo.tibiaUpperClearRadius = 0.0075 + bpaClearRadius;

% Original common tangent x = 52.36 mm.
geo.tibiaWallX = 0.05236 + bpaClearRadius;

% Original femur exposed line = profile center + 30 mm.
geo.femurLineX = geo.femurProfileCenter(1) + 0.030 + bpaClearRadius;

fprintf('\n============================================================\n')
fprintf('DYNAMIC ROUTE-POINT TEST -- ANALYTIC COLLISION VERSION\n')
fprintf('Design: %s\n', designName)
fprintf('p1 = [%.6f %.6f %.6f] m, femur frame\n', p1)
fprintf('p8 = [%.6f %.6f %.6f] m, t1 frame\n', p8)
fprintf('Inflated BPA diameter = %.3f mm\n', 1000*bpaInflatedDiameter)
fprintf('BPA clearance radius  = %.3f mm\n', 1000*bpaClearRadius)
fprintf('============================================================\n')

%% Boundary resolution
% These points are only used to provide possible tangent/wrap vertices.
% Collision detection itself is analytic and does NOT depend on this
% sampling density.
nCircle  = 24;
nEllipse = 36;
nWall    = 7;

geomTol = 1e-9;

plotAnglesD = [10 0 -10 -30 -60 -90 -108 -120];

% Set false if you want the diagnostic to run with minimum extra work.
overlayCurrentBuilder = true;

%% Optional current-builder overlay
predCurrent = [];

if overlayCurrentBuilder
    try
        predCurrent = predictKneeExt20mm(xRouteTest,ctx);
        if ~predCurrent.ok
            predCurrent = [];
        end
    catch
        predCurrent = [];
    end
end

%% Allocate outputs
N = ctx.N;

routeDynamic = repmat(struct( ...
    'phiD',NaN, ...
    'ok',false, ...
    'rawPointsFemur',[], ...
    'rawPointFrame',strings(0,1), ...
    'rawPointObstacle',strings(0,1), ...
    'pointsFemur',[], ...
    'pointsNative',[], ...
    'pointFrame',strings(0,1), ...
    'pointObstacle',strings(0,1), ...
    'pathLength',NaN, ...
    'nFemurContacts',NaN, ...
    'nTibiaContacts',NaN, ...
    'nGeneratedContacts',NaN, ...
    'CrossPoint',NaN, ...
    'frameSwitches',NaN, ...
    'obstacles',[]), N,1);

okRoute = false(N,1);
nFemur = nan(N,1);
nTibia = nan(N,1);
nGenerated = nan(N,1);
crossPointDynamic = nan(N,1);
pathLength = nan(N,1);
frameSwitches = nan(N,1);

fprintf('\nSolving %d knee positions',N)

for i = 1:N

    if i == 1 || mod(i,10) == 0 || i == N
        fprintf('\n  position %3d / %3d, phi = %8.3f deg', ...
            i,N,ctx.phiD(i))
    end

    obs = buildAnalyticObstacles( ...
        geo,p1,p8, ...
        ctx.T_Pam(:,:,i), ...
        ctx.T_ICR_t1(:,:,i), ...
        nCircle,nEllipse,nWall);

    % p8 expressed in the femur frame.
    p8_icr = RowVecTrans(ctx.T_ICR_t1(:,:,i),p8);
    p8_fem = RowVecTrans(ctx.T_Pam(:,:,i),p8_icr);

    [nodes,G] = buildVisibilityGraphAnalytic( ...
        p1,p8,p8_fem,obs,geomTol);

    [nodePath,distPath] = shortestpath(G,1,2);

    routeDynamic(i).phiD = ctx.phiD(i);
    routeDynamic(i).obstacles = obs;

    if isempty(nodePath) || ~isfinite(distPath)
        continue
    end

    % Remove unnecessary line-of-sight vertices, but preserve the sampled
    % boundary vertices that are needed to follow a curved obstacle.
    nodePath = simplifyDynamicPathAnalytic( ...
        nodePath,nodes,obs,geomTol);

    rawXY = nodes.xy(nodePath,:);
    rawNative = nodes.nativeXYZ(nodePath,:);
    rawFrames = nodes.frame(nodePath);
    rawNames = nodes.obstacle(nodePath);

    rawP3 = routePointsToFemur3D( ...
        rawXY,rawNative,rawFrames,p1, ...
        ctx.T_Pam(:,:,i),ctx.T_ICR_t1(:,:,i));

    % Graph distance includes curved boundary-arc weights.
    L = distPath;

    % Collapse each consecutive run on the same obstacle to ONE reported
    % physical contact vertex.
    [clusterIdx,clusterFrames,clusterNames] = ...
        clusterSameObstacleRuns(nodePath,nodes);

    clusterXY = nodes.xy(clusterIdx,:);
    clusterNative = nodes.nativeXYZ(clusterIdx,:);

    clusterP3 = routePointsToFemur3D( ...
        clusterXY,clusterNative,clusterFrames,p1, ...
        ctx.T_Pam(:,:,i),ctx.T_ICR_t1(:,:,i));

    isEndpoint = clusterNames == "p1" | clusterNames == "p8";

    nf = sum(clusterFrames == "femur" & ~isEndpoint);
    nt = sum(clusterFrames == "t1" & ~isEndpoint);
    ng = nf + nt;

    fs = sum(clusterFrames(2:end) ~= clusterFrames(1:end-1));

    firstT1 = find(clusterFrames == "t1",1,'first');

    if isempty(firstT1)
        cp = NaN;
    else
        cp = firstT1;
    end

    routeDynamic(i).ok = true;

    routeDynamic(i).rawPointsFemur = rawP3;
    routeDynamic(i).rawPointFrame = rawFrames;
    routeDynamic(i).rawPointObstacle = rawNames;

    routeDynamic(i).pointsFemur = clusterP3;
    routeDynamic(i).pointsNative = clusterNative;
    routeDynamic(i).pointFrame = clusterFrames;
    routeDynamic(i).pointObstacle = clusterNames;

    routeDynamic(i).pathLength = L;
    routeDynamic(i).nFemurContacts = nf;
    routeDynamic(i).nTibiaContacts = nt;
    routeDynamic(i).nGeneratedContacts = ng;
    routeDynamic(i).CrossPoint = cp;
    routeDynamic(i).frameSwitches = fs;

    okRoute(i) = true;
    nFemur(i) = nf;
    nTibia(i) = nt;
    nGenerated(i) = ng;
    crossPointDynamic(i) = cp;
    pathLength(i) = L;
    frameSwitches(i) = fs;
end

fprintf('\nDone.\n')

routeDynamicSummary = table( ...
    ctx.phiD(:),okRoute,nFemur,nTibia,nGenerated, ...
    crossPointDynamic,pathLength,frameSwitches, ...
    'VariableNames',{ ...
    'phiD','ok','FemurContacts','TibiaContacts','GeneratedContacts', ...
    'CrossPoint','PathLength_m','FrameSwitches'});

%% Print contact-count changes from extension toward flexion
fprintf('\n========== CLUSTERED CONTACT-COUNT CHANGES ==========\n')

lastState = [NaN NaN NaN NaN];

for i = N:-1:1

    if ~okRoute(i)
        fprintf('phi = %9.4f deg : NO COLLISION-FREE PATH FOUND\n',ctx.phiD(i))
        continue
    end

    state = [ ...
        nFemur(i), ...
        nTibia(i), ...
        crossPointDynamic(i), ...
        frameSwitches(i)];

    if any(state ~= lastState)
        fprintf(['phi = %9.4f deg : femur = %2d, t1 = %2d, ' ...
                 'total = %2d, CrossPoint = %g, frame switches = %d\n'], ...
            ctx.phiD(i),nFemur(i),nTibia(i),nGenerated(i), ...
            crossPointDynamic(i),frameSwitches(i));

        lastState = state;
    end
end

fprintf('=====================================================\n')

%% Explicitly print the region that failed in the first test
fprintf('\n========== -105 TO -120 DEG CHECK ==========\n')

idxTail = find(ctx.phiD <= -105);

for ii = idxTail(:)'
    fprintf(['phi = %9.4f deg | ok=%d | F=%g | T=%g | total=%g | ' ...
             'CrossPoint=%g | L=%.6f m\n'], ...
        ctx.phiD(ii),okRoute(ii),nFemur(ii),nTibia(ii), ...
        nGenerated(ii),crossPointDynamic(ii),pathLength(ii));
end

fprintf('============================================\n')

%% Contact-count plot
figure('Name','Dynamic clustered route point count','Color','w')
hold on

plot(ctx.phiD,nFemur,'o-','LineWidth',1.5)
plot(ctx.phiD,nTibia,'s-','LineWidth',1.5)
plot(ctx.phiD,nGenerated,'^-','LineWidth',1.5)

grid on
xlabel('Knee angle, deg')
ylabel('Clustered physical contact regions')
legend('Femur-frame contacts','t1-frame contacts','Total contacts', ...
    'Location','best')
title('Geometry-selected route complexity')

%% Dynamic CrossPoint plot
figure('Name','Dynamic CrossPoint','Color','w')

plot(ctx.phiD,crossPointDynamic,'o-','LineWidth',1.5)

grid on
xlabel('Knee angle, deg')
ylabel('CrossPoint = first t1 row')
title('Variable-length route CrossPoint')

%% Route snapshots
figure('Name','Dynamic route snapshots','Color','w')
tiledlayout(2,4,'TileSpacing','compact','Padding','compact')

for q = 1:numel(plotAnglesD)

    [~,i] = min(abs(ctx.phiD-plotAnglesD(q)));

    nexttile
    hold on

    obs = routeDynamic(i).obstacles;

    for m = 1:numel(obs)

        xy = obs(m).plotXY;

        plot([xy(:,1);xy(1,1)], ...
             [xy(:,2);xy(1,2)], ...
             '-','LineWidth',1,'HandleVisibility','off')
    end

    if routeDynamic(i).ok

        % Collision-free raw route. Multiple adjacent points on a curved
        % obstacle approximate one wrap/contact region.
        Pr = routeDynamic(i).rawPointsFemur;

        plot(Pr(:,1),Pr(:,2),'-','LineWidth',1.1, ...
            'HandleVisibility','off')

        % Reported clustered contact vertices. Plot markers only because a
        % straight line between cluster representatives could cut through
        % the wrapped obstacle.
        P = routeDynamic(i).pointsFemur;

        plot(P(:,1),P(:,2),'o','MarkerSize',6,'LineWidth',1.5, ...
            'HandleVisibility','off')

        for k = 1:size(P,1)

            text(P(k,1),P(k,2), ...
                sprintf(' %d:%s/%s', ...
                k,routeDynamic(i).pointFrame(k), ...
                routeDynamic(i).pointObstacle(k)), ...
                'FontSize',7,'Interpreter','none')
        end
    end

    if ~isempty(predCurrent)

        Pc = predCurrent.Location(:,:,i);

        for j = 5:8
            Pc(j,:) = RowVecTrans(ctx.T_Pam(:,:,i),Pc(j,:));
        end

        plot(Pc(:,1),Pc(:,2),'--','LineWidth',1.1, ...
            'HandleVisibility','off')
    end

    axis equal
    grid on

    xlabel('Femur-frame x, m')
    ylabel('Femur-frame y, m')

    title(sprintf('\\phi = %.2f deg | F=%g, T=%g, CP=%g', ...
        ctx.phiD(i),nFemur(i),nTibia(i),crossPointDynamic(i)))
end

fprintf('\nFull results are in routeDynamic and routeDynamicSummary.\n')
fprintf('This script has NOT changed buildDistalRingLocation20mm.m or the optimizer.\n')

%% =====================================================================
%% Local functions
%% =====================================================================

function obs = buildAnalyticObstacles( ...
    geo,p1,p8,T_Pam_i,T_ICR_t1_i,nCircle,nEllipse,nWall)

obs = struct( ...
    'name',{}, ...
    'frame',{}, ...
    'type',{}, ...
    'plotXY',{}, ...
    'candXY',{}, ...
    'candNativeXYZ',{}, ...
    'cyclic',{}, ...
    'center',{}, ...
    'radius',{}, ...
    'ellipseA',{}, ...
    'ellipseB',{}, ...
    'ellipseTheta',{}, ...
    'polyXY',{});

%% Femur cylinder
[plotXY,candXY] = circleBoundary( ...
    geo.femurCylCenter,geo.femurCylClearRadius,nCircle);

obs(end+1) = makeCircleObstacle( ...
    "femurCylinder","femur", ...
    geo.femurCylCenter,geo.femurCylClearRadius, ...
    plotXY,candXY,p1(3)); %#ok<AGROW>

%% Femur straight-body clearance
xLeftFemur = min([p1(1),geo.femurProfileCenter(1),-0.10])-0.15;

[polyXY,candXY] = verticalWallBoundary( ...
    xLeftFemur,geo.femurLineX,geo.femurLineY,nWall);

obs(end+1) = makePolygonObstacle( ...
    "femurWall","femur", ...
    polyXY,candXY,p1(3)); %#ok<AGROW>

%% Femoral condyle ellipse
[plotXY,candXY] = ellipseBoundary(geo,nEllipse);

obs(end+1) = makeEllipseObstacle( ...
    "femurEllipse","femur", ...
    geo.femurProfileCenter, ...
    geo.femurEllipseA,geo.femurEllipseB,geo.femurEllipseTheta, ...
    plotXY,candXY,p1(3)); %#ok<AGROW>

%% Tibia lower circle -- native t1
[plotT,candT] = circleBoundary( ...
    geo.tibiaLowerCenter,geo.tibiaLowerClearRadius,nCircle);

plotF = transformT1XYtoFemur( ...
    plotT,p8(3),T_Pam_i,T_ICR_t1_i);

candF = transformT1XYtoFemur( ...
    candT,p8(3),T_Pam_i,T_ICR_t1_i);

centerF = transformT1XYtoFemur( ...
    geo.tibiaLowerCenter,p8(3),T_Pam_i,T_ICR_t1_i);

obs(end+1) = makeCircleObstacleT1( ...
    "tibiaLower", ...
    centerF,geo.tibiaLowerClearRadius, ...
    plotF,candF,candT,p8(3)); %#ok<AGROW>

%% Tibia straight-body clearance -- native t1
xLeftTibia = min([p8(1),geo.tibiaLowerCenter(1),-0.10])-0.15;

[polyT,candT] = verticalWallBoundary( ...
    xLeftTibia,geo.tibiaWallX,geo.tibiaWallY,nWall);

polyF = transformT1XYtoFemur( ...
    polyT,p8(3),T_Pam_i,T_ICR_t1_i);

candF = transformT1XYtoFemur( ...
    candT,p8(3),T_Pam_i,T_ICR_t1_i);

obs(end+1) = makePolygonObstacleT1( ...
    "tibiaWall",polyF,candF,candT,p8(3)); %#ok<AGROW>

%% Tibia upper circle -- native t1
[plotT,candT] = circleBoundary( ...
    geo.tibiaUpperCenter,geo.tibiaUpperClearRadius,nCircle);

plotF = transformT1XYtoFemur( ...
    plotT,p8(3),T_Pam_i,T_ICR_t1_i);

candF = transformT1XYtoFemur( ...
    candT,p8(3),T_Pam_i,T_ICR_t1_i);

centerF = transformT1XYtoFemur( ...
    geo.tibiaUpperCenter,p8(3),T_Pam_i,T_ICR_t1_i);

obs(end+1) = makeCircleObstacleT1( ...
    "tibiaUpper", ...
    centerF,geo.tibiaUpperClearRadius, ...
    plotF,candF,candT,p8(3)); %#ok<AGROW>

end


function out = emptyObstacle
out = struct( ...
    'name',"", ...
    'frame',"", ...
    'type',"", ...
    'plotXY',[], ...
    'candXY',[], ...
    'candNativeXYZ',[], ...
    'cyclic',false, ...
    'center',[NaN NaN], ...
    'radius',NaN, ...
    'ellipseA',NaN, ...
    'ellipseB',NaN, ...
    'ellipseTheta',NaN, ...
    'polyXY',[]);
end


function out = makeCircleObstacle( ...
    name,frame,center,radius,plotXY,candXY,z)

out = emptyObstacle;

out.name = string(name);
out.frame = string(frame);
out.type = "circle";

out.plotXY = plotXY;
out.candXY = candXY;
out.candNativeXYZ = [candXY,repmat(z,size(candXY,1),1)];

out.cyclic = true;
out.center = center;
out.radius = radius;

end


function out = makeCircleObstacleT1( ...
    name,centerF,radius,plotF,candF,candT,z)

out = emptyObstacle;

out.name = string(name);
out.frame = "t1";
out.type = "circle";

out.plotXY = plotF;
out.candXY = candF;
out.candNativeXYZ = [candT,repmat(z,size(candT,1),1)];

out.cyclic = true;
out.center = centerF;
out.radius = radius;

end


function out = makeEllipseObstacle( ...
    name,frame,center,a,b,theta,plotXY,candXY,z)

out = emptyObstacle;

out.name = string(name);
out.frame = string(frame);
out.type = "ellipse";

out.plotXY = plotXY;
out.candXY = candXY;
out.candNativeXYZ = [candXY,repmat(z,size(candXY,1),1)];

out.cyclic = true;
out.center = center;
out.ellipseA = a;
out.ellipseB = b;
out.ellipseTheta = theta;

end


function out = makePolygonObstacle( ...
    name,frame,polyXY,candXY,z)

out = emptyObstacle;

out.name = string(name);
out.frame = string(frame);
out.type = "polygon";

out.plotXY = polyXY;
out.polyXY = polyXY;
out.candXY = candXY;
out.candNativeXYZ = [candXY,repmat(z,size(candXY,1),1)];

out.cyclic = false;

end


function out = makePolygonObstacleT1( ...
    name,polyF,candF,candT,z)

out = emptyObstacle;

out.name = string(name);
out.frame = "t1";
out.type = "polygon";

out.plotXY = polyF;
out.polyXY = polyF;
out.candXY = candF;
out.candNativeXYZ = [candT,repmat(z,size(candT,1),1)];

out.cyclic = false;

end


function xyFemur = transformT1XYtoFemur( ...
    xyT1,z,T_Pam_i,T_ICR_t1_i)

if isvector(xyT1)
    xyT1 = reshape(xyT1,1,2);
end

xyFemur = zeros(size(xyT1));

for k = 1:size(xyT1,1)

    qT1 = [xyT1(k,:),z];
    qICR = RowVecTrans(T_ICR_t1_i,qT1);
    qFemur = RowVecTrans(T_Pam_i,qICR);

    xyFemur(k,:) = qFemur(1:2);
end

end


function [nodes,G] = buildVisibilityGraphAnalytic( ...
    p1,p8,p8_fem,obs,tol)

% Node 1 = p1, Node 2 = p8.
xy = [p1(1:2);p8_fem(1:2)];
frame = ["femur";"t1"];
obstacle = ["p1";"p8"];
nativeXYZ = [p1;p8];

obsIndex = [0;0];
candIndex = [0;0];

% Add exposed boundary candidates that are not inside another obstacle.
for m = 1:numel(obs)

    C = obs(m).candXY;
    Cnative = obs(m).candNativeXYZ;

    for k = 1:size(C,1)

        if pointInsideOtherObstacle(C(k,:),m,obs,tol)
            continue
        end

        xy(end+1,:) = C(k,:); %#ok<AGROW>
        frame(end+1,1) = obs(m).frame; %#ok<AGROW>
        obstacle(end+1,1) = obs(m).name; %#ok<AGROW>
        nativeXYZ(end+1,:) = Cnative(k,:); %#ok<AGROW>

        obsIndex(end+1,1) = m; %#ok<AGROW>
        candIndex(end+1,1) = k; %#ok<AGROW>
    end
end

% Remove duplicate global XY candidates, preserving p1/p8.
keep = true(size(xy,1),1);

for a = 3:size(xy,1)

    if ~keep(a)
        continue
    end

    d = vecnorm(xy(1:a-1,:)-xy(a,:),2,2);

    if any(d < 1e-9)
        keep(a) = false;
    end
end

xy = xy(keep,:);
frame = frame(keep);
obstacle = obstacle(keep);
nativeXYZ = nativeXYZ(keep,:);
obsIndex = obsIndex(keep);
candIndex = candIndex(keep);

nodes.xy = xy;
nodes.frame = frame;
nodes.obstacle = obstacle;
nodes.nativeXYZ = nativeXYZ;
nodes.obsIndex = obsIndex;
nodes.candIndex = candIndex;

n = size(xy,1);

% Conservative upper bound: all visibility pairs plus boundary edges.
maxEdges = n*(n-1)/2 + 2*n;

I = zeros(maxEdges,1);
J = zeros(maxEdges,1);
W = zeros(maxEdges,1);
ne = 0;

%% Explicit boundary-following edges
% A true circle/ellipse cannot be followed by a straight chord because that
% chord lies inside the body. These graph edges represent following the
% curved clearance boundary between neighboring sampled vertices.

for m = 1:numel(obs)

    ids = find(obsIndex == m);

    if numel(ids) < 2
        continue
    end

    % Sort by original candidate order.
    [~,ord] = sort(candIndex(ids));
    ids = ids(ord);

    originalK = candIndex(ids);

    for q = 1:numel(ids)-1

        % Do not jump across candidates removed because they were hidden
        % inside another obstacle.
        if originalK(q+1)-originalK(q) ~= 1
            continue
        end

        ne = ne+1;
        I(ne) = ids(q);
        J(ne) = ids(q+1);
        W(ne) = boundaryStepLength( ...
            obs(m),originalK(q),originalK(q+1));
    end

    % Cyclic wrap from last candidate back to first.
    if obs(m).cyclic

        nCand = size(obs(m).candXY,1);

        if originalK(1) == 1 && originalK(end) == nCand

            ne = ne+1;
            I(ne) = ids(end);
            J(ne) = ids(1);
            W(ne) = boundaryStepLength(obs(m),nCand,1);
        end
    end
end

%% Straight visibility edges
for a = 1:n-1

    A = xy(a,:);

    for b = a+1:n

        % Boundary-neighbor pairs are already represented by the physical
        % boundary-following edge above.
        if obsIndex(a) > 0 && obsIndex(a) == obsIndex(b)

            ka = candIndex(a);
            kb = candIndex(b);
            nCand = size(obs(obsIndex(a)).candXY,1);

            adjacent = abs(ka-kb) == 1;

            if obs(obsIndex(a)).cyclic
                adjacent = adjacent || abs(ka-kb) == nCand-1;
            end

            if adjacent
                continue
            end
        end

        B = xy(b,:);

        if norm(B-A) < 1e-12
            continue
        end

        if segmentClearAnalytic(A,B,obs,tol)

            ne = ne+1;
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


function w = boundaryStepLength(ob,k1,k2)

P1 = ob.candXY(k1,:);
P2 = ob.candXY(k2,:);

switch ob.type

    case "circle"

        n = size(ob.candXY,1);

        % Adjacent samples are equally spaced in angle.
        dtheta = 2*pi/n;

        w = ob.radius*dtheta;

    case "ellipse"

        % Small adjacent ellipse arc approximated by the sampled chord.
        % At nEllipse = 36 this error is very small for this diagnostic.
        w = norm(P2-P1);

    otherwise

        w = norm(P2-P1);
end

end


function tf = pointInsideOtherObstacle(P,ownIndex,obs,tol)

tf = false;

for m = 1:numel(obs)

    if m == ownIndex
        continue
    end

    if pointStrictlyInsideObstacle(P,obs(m),tol)
        tf = true;
        return
    end
end

end


function inside = pointStrictlyInsideObstacle(P,ob,tol)

switch ob.type

    case "circle"

        inside = norm(P-ob.center) < ob.radius-tol;

    case "ellipse"

        local = ellipseLocalNormalized(P,ob);

        inside = dot(local,local) < 1-tol;

    case "polygon"

        [in,on] = inpolygon( ...
            P(1),P(2),ob.polyXY(:,1),ob.polyXY(:,2));

        inside = in && ~on;

    otherwise

        inside = false;
end

end


function clearFlag = segmentClearAnalytic(A,B,obs,tol)

clearFlag = true;

for m = 1:numel(obs)

    switch obs(m).type

        case "circle"

            if segmentPenetratesCircle( ...
                    A,B,obs(m).center,obs(m).radius,tol)

                clearFlag = false;
                return
            end

        case "ellipse"

            if segmentPenetratesEllipse(A,B,obs(m),tol)

                clearFlag = false;
                return
            end

        case "polygon"

            if segmentPenetratesPolygon(A,B,obs(m).polyXY,tol)

                clearFlag = false;
                return
            end
    end
end

end


function hit = segmentPenetratesCircle(A,B,C,R,tol)
% True only when the OPEN segment interior enters the circle interior.
% Tangency and a boundary endpoint leaving outward are allowed.

d = B-A;
dd = dot(d,d);

if dd < 1e-20
    hit = norm(A-C) < R-tol;
    return
end

t = dot(C-A,d)/dd;
t = min(max(t,0),1);

Q = A+t*d;

hit = norm(Q-C) < R-tol;

end


function hit = segmentPenetratesEllipse(A,B,ob,tol)
% Transform the rotated ellipse to a unit circle, then use the exact
% segment-to-origin minimum-distance test.

A0 = ellipseLocalNormalized(A,ob);
B0 = ellipseLocalNormalized(B,ob);

d = B0-A0;
dd = dot(d,d);

if dd < 1e-20
    hit = norm(A0) < 1-tol;
    return
end

t = -dot(A0,d)/dd;
t = min(max(t,0),1);

Q = A0+t*d;

hit = norm(Q) < 1-tol;

end


function local = ellipseLocalNormalized(P,ob)

th = ob.ellipseTheta;

R = [ ...
    cos(th), -sin(th); ...
    sin(th),  cos(th)];

q = (P-ob.center)*R;

local = [ ...
    q(1)/ob.ellipseA, ...
    q(2)/ob.ellipseB];

end


function hit = segmentPenetratesPolygon(A,B,poly,tol)
% Exact straight-segment / closed-polygon interior test.
% Polygon boundaries are allowed; only strict interior penetration rejects
% the graph edge.

% Endpoint strictly inside => penetration.
[inA,onA] = inpolygon(A(1),A(2),poly(:,1),poly(:,2));
[inB,onB] = inpolygon(B(1),B(2),poly(:,1),poly(:,2));

if (inA && ~onA) || (inB && ~onB)
    hit = true;
    return
end

% Collect all segment-edge crossing parameters.
tCuts = [0;1];

nP = size(poly,1);

for e = 1:nP

    C = poly(e,:);
    D = poly(mod(e,nP)+1,:);

    tHit = segmentSegmentIntersectionParams(A,B,C,D,tol);

    if ~isempty(tHit)
        tCuts = [tCuts;tHit(:)]; %#ok<AGROW>
    end
end

tCuts = sort(unique(round(tCuts,12)));

hit = false;

for k = 1:numel(tCuts)-1

    ta = tCuts(k);
    tb = tCuts(k+1);

    if tb-ta < 1e-12
        continue
    end

    tm = 0.5*(ta+tb);

    Pm = A+tm*(B-A);

    [in,on] = inpolygon( ...
        Pm(1),Pm(2),poly(:,1),poly(:,2));

    if in && ~on
        hit = true;
        return
    end
end

end


function tVals = segmentSegmentIntersectionParams(A,B,C,D,tol)

r = B-A;
s = D-C;

rxs = cross2(r,s);
qmp = C-A;
qmpxr = cross2(qmp,r);

tVals = [];

if abs(rxs) <= tol && abs(qmpxr) <= tol

    rr = dot(r,r);

    if rr < 1e-20
        return
    end

    t0 = dot(C-A,r)/rr;
    t1 = dot(D-A,r)/rr;

    lo = max(0,min(t0,t1));
    hi = min(1,max(t0,t1));

    if hi >= lo-tol
        tVals = [lo;hi];
    end

    return
end

if abs(rxs) <= tol
    return
end

t = cross2(qmp,s)/rxs;
u = cross2(qmp,r)/rxs;

if t >= -tol && t <= 1+tol && ...
        u >= -tol && u <= 1+tol

    tVals = min(max(t,0),1);
end

end


function z = cross2(a,b)
z = a(1)*b(2)-a(2)*b(1);
end


function nodePath = simplifyDynamicPathAnalytic( ...
    nodePath,nodes,obs,tol)

changed = true;

while changed && numel(nodePath) > 2

    changed = false;
    k = 2;

    while k <= numel(nodePath)-1

        idxPrev = nodePath(k-1);
        idxNow  = nodePath(k);
        idxNext = nodePath(k+1);

        % Do not collapse a vertex out of a same-obstacle boundary run.
        sameBoundaryRun = ...
            nodes.obsIndex(idxPrev) > 0 && ...
            nodes.obsIndex(idxPrev) == nodes.obsIndex(idxNow) && ...
            nodes.obsIndex(idxNow) == nodes.obsIndex(idxNext);

        if sameBoundaryRun
            k = k+1;
            continue
        end

        A = nodes.xy(idxPrev,:);
        B = nodes.xy(idxNext,:);

        if segmentClearAnalytic(A,B,obs,tol)

            nodePath(k) = [];
            changed = true;

        else

            k = k+1;
        end
    end
end

end


function [idxOut,frameOut,nameOut] = ...
    clusterSameObstacleRuns(nodePath,nodes)
% Consecutive raw graph vertices on one obstacle approximate one physical
% tangent/wrap contact region. Report one representative sampled boundary
% vertex for that whole run.

idxOut = zeros(0,1);
frameOut = strings(0,1);
nameOut = strings(0,1);

k = 1;

while k <= numel(nodePath)

    idx0 = nodePath(k);

    name0 = nodes.obstacle(idx0);
    frame0 = nodes.frame(idx0);

    if name0 == "p1" || name0 == "p8"

        idxOut(end+1,1) = idx0; %#ok<AGROW>
        frameOut(end+1,1) = frame0; %#ok<AGROW>
        nameOut(end+1,1) = name0; %#ok<AGROW>

        k = k+1;
        continue
    end

    k2 = k;

    while k2 < numel(nodePath)

        idxNext = nodePath(k2+1);

        if nodes.obstacle(idxNext) == name0 && ...
                nodes.frame(idxNext) == frame0

            k2 = k2+1;

        else

            break
        end
    end

    runIdx = nodePath(k:k2);
    runXY = nodes.xy(runIdx,:);

    % Representative vertex nearest the mean of the wrap run.
    qMean = mean(runXY,1);

    [~,jBest] = min(vecnorm(runXY-qMean,2,2));

    idxRep = runIdx(jBest);

    idxOut(end+1,1) = idxRep; %#ok<AGROW>
    frameOut(end+1,1) = frame0; %#ok<AGROW>
    nameOut(end+1,1) = name0; %#ok<AGROW>

    k = k2+1;
end

end


function P3 = routePointsToFemur3D( ...
    xyFemur,nativeXYZ,frames,p1,T_Pam_i,T_ICR_t1_i)

P3 = zeros(size(xyFemur,1),3);
P3(:,1:2) = xyFemur;

for k = 1:size(xyFemur,1)

    if frames(k) == "femur"

        P3(k,3) = p1(3);

    else

        qICR = RowVecTrans(T_ICR_t1_i,nativeXYZ(k,:));
        qFemur = RowVecTrans(T_Pam_i,qICR);

        P3(k,3) = qFemur(3);
    end
end

end


function [xy,cand] = circleBoundary(C,R,n)

th = linspace(0,2*pi,n+1).';
th(end) = [];

xy = C+R*[cos(th),sin(th)];
cand = xy;

end


function [xy,cand] = ellipseBoundary(geo,n)

th = linspace(0,2*pi,n+1).';
th(end) = [];

local = [ ...
    geo.femurEllipseA*cos(th), ...
    geo.femurEllipseB*sin(th)];

Rell = [ ...
    cos(geo.femurEllipseTheta), -sin(geo.femurEllipseTheta); ...
    sin(geo.femurEllipseTheta),  cos(geo.femurEllipseTheta)];

xy = geo.femurProfileCenter+local*Rell';
cand = xy;

end


function [poly,cand] = verticalWallBoundary( ...
    xLeft,xRight,yRange,nWall)

y0 = min(yRange);
y1 = max(yRange);

% Mathematical blocked region, not a literal rectangular part.
poly = [ ...
    xRight,y0; ...
    xRight,y1; ...
    xLeft, y1; ...
    xLeft, y0];

y = linspace(y0,y1,nWall).';

cand = [ ...
    repmat(xRight,nWall,1), ...
    y];

end
