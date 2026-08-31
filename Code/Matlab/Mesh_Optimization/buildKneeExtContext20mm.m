function ctx = buildKneeExtContext20mm()
    % all the transform setup code here



%% Joint rotation transformation matrices
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);
R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
T_Pam_inv = zeros(4, 4, positions);
t1toICR = zeros(1,3,positions);
T_t1_ICR = zeros(4, 4, positions);
T_ICR_t1 = zeros(4, 4, positions);

c = pi/180; %Convert from degrees to radians

%Knee Extension and Flexion
%Human
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
%Robot
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

%Theta1 to ICR
t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions); %RoM

%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;
phiD = phi*180/pi;  %Knee angle in degrees

for i = 1:positions
    hipToKnee = [fcn1(phi(i)), fcn2(phi(i)), 0];
    R(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T(:, :, i) = RpToTrans(R(:, :, i), hipToKnee');
    
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam(:, :, i) = [cos(phi(i)), -sin(phi(i)), 0;   %Rotation matrix for robot
                    sin(phi(i)), cos(phi(i)), 0;
                    0, 0, 1];
    
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot
    T_Pam_inv(:, :, i) = TransInv(T_Pam(:, :, i));
    
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0]; %distance from theta1 to ICR
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR(1,:,i)');    %transform from the ICR frame to theta1
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR(1,:,i)');    %transform from t1 frame to ICR
end

%% Build Optimization Context
ctx.Name       = 'Vastus Medialis Proximal Ring 20mm BPA';
ctx.CrossPoint = 6;
ctx.routeRows = 9;
ctx.femurRows = 1:5;
ctx.t1Rows = 6:9;
ctx.Dia        = 20;
ctx.BPAcount   = 2;

ctx.T          = T;
ctx.T_Pam      = T_Pam;
ctx.T_Pam_inv  = T_Pam_inv;
ctx.T_ICR_t1   = T_ICR_t1;
ctx.T_t1_ICR   = T_t1_ICR;
ctx.t1toICR    = t1toICR;

ctx.phi        = phi(:);
ctx.phiD       = phiD(:);
ctx.pos        = pos;
ctx.N          = numel(phiD);

ctx.fitting    = 0.021;
ctx.wraps      = 3;



[~, ctx.idxMaxFlex] = min(ctx.phiD);

%% OpenSim target: column 2 = knee angle, columns 3-5 = vastus muscles.
H = readmatrix('OpenSim_Vasti_Results.txt', ...
    'FileType', 'text', ...
    'NumHeaderLines', 7);

ctx.humanAngleD = H(:,2);

ctx.targetName = "vas_med_r";  % "vas_int_r", "vas_lat_r", or "vas_med_r"

switch ctx.targetName
    case "vas_int_r"
        ctx.humanTorque = H(:,3);
    case "vas_lat_r"
        ctx.humanTorque = H(:,4);
    case "vas_med_r"
        ctx.humanTorque = H(:,5);
    otherwise
        error('Unknown Vasti targetName: %s', ctx.targetName)
end

ctx.humanTorqueAbs = abs(ctx.humanTorque);
ctx.torqueScale = max(1, max(ctx.humanTorqueAbs));

%% BPA / stiffness constants
ctx.targetPressure = 620;
ctx.Pbins = 620;

ctx.KMAX = 0.255;          % KMAX = (rest - kmax)/rest at 620 kPa
ctx.maxRelStrain = 1.0;    % allow relative strain up to KMAX
ctx.minStrain = -0.03;

load minimizeExtPin10_results_20260819_2transforms_Z2.mat filtered_results xCols

pick = 1;     % whatever candidate you decided to use

g = filtered_results(pick,xCols);

ctx.Xi0 = g(1);
ctx.Xi1 = g(2);
ctx.Xi2 = g(3);
ctx.Xi3 = g(4);

ctx.wraps = 3;

%% 20 mm extensor routing geometry
% All geometry below is in meters.
%
% Femur-side geometry is in femur frame.
% Tibia-side geometry is in theta1 / t1 frame.

geo = struct;

% Fully inflated 20 mm Festo BPA diameter = 38.5 mm.
% Clearance geometry is based on the inflated outside radius.
geo.bpaRadius = 0.0385/2;      % 0.01925 m

% Ignore extremely small direction changes in Xi3 geometry.
geo.alphaTol = deg2rad(1.0);

%% Femur-side geometry

% 10 mm radius cylindrical object.
geo.femurCylCenter = [0.05390, -0.27476];
geo.femurCylRadius = 0.010;

% Required BPA centerline radius.
geo.femurCylClearRadius = geo.femurCylRadius + geo.bpaRadius;
% = 0.02925 m


% Femoral lower profile center.
geo.femurProfileCenter = [0.01817, -0.41031];

% Original vertical line is x = 30 mm relative to profile center.
% BPA centerline must lie another inflated BPA radius in +x.
geo.femurLineX = geo.femurProfileCenter(1) + 0.030 + geo.bpaRadius;

% Original line extends from local y = 0 to +48.74 mm.
geo.femurLineY = geo.femurProfileCenter(2) + [0, 0.04874];

% Physical femoral-condyle ellipse.
geo.femurEllipsePhysicalA = 0.06640/2;
geo.femurEllipsePhysicalB = 0.05439/2;
geo.femurEllipseTheta = deg2rad(25.11);

% Backward-compatible approximate expanded-axis values. The new route
% builder does NOT use these for condyle collision; it uses the true
% normal-offset boundary calculated below.
geo.femurEllipseA = geo.femurEllipsePhysicalA + geo.bpaRadius;
geo.femurEllipseB = geo.femurEllipsePhysicalB + geo.bpaRadius;

% Precompute the centerline clearance envelope:
% physical ellipse point + BPA_radius * outward unit normal.
%
% The robot does not have condyle material to the right of the femur
% vertical wall. Clip the inflated ellipse at femurLineX and close it with a
% straight vertical segment; no connector fillet/radius is modeled here.
nOffset = 180;
geo.femurOffsetTheta = linspace(-pi, pi, nOffset+1).';
geo.femurOffsetTheta(end) = [];

th = geo.femurEllipseTheta;
Rell = [cos(th), -sin(th); sin(th), cos(th)];

ct = cos(geo.femurOffsetTheta);
st = sin(geo.femurOffsetTheta);

ellipseLocal = [ ...
    geo.femurEllipsePhysicalA*ct, ...
    geo.femurEllipsePhysicalB*st];

normalLocal = [ ...
    ct/geo.femurEllipsePhysicalA, ...
    st/geo.femurEllipsePhysicalB];

normalLocal = normalLocal ./ vecnorm(normalLocal,2,2);

offsetLocal = ellipseLocal + geo.bpaRadius*normalLocal;

geo.femurOffsetBoundaryFull = ...
    geo.femurProfileCenter + offsetLocal*Rell';

[geo.femurOffsetBoundary, femurClipY] = ...
    clipClosedPolygonLeftOfVerticalLine( ...
        geo.femurOffsetBoundaryFull, geo.femurLineX);

geo.femurCondyleClipX = geo.femurLineX;
geo.femurCondyleClipY = [min(femurClipY), max(femurClipY)];

% Slightly denser closed grid used only for tangent root bracketing.
geo.femurOffsetThetaSearch = linspace(-pi, pi, 361).';

% Corrected avoidance references on the true normal-offset curve.
minorLocal = [0, -(geo.femurEllipsePhysicalB + geo.bpaRadius)];
geo.femurMinorAvoid = geo.femurProfileCenter + minorLocal*Rell';

majorLocal = [-(geo.femurEllipsePhysicalA + geo.bpaRadius), 0];
geo.femurMajorAvoid = geo.femurProfileCenter + majorLocal*Rell';


%% Tibia / t1 geometry

% Lower 5 mm radius cylinder.
geo.tibiaLowerCenter = [0.04736, -0.00387];
geo.tibiaLowerRadius = 0.005;
geo.tibiaLowerClearRadius = geo.tibiaLowerRadius + geo.bpaRadius;
% = 0.02425 m

% Upper 7.5 mm radius cylinder.
geo.tibiaUpperCenter = [0.04486, 0.02064];
geo.tibiaUpperRadius = 0.0075;
geo.tibiaUpperClearRadius = geo.tibiaUpperRadius + geo.bpaRadius;
% = 0.02675 m

% Original common tangent x = 52.36 mm.
geo.tibiaWallX = 0.05236 + geo.bpaRadius;
% = 0.07161 m

geo.tibiaWallY = [ ...
    geo.tibiaLowerCenter(2), ...
    geo.tibiaUpperCenter(2)];

% Route p8 must stay to the +X side of this sloped t1-frame reference,
% allowing only 7.5 mm of horizontal travel in the -X direction.  This is a
% solver feasibility limit, not a collision gate in the route-state logic.
geo.p8MinXLineT1 = [ ...
    34.60,  -18.87; ...
    31.38, -181.83]/1000;
geo.p8MinXMinusOffset = 7.5/1000;

% Temporary tibia wrap-loss curve through p8, p7, and p6.  The center is
% from the current CAD/debug model; the effective radius is averaged from
% the three route seed points after p8 is adjusted for the candidate pEnd.
geo.tibiaCurveCenter = [33.85, 19.27]/1000;


%% Tendon limits

geo.tendonMin = 0.025;
% No hard physical tendon maximum is stored here. tendonLimit20mm computes
% the geometry-dependent maximum directly from pEnd.


%% 20 mm extensor route initial geometry

%Origin and routing location for 20mm
% p1:p5 are femur frame. p6:p9 are t1 frame.
% These SolidWorks points are seed / branch-reference geometry.
%Original
% ctx.routeSeed = [ ...
%     0.04000,  0.03500,  0.00000; ... % p1 Origin
%     0.08315, -0.27476,  0.00000; ... % p2 BPA contacts mounting base
%     0.06536, -0.42389,  0.00000; ... % p3 femoral contact, updated
%     0.05438, -0.44085,  0.00000; ... % p4 femoral condyle contact, updated
%     0.03794, -0.45233,  0.00000; ... % p5 femoral condyle contact, updated
%     0.06377,  0.03955,  0.00000; ... % p6 Tibia contact initial
%     0.07161,  0.00765,  0.00000; ... % p7 Tibia contact, updated
%     0.06450, -0.02101,  0.00000; ... % p8 Tibia contact, updated
%     0.03406, -0.06687,  0.00000];    % p9 patellar ligament ring, distal

%Does not violate 20mm restriction
% ctx.routeSeed = [ ...
%     0.04000,  0.03500,  0.00000; ... % p1 Origin
%     0.08390, -0.27476,  0.00000; ... % p2 BPA contacts mounting base
%     0.06667, -0.42433,  0.00000; ... % p3 femoral contact, updated
%     0.05524, -0.44361,  0.00000; ... % p4 femoral condyle contact, updated
%     0.03826, -0.45301,  0.00000; ... % p5 femoral condyle contact, updated
%     0.07236,  0.03173,  0.00000; ... % p6 Tibia contact initial
%     0.07236, -0.01169,  0.00000; ... % p7 Tibia contact, updated
%     0.03406, -0.06687,  0.00000];    % p8 patellar ligament ring, distal

%Improved, keeps 9 points.
ctx.routeSeed = [ ...
    0.04000,  0.03500,  0.00000; ... % p1 Origin
    0.08390, -0.27476,  0.00000; ... % p2 BPA contacts mounting base
    0.06634, -0.42724,  0.00000; ... % p3 femoral contact, updated
    0.05152, -0.44674,  0.00000; ... % p4 femoral condyle contact, updated
    0.03826, -0.45301,  0.00000; ... % p5 femoral condyle contact, updated
    0.06985,  0.03433,  0.00000; ... % p6 Tibia contact initial
    0.07350,  0.01913,  0.00000; ... % p7 Tibia contact, updated
    0.07198, -0.01223,  0.00000; ... % p8 Tibia contact, updated
    0.03406, -0.06687,  0.00000];    % p9 patellar ligament ring, distal

ctx.geo = geo;

% p1 is femur-side BPA attachment.
p1_0 = ctx.routeSeed(1,:);

% pEnd is the movable distal tendon/ring attachment in t1 and is row p9.
pEnd_0 = ctx.routeSeed(9,:);

% Pick a reasonable initial tendon length above the 25 mm minimum,
% but make sure it does not exceed the geometry-dependent limit.
[~, tendonGeom0] = tendonLimit20mm(pEnd_0, ctx.geo);

if tendonGeom0 < ctx.geo.tendonMin
    error('Initial pEnd does not permit the minimum 25 mm tendon length.')
end

tendon0 = min(0.040, 0.95*tendonGeom0);
tendon0 = max(ctx.geo.tendonMin, tendon0);

fprintf('Initial tendon = %.6f m\n', tendon0)
fprintf('Initial geometry-dependent tendon maximum = %.6f m\n', tendonGeom0)


%% Build initial 20 mm route

[Location0, ~, route0] = buildDistalRingLocation20mm( ...
    p1_0, ...
    pEnd_0, ...
    tendon0, ...
    ctx);

Lmt0 = muscleLengthNormal(Location0, ctx.CrossPoint, ctx.T_Pam);


%% Initial resting length
%
% This is ONLY the initial design guess based on the longest geometric path.
% It does NOT force the BPA to be fully contracted at +10 deg extension.
% KMAX remains an upper strain/contraction limit, not an equality condition.
% The optimizer is free to choose rest within lb(7):ub(7).

rest0 = max(Lmt0) - 2*ctx.fitting - tendon0 - ctx.Xi0;

if rest0 <= 0
    error('Initial rest0 is non-positive. Check 20 mm route geometry.')
end


%% Optimization variables
% x = [p1(1:3), pEnd(1:3), rest, tendon]

x0 = [p1_0, pEnd_0, rest0, tendon0];

lb = x0;
ub = x0;

% p1 search region
lb(1:3) = p1_0 + [-0.020, -0.060, 0];
ub(1:3) = p1_0 + [ 0.125,  0.125,  0.060];

% pEnd search region
lb(4:6) = pEnd_0 + [-0.030,    -0.200, 0];
ub(4:6) = pEnd_0 + [0.075, 0.024,  0.030];

% Keep the distal t1 insertion y search inside the current CAD envelope.
% lb(5) = -114/1000;
% ub(5) =   24/1000;

% BPA rest-length bounds
lb(7) = 0.35;  %distance travel required for RoM divided by kmax 0.121/.255 = 0.475; Judging by Opt_sanity_Ext, we probably need like travel/(kmax/2) which would be 0.949 m. So set upper limit to 1.2 m?
ub(7) = 1.2;

% Global tendon bounds.
% The position-dependent upper bound is handled by nonlconExt20mm.
lb(8) = ctx.geo.tendonMin;
ub(8) = 0.3; %loose optimizer ceiling

ctx.x0 = x0;
ctx.lb = lb;
ctx.ub = ub;


%% Initial-design diagnostics

ctx.initialRawLocation = route0.raw(:,:,ctx.idxMaxFlex);
ctx.initialLocation = Location0;
ctx.initialLmt0 = Lmt0;
ctx.initialMaxLmt0 = max(Lmt0);
ctx.initialIdxMaxLmt0 = find(Lmt0 == max(Lmt0), 1, 'first');
ctx.initialMaxLmtAngleD = ctx.phiD(ctx.initialIdxMaxLmt0);
ctx.initialRest0 = rest0;
ctx.initialTendon0 = tendon0;
ctx.initialTendonMax = tendonGeom0;

end


%% Local helper functions

% function raw = distalRingRawLocation(phiD_i)
% % Location matrix from Knee_Extensor_10mm.
% %
% % Rows 1:4 are femur / parent frame.
% % Rows 5:8 are theta1 / tibia frame before conversion to ICR.
% %
% % Original Knee_Extensor_10mm points:
% % p1 = [0.040, 0.035, 0];                 %Origin
% % p2 = [0.0759, -0.27476, 0];             %BPA contacts mounting base
% % p3 = [0.05982, -0.37427, 0.000];        %femur channel contact, updated
% % p4 = [0.03955, -0.42183, 0.000];        %femoral condyle contact, updated
% % p5 = [0.05871, 0.025, 0];               %Tibia contact initial
% % p6 = [0.05871, 0.01228, 0];             %Tibia contact, updated
% % p7 = [0.054, -0.00612, 0];              %Tibia tendon contact, updated
% % p8 = [0.03604, -0.02844, 0];            %patellar ligament ring
% 
% p1 = [0.040,   0.035,    0.000];          %Origin
% p2 = [0.0759, -0.27476,  0.000];          %BPA contacts mounting base
% p3 = [0.05982,-0.37427,  0.000];          %femur channel contact, updated
% p4 = [0.03955,-0.42183,  0.000];          %femoral condyle contact, updated
% p5 = [0.05871, 0.025,    0.000];          %Tibia contact initial
% p6 = [0.05871, 0.01228,  0.000];          %Tibia contact, updated
% p7 = [0.054,  -0.00612,  0.000];          %Tibia tendon contact, updated
% p8 = [0.03604,-0.02844,  0.000];          %patellar ligament ring
% 
% %Set up angle limits (degrees)
% % AA = -68;
% BB = -69.8;   %add p4 when flexion reaches this value
% CC = -9.19;   %add p3 and p5 when flexion reaches this value
% DD = 10;
% 
%     if phiD_i > CC
%         raw = [p1;
%                p2;
%                p2;
%                p2;
%                p6;
%                p6;
%                p7;
%                p8];
% 
%     elseif phiD_i <= CC && phiD_i > BB
%         raw = [p1;
%                p2;
%                p3;
%                p3;
%                p5;
%                p6;
%                p7;
%                p8];
% 
%     elseif phiD_i <= BB
%         raw = [p1;
%                p2;
%                p3;
%                p4;
%                p5;
%                p6;
%                p7;
%                p8];
% 
%     else
%         raw = zeros(8,3);
%     end
% 
% end
% 
% function Location = buildDistalRingLocation(p1, p2, phiD, T_ICR_t1)
% % Build the extensor Location matrix using the Knee_Extensor_10mm route.
% %
% % p1 replaces row 1.
% % p2 replaces row 8.
% % Rows 5:8 are converted from theta1 -> ICR with T_ICR_t1(:,:,i).
% 
% N = numel(phiD);
% Location = zeros(8, 3, N);
% 
%     for i = 1:N
%         raw = distalRingRawLocation(phiD(i));
% 
%         raw(1,:) = p1;
%         raw(8,:) = p2;
% 
%         for j = 5:8
%             raw(j,:) = RowVecTrans(T_ICR_t1(:,:,i), raw(j,:));
%         end
% 
%         Location(:,:,i) = raw;
%     end
% 
% end

function Lmt = muscleLengthNormal(Location, CrossPoint, T)
% Normal undeformed musculotendon length calculation, matching
% MonoPamDataExplicit SegmentLengths/MuscleLength behavior.

N = size(Location, 3);
M = size(Location, 1);
Lmt = zeros(N, 1);

    for ii = 1:N
        for j = 1:M-1
            pointA = Location(j,:,ii);
            pointB = Location(j+1,:,ii);
    
            if j+1 == CrossPoint
                pointB = RowVecTrans(T(:,:,ii), pointB);
            end
    
            Lmt(ii) = Lmt(ii) + norm(pointA - pointB);
        end
    end

end


function [polyClip, yCuts] = clipClosedPolygonLeftOfVerticalLine(poly, xMax)
% Sutherland-Hodgman clip of a closed polygon to x <= xMax.

tol = 1e-12;
n = size(poly, 1);
polyClip = zeros(0, 2);
yCuts = zeros(0, 1);

for i = 1:n
    A = poly(i,:);
    B = poly(mod(i,n)+1,:);

    insideA = A(1) <= xMax + tol;
    insideB = B(1) <= xMax + tol;

    if insideA && insideB
        polyClip(end+1,:) = B; %#ok<AGROW>

    elseif insideA && ~insideB
        I = verticalLineIntersection(A, B, xMax);
        polyClip(end+1,:) = I; %#ok<AGROW>
        yCuts(end+1,1) = I(2); %#ok<AGROW>

    elseif ~insideA && insideB
        I = verticalLineIntersection(A, B, xMax);
        polyClip(end+1,:) = I; %#ok<AGROW>
        yCuts(end+1,1) = I(2); %#ok<AGROW>
        polyClip(end+1,:) = B; %#ok<AGROW>
    end
end

polyClip = removeAdjacentDuplicatePoints(polyClip, tol);
yCuts = uniqueWithinTolerance(sort(yCuts), 1e-9);

if size(polyClip,1) < 3 || numel(yCuts) < 2
    error('buildKneeExtContext20mm:BadFemurCondyleClip', ...
        'Unable to clip femur condyle boundary at the vertical wall.')
end

end


function I = verticalLineIntersection(A, B, xValue)

dx = B(1) - A(1);

if abs(dx) < 1e-14
    I = [xValue, A(2)];
    return
end

t = (xValue - A(1))/dx;
t = min(max(t, 0), 1);
I = A + t*(B-A);
I(1) = xValue;

end


function P = removeAdjacentDuplicatePoints(P, tol)

if isempty(P)
    return
end

keep = true(size(P,1),1);

for i = 2:size(P,1)
    if norm(P(i,:) - P(i-1,:)) <= tol
        keep(i) = false;
    end
end

P = P(keep,:);

if size(P,1) > 1 && norm(P(end,:) - P(1,:)) <= tol
    P(end,:) = [];
end

end


function x = uniqueWithinTolerance(x, tol)

if isempty(x)
    return
end

keep = true(size(x));

for i = 2:numel(x)
    if abs(x(i)-x(i-1)) <= tol
        keep(i) = false;
    end
end

x = x(keep);

end
