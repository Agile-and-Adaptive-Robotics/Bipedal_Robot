function ctx = buildKneeExtContext20mm()
    % all the transform setup code here



%% Joint rotation transformation matrices
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);
R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
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
    
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0]; %distance from theta1 to ICR
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR(1,:,i)');    %transform from the ICR frame to theta1
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR(1,:,i)');    %transform from t1 frame to ICR
end

%% Build Optimization Context
ctx.Name       = 'Vastus Medialis Proximal Ring 20mm BPA';
ctx.CrossPoint = 5;
ctx.Dia        = 20;
ctx.BPAcount   = 2;

ctx.T          = T;
ctx.T_Pam      = T_Pam;
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


%% Distal-ring route initial geometry
% These are the row values from the Distal Ring Location matrix in
% Knee_Extensor_40mm. Rows 1 and 6 are the optimizer attachment variables.
%
% p1_0: row 1, femur-side BPA attachment
% p2_0: row 6, tibia/theta1-side distal ring attachment
[~, idxLongest] = min(phiD);

raw = distalRingRawLocation(phiD(idxLongest));
p1_0 = raw(1,:);
p2_0 = raw(8,:);

% Tendon initial guess:
% Use the norm between rows 5 and 6 of the distal ring Location matrix.
tendon0 = norm(raw(7,:) - raw(8,:));
fprintf('Initial tendon from norm(row7-row8) = %.6f m\n', tendon0)

% Resting length initial guess:
% Use undeformed normal musculotendon length:
%   rest0 = max(Lmt_undeformed) - 2*fitting - tendon0 - Xi0
Location0 = buildDistalRingLocation(p1_0, p2_0, ctx.phiD, ctx.T_ICR_t1);
Lmt0 = muscleLengthNormal(Location0, ctx.CrossPoint, ctx.T_Pam);

rest0 = max(Lmt0) - 2*ctx.fitting - tendon0 - ctx.Xi0;

if rest0 <= 0
    error('Initial rest0 is non-positive. Check distal ring path geometry.')
end

%% Optimization variables:
% [p1(1:3), p2(1:3), rest, kmaxFrac, tendon]
x0 = [p1_0, p2_0, rest0, tendon0];

lb = x0;
ub = x0;

% Attachment search box, meters
lb(1:3) = p1_0 + [-0.060, -0.060, -0.060];
ub(1:3) = p1_0 + [ 0.100,  0.100,  0.060];

lb(4:6) = p2_0 + [0, -0.060, -0.060];
ub(4:6) = p2_0 + [ 0.06,  0.06,  0.060];

% BPA / tendon bounds
lb(7) = 0.360;    ub(7) = 0.520;    % rest length, m
lb(8) = max(0.025,tendon0-0.040);    ub(8) = tendon0+0.1;    % tendon length, m

ctx.x0 = x0;
ctx.lb = lb;
ctx.ub = ub;

%% Initial-design diagnostics
ctx.initialRawLocation = raw;
ctx.initialLocation = Location0;
ctx.initialLmt0 = Lmt0;
ctx.initialMaxLmt0 = max(Lmt0);
ctx.initialIdxMaxLmt0 = find(Lmt0 == max(Lmt0), 1, 'first');
ctx.initialMaxLmtAngleD = ctx.phiD(ctx.initialIdxMaxLmt0);
ctx.initialRest0 = rest0;
ctx.initialTendon0 = tendon0;

end


%% Local helper functions

function raw = distalRingRawLocation(phiD_i)
% Location matrix from Knee_Extensor_10mm.
%
% Rows 1:4 are femur / parent frame.
% Rows 5:8 are theta1 / tibia frame before conversion to ICR.
%
% Original Knee_Extensor_10mm points:
% p1 = [0.040, 0.035, 0];                 %Origin
% p2 = [0.0759, -0.27476, 0];             %BPA contacts mounting base
% p3 = [0.05982, -0.37427, 0.000];        %femur channel contact, updated
% p4 = [0.03955, -0.42183, 0.000];        %femoral condyle contact, updated
% p5 = [0.05871, 0.025, 0];               %Tibia contact initial
% p6 = [0.05871, 0.01228, 0];             %Tibia contact, updated
% p7 = [0.054, -0.00612, 0];              %Tibia tendon contact, updated
% p8 = [0.03604, -0.02844, 0];            %patellar ligament ring

p1 = [0.040,   0.035,    0.000];          %Origin
p2 = [0.0759, -0.27476,  0.000];          %BPA contacts mounting base
p3 = [0.05982,-0.37427,  0.000];          %femur channel contact, updated
p4 = [0.03955,-0.42183,  0.000];          %femoral condyle contact, updated
p5 = [0.05871, 0.025,    0.000];          %Tibia contact initial
p6 = [0.05871, 0.01228,  0.000];          %Tibia contact, updated
p7 = [0.054,  -0.00612,  0.000];          %Tibia tendon contact, updated
p8 = [0.03604,-0.02844,  0.000];          %patellar ligament ring

%Set up angle limits (degrees)
% AA = -68;
BB = -69.8;   %add p4 when flexion reaches this value
CC = -9.19;   %add p3 and p5 when flexion reaches this value
DD = 10;

    if phiD_i > CC
        raw = [p1;
               p2;
               p2;
               p2;
               p6;
               p6;
               p7;
               p8];
    
    elseif phiD_i <= CC && phiD_i > BB
        raw = [p1;
               p2;
               p3;
               p3;
               p5;
               p6;
               p7;
               p8];
    
    elseif phiD_i <= BB
        raw = [p1;
               p2;
               p3;
               p4;
               p5;
               p6;
               p7;
               p8];
    
    else
        raw = zeros(8,3);
    end

end

function Location = buildDistalRingLocation(p1, p2, phiD, T_ICR_t1)
% Build the extensor Location matrix using the Knee_Extensor_10mm route.
%
% p1 replaces row 1.
% p2 replaces row 8.
% Rows 5:8 are converted from theta1 -> ICR with T_ICR_t1(:,:,i).

N = numel(phiD);
Location = zeros(8, 3, N);

    for i = 1:N
        raw = distalRingRawLocation(phiD(i));
    
        raw(1,:) = p1;
        raw(8,:) = p2;
    
        for j = 5:8
            raw(j,:) = RowVecTrans(T_ICR_t1(:,:,i), raw(j,:));
        end
    
        Location(:,:,i) = raw;
    end

end

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