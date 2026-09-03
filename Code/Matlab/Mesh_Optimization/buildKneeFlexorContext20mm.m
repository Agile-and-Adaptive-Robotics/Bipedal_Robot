function ctx = buildKneeFlexorContext20mm()
    % all the transform setup code here



%% Joint rotation transformation matrices
positions = 100;
fprintf('The algorithm will be calculating Torque at %d different joint positions.\n', positions)

R = zeros(3, 3, positions);
T = zeros(4, 4, positions);
R_Pam = zeros(3, 3, positions);
T_Pam = zeros(4, 4, positions);
T_Pam_inv = zeros(size(T_Pam));
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
phi = linspace(kneeMin, kneeMax, positions);
%We want one of our positions to be home position, so let's make the
%smallest value of phi equal to 0
[val, pos] = min(abs(phi));
phi(pos) = 0;

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
    
    T_Pam(:, :, i) = RpToTrans(R_Pam(:, :, i), hipToKnee_Pam');     %Transformation matrix for robot, ICR to femur
    T_Pam_inv(:,:,i) = T_Pam(:,:,i)\eye(4);   % T_Pam_inv maps femur -> ICR;
    t1toICR(1,:,i) = [fcn13(phi(i)), fcn14(phi(i)), 0]; %distance from theta1 to ICR
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR(1,:,i)');    %transform from the ICR frame to theta1
    T_ICR_t1(:, :, i) = RpToTrans(eye(3), -t1toICR(1,:,i)');    %transform from t1 frame to ICR
end

% Precompute the femur-to-t1 transformation used by the route builder.
% T_f_t1 maps t1 -> femur, T_t1_f maps femur frame to t1
T_t1_f = pagemtimes(T_t1_ICR, T_Pam_inv);
T_f_t1 = zeros(size(T_t1_f));
for i = 1:positions
    T_f_t1(:,:,i) = T_t1_f(:,:,i)\eye(4);
end

ctx.T_Pam      = T_Pam;
ctx.T_Pam_inv  = T_Pam_inv;
ctx.T_ICR_t1   = T_ICR_t1;
ctx.T_t1_ICR   = T_t1_ICR;
ctx.T_t1_f     = T_t1_f;
ctx.T_f_t1     = T_f_t1;
ctx.phi        = phi(:);
ctx.phiD       = phi(:)*180/pi;
ctx.pos        = pos;
ctx.N          = numel(ctx.phiD);

%% Build optimization context

ctx.Name       = 'Bicep Femoris (Short Head)';
ctx.CrossPoint = 2;
ctx.Dia        = 20;

ctx.fitting    = 0.021;
ctx.geo = buildGeoExclusion();
ctx.wraps      = 6;
ctx.BPAcount   = 1;

% Intermediate t1-frame routing contact.  The old human-model wrap seed is
% retained in x-y.  buildKneeFlexorRoute20mm places its z coordinate on the
% full-extension p1-pEnd line projected into the t1 y-z plane.
ctx.wrapPointT1XY = [-0.0026, -0.00587];
ctx.wrapBendRadius = 0.025;       % local bend radius used by Xi3, m
ctx.wrapAngleToleranceD = 0;      % direct atan2 comparison; no added margin

% Choose extension angle for the resting-length constraint.
% Use 10 if you trust kneeMax; use 5 if the CAD hard stop is real.
ctx.extensionAngleD = 5;
[~, ctx.idxExtension] = min(abs(ctx.phiD - ctx.extensionAngleD));

% OpenSim target: column 2 = knee angle, column 4 = bifemsh_r.
H = readmatrix('OpenSim_Bifem_Results.txt', ...
    'FileType', 'text', ...
    'NumHeaderLines', 7);

ctx.humanAngleD = H(:,2);
ctx.humanTorque = H(:,4);
ctx.humanTorqueAbs = abs(ctx.humanTorque);

% Use 620 kPa as the "can this BPA meet the human target?" condition.
ctx.targetPressure = 620;
ctx.Pbins = 620;

ctx.KMAX = 0.255;          % KMAX = (rest - kmax)/rest at 620 kPa
ctx.maxRelStrain = 1.0;    % allow relative strain up to 1
ctx.minStrain = -0.03;
ctx.requiredTorqueMargin = 0.01;
ctx.offAxisPenaltyWeight = 1e3;
ctx.wrapReleasePenaltyWeight = 1e3;

ctx.torqueScale = max(1, max(ctx.humanTorqueAbs));

% Replace these with your best current estimates if desired
load minimizeFlxPin10_results_20260730_2transforms_Z2.mat filtered_results xCols
pick = 1;
g = filtered_results(pick,xCols);
Xi0 = g(1);
Xi1 = g(2);
Xi2 = g(3);

clear filtered_results xCols pick
load minimizeExtPin10_results_20260819_2transforms_Z2.mat filtered_results xCols
pick = 1;
g(4) = filtered_results(pick,xCols(4));
Xi3 = g(4);

% Fixed stiffness / compliance parameters.
% These are already identified and are NOT optimization variables.
ctx.Xi0 = Xi0;
ctx.Xi1 = Xi1;
ctx.Xi2 = Xi2;
ctx.Xi3 = Xi3;

ctx.wraps = 6;

% Fixed geometry from the previous two-point, no-wrap test.  These values
% define the original BPA curves and the separate Tx/Ty limits.  Do not
% change them when moving the initial point for a new wrapped search.
p10_previous = [-0.050,   0.035,   0.050];
p20_previous = [-0.01224, -0.00887, 0.02787];

% Initial geometry for this wrapped optimization.  Change p2_0 here to
% move the new test's search point without changing the previous test.
p1_0 = p10_previous;
p2_0 = [-0.01224, -0.00887, 0.02787];  % NEW WRAPPED-TEST SEARCH POINT

rest0     = 0.415;
tendon0   = 0.025;

%from human model
% p1_0 = [0.005, -0.211, 0.023];
% % p2_0 = [-0.0026   -0.0105    0.0290];         %wrapping point
% % p2_0 = [0.0044   -0.0305    0.0340];           %insertion point
% p2_0 = [-0.0026   -0.0305    0.0340];           %insertion point, modified X value
% 
% 
% rest0     = 0.225;
% tendon0   = 0.025;

% Optimization variables:
% [p1(1:3), pEnd(1:3), rest, tendon]
x0 = [p1_0, p2_0, rest0, tendon0];

lb = x0;
ub = x0;

% Attachment search box, meters
lb(1:3) = p1_0 + [-0.075, -0.100, -0.075];
ub(1:3) = p1_0 + [0.030,  0.100,  0];

lb(4:6) = p2_0 + [-0.100, -0.2, -0.05];
ub(4:6) = p2_0 + [ 0.025, -0.037,  0.015];

% BPA / tendon bounds
lb(7) = 0.415;    ub(7) = 1.200;    % rest length, m
lb(8) = 0.025;    ub(8) = 0.20;    % tendon length, m

ctx.x0 = x0;
ctx.lb = lb;
ctx.ub = ub;

% Original two-point BPA used for plots and the separate x/y off-axis
% limits.  It uses the X3 class with a zero bend measure, so no wrapping
% point or X3 length loss is introduced into the original calculation.
ctx.originalP1 = p10_previous;
ctx.originalPEnd = p20_previous;
ctx.originalRest = 0.415;
ctx.originalTendon = 0.015;
ctx.originalWraps = 6;

predOriginal = predictOriginalKneeFlexor20mm(ctx);
if ~predOriginal.ok || any(~isfinite(predOriginal.Torque(:)))
    error('buildKneeFlexorContext20mm:X0Prediction', ...
        'The original prediction failed or produced nonfinite torque: %s', ...
        predOriginal.failReason)
end
ctx.originalTorqueX = predOriginal.TorqueX(:);
ctx.originalTorqueY = predOriginal.TorqueY(:);
ctx.offAxisTorqueScaleX = max(1, max(abs(ctx.originalTorqueX)));
ctx.offAxisTorqueScaleY = max(1, max(abs(ctx.originalTorqueY)));

end
