% Robot PAM Calculation Optimization
% Author: Connor Morrow
% Date: 09/30/2020
% Description: This script calculate the torque generated across the back joint by
% the PAMs that cross over them. It is intended to be used by the 
% optimization script.

if exist('divisions', 'var') == 0
    divisions = 100;
end

if exist('ChooseJoint', 'var') == 0
    ChooseJoint = 'Bi_Hip';
end

%Flag that tells the code to not run pieces that have been calculated
%previously
if exist('beginOptimization', 'var') == 0
    beginOptimization = 0;
end

if beginOptimization == 0
    %Hip Joint z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 110*pi/180;
    MinTheta = -15*pi/180;
    Home = [-0.0707; -0.0661; 0.0835];                   %Home position from the Tibia
    Joint1 = JointData('HipZ', Raxis, MaxTheta, MinTheta, Home, divisions);

    T_h = Joint1.TransformationMat;

    %Knee z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 0.17453293;
    MinTheta = -2.0943951;
    Home = [-0.0707; -0.0661; 0.0835];      %Not actually the home position, but a placeholder
    Joint2 = JointData('Knee', Raxis, MaxTheta, MinTheta, Home, divisions);

    knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
    knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
    KneeXFunc = fit(knee_angle_x,knee_x,'cubicspline');
    knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
    knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
    KneeYFunc = fit(knee_angle_y,knee_y,'cubicspline');

    kneeHome = zeros(3, divisions);
    for j = 1:divisions
        kneeHome(:, j) = [KneeXFunc(Joint2.Theta(j)); KneeYFunc(Joint2.Theta(j)); 0];
    end

    T = zeros(4, 4, 2, divisions*divisions);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_h(:,:, i);
            T(:, :, 2, ii) = [Joint2.RotMat(:, :, j), kneeHome(:, j); 0, 0, 0, 1]; %Change the knee position
            ii = ii + 1;
        end
    end
end

%Bicep Femoris, Long Head p5 -> t5  
if beginOptimization == 0
    %Bifemlh
    Location1 = [-0.126, -0.03, -0.023, 0.005;
                    -0.103, -0.036, -0.056, -0.378;
                    0.069, 0.029, 0.034, 0.035];
    BFCrossPoints = [2 2];
    BFMIF = 896;  %max isometric force
    Axis1 = [10 20 30;
            0, 0, 30];                           %The axis of interest when calculating the moment arm about each joint. Looking at x, y, z for hip, and only z for knee

   %For Semimem, p5 -> semwr1 -> semwr2 -> twr1 -> t1 -> t6
   LocationSem = [-0.126, -0.041, -0.006, 0.022, 0.047, 0.038; 
                  -0.103, -0.237, -0.404, -0.052, -0.113, -0.382;
                  0.069, -0.033, -0.057, -0.03, -0.002, 0.004];
   SemCrossPoints = [2, 4];
   SemMIF = 1288+410;
   SemAxis = [10, 20, 30;
              0, 0, 30];

   %For Sar, p3 -> sarwr1 -> sarwr3 -> twr1 -> t1
%    LocationSar = [-0.026, 0.091, 0.098, 0.038, 0.022, 0.047;
%                   0.003, -0.105, -0.276, -0.017, -0.052, -0.113;
%                   0.138, 0.037, -0.4, -0.044, -0.03, -0.002];

    LocationSar = [-0.026, 0.091, 0.098, 0.037, 0.022, 0.047, 0.038;
                  0.003, -0.105, -0.276, -0.4, -0.052, -0.113, -0.382;
                  0.138, 0.037, -0.017, -0.044, -0.03, -0.002, -0.004];
    SarCrossPoints = [2, 5];
    SarMIF = 156;
    SarAxis = [10, 20, 30;
              0, 0, 30];

  %For tfl, p3 -> tflwr1 -> tf1wr2 -> t2 -> t5
   LocationTfl = [-0.026, 0.037, 0.005, 0.006, 0.005;
                  0.003, -0.037, -0.405, -0.062, -0.378;
                  0.138, 0.071, 0.047, 0.047, 0.035];
   TflCrossPoints = [2, 4];
   TflMIF = 233;
   TflAxis = [10, 20, 30;
              0, 0, 30];

  %------------- going to try to make Rectus Femoris into a function
  %------------- god have mercy on us all.
  %For Rect_fem, p10 -> rectwr1 -> rectwr3 -> pal(function [fcn3, fcn4,
  %0.002];
  rect_fem_x = [0.0156367; 0.0179948; 0.0274274; 0.029683; 0.0306; 0.0366; 0.0422; 0.0451; 0.0484; 0.0533; 0.0617; 0.0634; 0.067; 0.0733];
  rect_fem_xD = pi/180*[-120.118; -114.871; -90.068; -83.532; -80; -60; -40; -30; -20; -10; 0; 1.6; 5; 10];
  fcn3 = fit(rect_fem_xD,rect_fem_x,'smoothingspline');
  rect_fem_y =[0.0234; 0.0238; 0.0251; 0.0253; 0.025284; 0.0249; 0.0243; 0.0239; 0.0234; 0.0228; 0.0210; 0.0206; 0.0192; 0.0160];
  rect_fem_yD = pi/180*[-120; -114.6; -90; -83.5; -80.01; -60; -40; -30; -20; -10; 0; 1.6; 5; 10 ];
  fcn4 = fit(rect_fem_yD,rect_fem_y,'smoothingspline');

  LocationRec = zeros(3, 5, size(T, 4));
  iii = 1;
  for i = 1:size(Joint1.Theta, 2)
      for ii = 1:size(Joint2.Theta, 2)            %For Every joint1 angle, look at every possible knee angle
          LocationRec(:, :, iii) = [-0.072, 0.072, 0.065, 0.058, fcn3(Joint2.Theta(ii));
                                  0.125, -0.117, -0.156, -0.329, fcn4(Joint2.Theta(ii));
                                  0.079, 0.015, 0.005, -0.003, 0.002];
          iii = iii+1;
      end
  end
    RecCrossPoints = [2, 5];
    RecMIF = 1169;
    RecAxis = [10, 20, 30;
              0, 0, 30];

    %For Grac tqKZ, p8 -> gracwr1 -> gracwr2 -> twr1 -> t1 -> t6
    LocationGrac = [-0.074, -0.001, 0.018, 0.022, 0.047, 0.038;
                   0.125, -0.156, -0.401, -0.052, -0.113, -0.382;
                   0.079, -0.036, -0.053, -0.03, -0.002, 0.004];
    GracCrossPoints = [2, 4];
    GracMIF = 162;
    GracAxis = [10, 20, 30;
               0, 0, 30];


    %For bifemsh, f2 -> twr2 -> twr3
    LocationBif = [-0.007, -0.03, -0.023; 
                  -0.096, -0.036, -0.056;
                  0.032, 0.029, 0.029];
    BifCrossPoints = [0, 2];
    BifMIF = 804;
    BifAxis = [0 0 0;
              0, 0, 30];

    %Another function trial
    %For Vas_int, f1 -> vasintwr1 -> vasintwr2 -> pa2 (function [fcn7,
    %fcn8, 0.002])
    fcn7 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.008285),'smoothingspline');
    fcn8 = fit(rect_fem_yD,rect_fem_y+(0.0256239-0.02346),'smoothingspline');
    LocationVas = zeros(3, 4, size(T, 4));
    iii = 1;
    for i = 1:size(Joint1.Theta, 2)
      for ii = 1:size(Joint2.Theta, 2)            %For Every joint1 angle, look at every possible knee angle
          LocationVas(:, :, iii) = [0.0, 0.057, 0.054, fcn7(Joint2.Theta(ii));
                                    -0.043, -0.216, -0.332, fcn8(Joint2.Theta(ii));
                                    0.042, 0.029, 0.013, 0.002];
          iii = iii+1;
      end
    end
    VasCrossPoints = [0, 4];
    VasMIF = 1365+1294+1871;
    VasAxis = [0 0 0;
              0, 0, 30];
          
    Location = {Location1, LocationSem, LocationSar, LocationTfl, LocationRec, LocationGrac, LocationBif, LocationVas};
end

Muscle1 = PamData('Bicep Femoris, Long Head', Location{1}, BFCrossPoints, BFMIF, T, Axis1);
Muscle2 = PamData('Semimembranosus', Location{2}, SemCrossPoints, SemMIF, T, SemAxis);
Muscle3 = PamData('Sartorius', Location{3}, SarCrossPoints, SarMIF, T, SarAxis);
Muscle4 = PamData('Tensor Fasciae Latae', Location{4}, TflCrossPoints, TflMIF, T, TflAxis);
Muscle5 = PamData('Rectus Femoris', Location{5}, RecCrossPoints, RecMIF, T, RecAxis);
Muscle6 = PamData('Gracilis', Location{6}, GracCrossPoints, GracMIF, T, GracAxis);
Muscle7 = PamData('Bicep Femoris, Short Head', Location{7}, BifCrossPoints, BifMIF, T, BifAxis);
Muscle8 = PamData('Vastus Intermedius', Location{8}, VasCrossPoints, VasMIF, T, VasAxis);

%Sort Torques, to make equation management easier
T1 = Muscle1.Torque;
T2 = Muscle2.Torque;
T3 = Muscle3.Torque;
T4 = Muscle4.Torque;
T5 = Muscle5.Torque;
T6 = Muscle6.Torque;
T7 = Muscle7.Torque;
T8 = Muscle8.Torque;

Torque1x = T1(1, :, 1)+T2(1, :, 1)+T3(1, :, 1)+T4(1, :, 1)+T5(1, :, 1)+T6(1, :, 1);
Torque1y = T1(1, :, 2)+T2(1, :, 2)+T3(1, :, 2)+T4(1, :, 2)+T5(1, :, 2)+T6(1, :, 2);
Torque1z = T1(1, :, 3)+T2(1, :, 3)+T3(1, :, 3)+T4(1, :, 3)+T5(1, :, 3)+T6(1, :, 3);
Torque2a = T1(2, :, 1)+T2(2, :, 3)+T3(2, :, 3)+T4(2, :, 3)+T5(2, :, 3)+T6(2, :, 3);
Torque2b = T7(2, :, 3)+T8(2, :, 3);

%Create the Mesh of Torques to corespond with the joint angles
Torque1xM = zeros(divisions, divisions); Torque1yM = zeros(divisions, divisions); Torque1zM = zeros(divisions, divisions); Torque2aM = zeros(divisions, divisions); Torque2bM = zeros(divisions, divisions); 
for i = 1:divisions
    Torque1xM(:, i) = Torque1x(((i-1)*divisions)+1:i*divisions);
    Torque1yM(:, i) = Torque1y(((i-1)*divisions)+1:i*divisions);
    Torque1zM(:, i) = Torque1z(((i-1)*divisions)+1:i*divisions);
    Torque2aM(:, i) = Torque2a(((i-1)*divisions)+1:i*divisions);
    Torque2bM(:, i) = Torque2b(((i-1)*divisions)+1:i*divisions);
end

%Including Generic variable name for plotting
RobotAxis1 = Joint1.Theta;
RobotAxis1Label = 'Hip Flexion, Degrees';
RobotAxis2 = Joint2.Theta;
RobotAxis2Label = 'Knee Flexion, Degrees';
RobotTorque1 = Torque1xM;
RobotTorque2 = Torque1yM;
RobotTorque3 = Torque1zM;
RobotTorque4 = Torque2aM;
RobotTorque5 = Torque2bM;
RobotTitle1 = 'Robot Hip Torque, Biarticulate Muscles, X axis';
RobotTitle2 = 'Robot Hip Torque, Biarticulate Muscles, Y axis';
RobotTitle3 = 'Robot Hip Torque, Biarticulate Muscles, Z axis';
RobotTitle4 = 'Robot Knee Torque, Biarticulate Muscles, Z axis';
RobotTitle5 = 'Robot Knee Torque, Uniarticulate Muscles, Z axis';


%Store all necessary variables into a cell array for the optimization
%script
Muscles = {Muscle1, Muscle2, Muscle3, Muscle4, Muscle5, Muscle6, Muscle7, Muscle8};