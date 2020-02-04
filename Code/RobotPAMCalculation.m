% Robot PAM Calculation
% Author: Connor Morrow
% Date: 12/11/2019
% Description: This script calculate the torque generated across joints by
% the PAMs that cross over them. 

if exist('divisions', 'var') == 0
    divisions = 100;
end

if exist('ChooseJoint', 'var') == 0
    ChooseJoint = 'Back';
end

if isequal(ChooseJoint, 'Back')
    %Back x axis rotation
    Raxis = [1; 0; 0];
    MaxTheta = 15*pi/180;
    MinTheta = -15*pi/180;
    Home = [-0.1007; 0.0815; 0];                   
    BackX = JointData('BackX', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    %Back y rotation
    Raxis = [0; 1; 0];
    MaxTheta = 0;
    MinTheta =0;               %Ben's script lists the min and max, however only 0 is used. Why?
    Home = [-0.1007; 0.0815; 0];
    BackY = JointData('BackY', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    %Back z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 15*pi/180;
    MinTheta = -45*pi/180;
    Home = [-0.1007; 0.0815; 0];
    BackZ = JointData('BackZ', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    R = zeros(3, 3, 1, divisions^2);
    T = zeros(4, 4, 1, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            R(:, :, 1, ii) = BackX.RotMat(:, :, i)*BackY.RotMat(:, :, i)*BackZ.RotMat(:, :, j);
            T(:, :, 1, ii) = [R(:, :, 1, ii), Home; 0 0 0 1];
            ii = ii+1;
        end
    end
    
    T(:, :, 2, :) = T(:, :, 1, :);
    T(:, :, 3, :) = T(:, :, 1, :);
    
    %Erector Spinae p7 -> b2    
    ESLocation = [-0.122, -0.158, -0.072, -0.039, 0.005, -0.032, -0.072, -0.04, -0.019;
                -0.051, 0.054, 0.172, 0.179, 0.027, 0.027, 0.175, 0.182, 0.048;
                0.078, 0.041, 0.031, 0.029, 0.047, 0.048, 0.048, 0.045, 0.059];
    ESCrossPoints = 3;                    %Via points are the points where a transformation matrix is needed. Typically wrap point + 1
    ESMIF = 2500;
    Axis1 = [10, 10];                           %The axis of interest when calculating the moment arm about each joint. The axis is 1, but is listed as 10 so that the cross product doesn't rotate the resulting vector. See PamData > CrossProd
    Muscle1 = PamData('Erector Spinae', ESLocation, ESCrossPoints, ESMIF, T, Axis1);
    
        Axis1 = [20, 10];                           %The axis of interest when calculating the moment arm about each joint. The axis is 1, but is listed as 10 so that the cross product doesn't rotate the resulting vector. See PamData > CrossProd
    Muscle2 = PamData('Erector Spinae', ESLocation, ESCrossPoints, ESMIF, T, Axis1);
    
        Axis1 = [30, 10];                           %The axis of interest when calculating the moment arm about each joint. The axis is 1, but is listed as 10 so that the cross product doesn't rotate the resulting vector. See PamData > CrossProd
    Muscle3 = PamData('Erector Spinae', ESLocation, ESCrossPoints, ESMIF, T, Axis1);
    

    TorqueX = Muscle1.Torque(1, :, :);
    TorqueY = Muscle2.Torque(1, :, :);
    TorqueZ = Muscle3.Torque(1, :, :);
    
    %Create the Mesh of Torques to corespond with the joint angles
    TorqueMX = zeros(divisions, divisions); TorqueMY = zeros(divisions, divisions); TorqueMZ = zeros(divisions, divisions);
    for i = 1:divisions
        TorqueMX(:, i) = TorqueX(((i-1)*divisions)+1:i*divisions);
        TorqueMY(:, i) = TorqueY(((i-1)*divisions)+1:i*divisions);
        TorqueMZ(:, i) = TorqueZ(((i-1)*divisions)+1:i*divisions);
    end
    
    %Transposing to match Ben's results
    TorqueMX = TorqueMX';
    TorqueMY = TorqueMY';
    TorqueMZ = TorqueMZ';
    
    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    RobotAxis1 = BackZ.Theta;
    RobotAxis1Label = 'Back Flexion, Degrees';
    RobotAxis2 = BackX.Theta;
    RobotAxis2Label = 'Lumbar Bending, Degrees';
    RobotTorque1 = TorqueMX;
    RobotTorque2 = TorqueMY;
    RobotTorque3 = TorqueMZ;
    RobotTitle1 = 'Robot Back Torque, Erector Spinae, X Axis';
    RobotTitle2 = 'Robot Back Torque, Erector Spinae, Y Axis';
    RobotTitle3 = 'Robot Back Torque, Erector Spinae, Z Axis';
    
elseif isequal(ChooseJoint, 'Bi_Hip')
    %Hip Joint z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 110*pi/180;
    MinTheta = -15*pi/180;
    Home = [-0.0707; -0.0661; 0.0835];                   %Home position from the Tibia
    HipZ = JointData('HipZ', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    T_h = HipZ.TransformationMat;
    
    %Knee z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 0.17453293;
    MinTheta = -2.0943951;
    Home = [-0.0707; -0.0661; 0.0835];      %Not actually the home position, but a placeholder
    Knee = JointData('Knee', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
    knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
    KneeXFunc = fit(knee_angle_x,knee_x,'cubicspline');
    knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
    knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
    KneeYFunc = fit(knee_angle_y,knee_y,'cubicspline');
    
    kneeHome = zeros(3, divisions);
    for j = 1:divisions
        kneeHome(:, j) = [KneeXFunc(Knee.Theta(j)); KneeYFunc(Knee.Theta(j)); 0];
    end
    
    T = zeros(4, 4, 2, divisions*divisions);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_h(:,:, i);
            T(:, :, 2, ii) = [Knee.RotMat(:, :, j), kneeHome(:, j); 0, 0, 0, 1]; %Change the knee position
            ii = ii + 1;
        end
    end
    
    %Bicep Femoris, Long Head p5 -> t5
    p5 = [-0.126; -0.103; 0.069]; %origin: semimem/semiten, bifemlh
    twr2 = [-0.03; -0.036; 0.029]; %common wrapping point for bifemlh, bifemsh
    twr3 = [-0.023; -0.056; 0.034]; %insertion: bifemsh. wrapping point: bifemlh
    t5 = [0.005; -0.378; 0.035]; %insertion: tfl, bifemlh
    
    BFLocation = [-0.126, -0.03, -0.023, 0.005;
                    -0.103, -0.036, -0.056, -0.378;
                    0.069, 0.029, 0.034, 0.035];
    BFCrossPoints = [2 2];
    BFMIF = 896;  %max isometric force
    
    Axis1 = [10; 30];                           %The axis of interest when calculating the moment arm about each joint
    Muscle1 = PamData('Bicep Femoris X', BFLocation, BFCrossPoints, BFMIF, T, Axis1);
    
    Axis2 = [20; 30];                           %The axis of interest when calculating the moment arm about each joint
    Muscle2 = PamData('Bicep Femoris Y', BFLocation, BFCrossPoints, BFMIF, T, Axis2);
    
    Axis3 = [30; 30];                           %The axis of interest when calculating the moment arm about each joint
    Muscle3 = PamData('Bicep Femoris Z', BFLocation, BFCrossPoints, BFMIF, T, Axis3);
    
    HipXTorqueR = Muscle1.Torque(1, :, :);
    HipYTorqueR = Muscle2.Torque(1, :, :);
    HipZTorqueR = Muscle3.Torque(1, :, :);
    KneeTorqueR = Muscle1.Torque(2, :, :);
    
    %Create the Mesh of Torques to corespond with the joint angles
    HipXTorqueMR = zeros(divisions, divisions); HipYTorqueMR = zeros(divisions, divisions); HipZTorqueMR = zeros(divisions, divisions); KneeTorqueMR = zeros(divisions, divisions); 
    for i = 1:divisions
        HipXTorqueMR(:, i) = HipXTorqueR(((i-1)*divisions)+1:i*divisions);
        HipYTorqueMR(:, i) = HipYTorqueR(((i-1)*divisions)+1:i*divisions);
        HipZTorqueMR(:, i) = HipZTorqueR(((i-1)*divisions)+1:i*divisions);
        KneeTorqueMR(:, i) = KneeTorqueR(((i-1)*divisions)+1:i*divisions);
    end
    
    %Including Generic variable name for plotting
    RobotAxis1 = HipZ.Theta;
    RobotAxis1Label = 'Hip Flexion, Degrees';
    RobotAxis2 = Knee.Theta;
    RobotAxis2Label = 'Knee Flexion, Degrees';
    RobotTorque1 = HipXTorqueMR;
    RobotTorque2 = HipYTorqueMR;
    RobotTorque3 = HipZTorqueMR;
    RobotTorque4 = KneeTorqueMR;
    RobotTitle1 = 'Robot Hip Torque, bifemlh, X axis';
    RobotTitle2 = 'Robot Hip Torque, bifemlh, Y axis';
    RobotTitle3 = 'Robot Hip Torque, bifemlh, Z axis';
    RobotTitle4 = 'Robot Knee Torque, bifemlh, Z axis';
    
elseif isequal(ChooseJoint, 'Calves')
    %Create a Joint object that calculates things like transformation
    %matrix
    %Note: Currently only goes for the matrices at minimum theta. Will
    %implement some iteration process to go from minimum to maximum. 
        
    Raxis = [-0.10501355; -0.17402245; 0.97912632];
    MaxTheta = 20*pi/180;
    MinTheta = -50*pi/180;
    Home = [0; -0.43; 0];                   %Home position from the Tibia
    Ankle = JointData('Ankle', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    Raxis = [0.7871796; 0.60474746; -0.12094949];
    MaxTheta = 0;
    MinTheta = 0;               %Ben's script lists the min and max, however only 0 is used. Why?
    Home = [-0.04877; -0.04195; 0.00792];   %Home position from the Talus
    Subtalar = JointData('Subtalar', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    %Knee z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 0.17453293;
    MinTheta = -2.0943951;
    Home = [-0.0707; -0.0661; 0.0835];      %Not actually the home position, but a placeholder
    Knee = JointData('Knee', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
    knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
    KneeXFunc = fit(knee_angle_x,knee_x,'cubicspline');
    knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
    knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
    KneeYFunc = fit(knee_angle_y,knee_y,'cubicspline');
    
    kneeHome = zeros(3, divisions);
    for j = 1:divisions
        kneeHome(:, j) = [KneeXFunc(Knee.Theta(j)); KneeYFunc(Knee.Theta(j)); 0];
    end
    
    T_a = Ankle.TransformationMat(:, :, :);
    T_s = Subtalar.TransformationMat(:, :, :);
    
    
    %Create a Mesh of the Ankle and MTP. First two dimensions are the
    %values of the transformation matrics. The third dimension determines
    %the joint being observed, and the fourth dimension is its location
    %in the mesh, which will need to be converted from a one dimensional
    %array to two dimensions
    T = zeros(4, 4, 3, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_a(:, :, i);            %Keep the ankle at one point as the MTP changes
            T(:, :, 2, ii) = T_s(:, :, i);            %Subtalar doesn't change
            ii = ii + 1;
        end
    end
    
    %Soleus a1 -> t1
    SLocation = [0.047, 0;
                -0.113, 0.031;
                -0.002, -0.005];
    SCrossPoints = [2 2];
    SMIF = 3549;  %max isometric force
    Axis3 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle3 = PamData('Soleus', SLocation, SCrossPoints, SMIF, T, Axis3);
    
    T(:, :, 2:3, :) = T(:, :, 1:2, :);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = [Knee.RotMat(:, :, j), kneeHome(:, j); 0, 0, 0, 1]; %Change the knee position
            ii = ii + 1;
        end
    end
    
    %Medial Gastroocnemius a1 -> f4            
    MGLocation = [  0.015, -0.03, 0;
                    -0.374, -0.402, 0.031;
                    -0.025, -0.026, -0.005];
        
    MGCrossPoints = [3 3 3];
    MGMIF = 1558;  %max isometric force
    Axis1 = [3; 3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle1 = PamData('Medial Gastroocnemius', MGLocation, MGCrossPoints, MGMIF, T, Axis1);
    
    %Lateral Gastrocenemius a1 -> f5
    LGLocation = [0.009, -0.022, 0;
                    -0.378, -0.395, 0.031;
                    0.027, 0.027, -0.005];
                
    LGCrossPoints = [3 3 3];
    LGMIF = 683;  %max isometric force
    Axis2 = [3; 3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle2 = PamData('Lateral Gastrocenemius', LGLocation, LGCrossPoints, LGMIF, T, Axis2);
    
    %Torque Calcs, "R" for robot
    AnkleTorqueR = Muscle1.Torque(2, :, :) + Muscle2.Torque(2, :, :);
    AnkleSoleusTorqueR = Muscle3.Torque(1, :, :);
    KneeTorqueR = Muscle1.Torque(1, :, :) + Muscle2.Torque(1, :, :);
    
    %Create the Mesh of Torques to corespond with the joint angles
    AnkleTorqueMR = zeros(divisions, divisions); AnkleSoleusTorqueMR = zeros(divisions, divisions); KneeTorqueMR = zeros(divisions, divisions);
    for i = 1:divisions
        AnkleTorqueMR(:, i) = AnkleTorqueR(((i-1)*divisions)+1:i*divisions);
        AnkleSoleusTorqueMR(:, i) = AnkleSoleusTorqueR(((i-1)*divisions)+1:i*divisions);
        KneeTorqueMR(:, i) = KneeTorqueR(((i-1)*divisions)+1:i*divisions);
    end
    
    %Including Generic variable name for plotting
    RobotAxis1 = Ankle.Theta;
    RobotAxis1Label = 'Ankle Flexion, Degrees';
    RobotAxis2 = Knee.Theta;
    RobotAxis2Label = 'Knee Flexion, Degrees';
    RobotTorque1 = KneeTorqueMR;
    RobotTorque2 = AnkleTorqueMR;
    RobotTorque3 = AnkleSoleusTorqueMR;
    RobotTitle1 = 'Robot Knee Torque, Gastrocnemius, Z axis';
    RobotTitle2 = 'Robot Ankle Torque, Gastrocnemius, Z` axis';
    RobotTitle3 = 'Robot Ankle Torque, Soleus, Z` axis';
    
elseif isequal(ChooseJoint, 'Foot')
    Raxis = [-0.10501355; -0.17402245; 0.97912632];
    MaxTheta = 20*pi/180;
    MinTheta = -50*pi/180;
    Home = [0; -0.43; 0];                   %Home position from the Tibia
    Ankle = JointData('Ankle', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    Raxis = [0.7871796; 0.60474746; -0.12094949];
    MaxTheta = 35*pi/180;
    MinTheta = -25*pi/180;
    Home = [-0.04877; -0.04195; 0.00792];   %Home position from the Talus
    Subtalar = JointData('Subtalar', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    T_a = Ankle.TransformationMat(:, :, :);
    T_s = Subtalar.TransformationMat(:, :, :);

    %Create a Mesh of the Ankle and MTP. First two dimensions are the
    %values of the transformation matrics. The third dimension determines
    %the joint being observed, and the fourth dimension is its location
    %in the mesh, which will need to be converted from a one dimensional
    %array to two dimensions
    T = zeros(4, 4, 3, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_a(:, :, i);            %Keep the ankle at one point as the MTP changes
            T(:, :, 2, ii) = T_s(:, :, j);            
            ii = ii + 1;
        end
    end
    
    %Tibialis Posterior t4 -> a4
    TPLocation = [  -0.02, -0.014, 0.042, 0.077;
                    -0.085, -0.405, 0.033, 0.016;
                    0.002, -0.023, -0.029, -0.028];
        
    TPCrossPoints = [3 3];
    TPMIF = 3549;  %max isometric force
    Axis1 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle1 = PamData('Tibialis Posterior', TPLocation, TPCrossPoints, TPMIF, T, Axis1);
    
    %Tibialis Anterior t2 -> a3
    TALocation = [0.006, 0.033, 0.117;
                -0.062, -0.394, 0.014;
                0.047, 0.007, -0.036];
    TACrossPoints = [3 3];
    TAMIF = 905;  %max isometric force
    Axis2 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle2 = PamData('Tibialis Posterior', TALocation, TACrossPoints, TAMIF, T, Axis2);
    
    %Peroneus Brevis t2 -> a2
    PBLocation = [0.006, -0.02, -0.014, 0.047, 0.079;
                -0.062, -0.418, -0.429, 0.027, 0.022;
                0.047, 0.028, 0.029, 0.023, 0.034];
    PBCrossPoints = [4 4];
    PBMIF = 435;  %max isometric force
    Axis3 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle3 = PamData('Peroneus Brevis', PBLocation, PBCrossPoints, PBMIF, T, Axis3);
    
    %Peroneus Longus t2 -> a3
    PLLocation = [0.005, -0.021, -0.016, 0.044, 0.048, 0.085, 0.117;
                -0.062, -0.42, -0.432, 0.023, 0.011, 0.007, 0.014;
                0.047, 0.029, 0.029, 0.022, 0.028, 0.012, -0.036];
    PLCrossPoints = [4 4];
    PLMIF = 943;  %max isometric force
    Axis4 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle4 = PamData('Peroneus Longus', PLLocation, PLCrossPoints, PLMIF, T, Axis4);
    
    %Peroneus Tertius t2 -> a2
    PTLocation = [0.006, 0.023, 0.079;
                  -0.062, -0.407, 0.016;
                  0.079, 0.022, 0.034];
            
    PTCrossPoints = [3 3];
    PTMIF = 180;  %max isometric force
    Axis5 = [3; 1];                           %The axis of interest when calculating the moment arm about each joint
    Muscle5 = PamData('Peroneus Tertius', PTLocation, PTCrossPoints, PTMIF, T, Axis5);
    
    %Torque Calcs, "R" for robot
    AnkleTorqueR = Muscle1.Torque(1, :, :) + Muscle2.Torque(1, :, :) + Muscle3.Torque(1, :, :) + Muscle4.Torque(1, :, :) + Muscle5.Torque(1, :, :);
    SubtalarTorqueR = Muscle1.Torque(2, :, :) + Muscle2.Torque(2, :, :) + Muscle3.Torque(2, :, :) + Muscle4.Torque(2, :, :) + Muscle5.Torque(2, :, :);

    AnkleTorqueMR = zeros(divisions, divisions); SubtalarTorqueMR = zeros(divisions, divisions);
    %Create the Mesh of Torques to corespond with the joint angles
    for i = 1:divisions
        AnkleTorqueMR(:, i) = AnkleTorqueR(((i-1)*divisions)+1:i*divisions);
        SubtalarTorqueMR(:, i) = SubtalarTorqueR(((i-1)*divisions)+1:i*divisions);
    end
    
    %Including Generic variable name for plotting
    RobotAxis1 = Ankle.Theta;
    RobotAxis1Label = 'Ankle Flexion, Degrees';
    RobotAxis2 = Subtalar.Theta;
    RobotAxis2Label = 'MTP Flexion, Degrees';
    RobotTorque1 = AnkleTorqueMR;
    RobotTorque2 = SubtalarTorqueMR;
    RobotTitle1 = 'Robot Ankle Torque, Z Axis';
    RobotTitle2 = 'Robot Subtalar Torque, X Axis';
    
elseif isequal(ChooseJoint, 'Toe')
     
    Raxis = [-0.10501355; -0.17402245; 0.97912632];
    MaxTheta = 20*pi/180;
    MinTheta = -50*pi/180;
    Home = [0; -0.43; 0];                   %Home position from the Tibia
    Ankle = JointData('Ankle', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    Raxis = [0.7871796; 0.60474746; -0.12094949];
    MaxTheta = 0;
    MinTheta = 0;               %Ben's script lists the min and max, however only 0 is used. Why?
    Home = [-0.04877; -0.04195; 0.00792];   %Home position from the Talus
    Subtalar = JointData('Subtalar', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    Raxis = [-0.5809544; 0; 0.81393611];
    MaxTheta = 80*pi/180;
    MinTheta = -30*pi/180;
    Home = [0.1788; -0.002; 0.00108];       %Home Position from the Calcn
    MTP = JointData('MTP', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    T_a = Ankle.TransformationMat(:, :, :);
    T_s = Subtalar.TransformationMat(:, :, :);
    T_m = MTP.TransformationMat(:, :, :);
    
    %Create a Mesh of the Ankle and MTP. First two dimensions are the
    %values of the transformation matrics. The third dimension determines
    %the joint being observed, and the fourth dimension is its location
    %in the mesh, which will need to be converted from a one dimensional
    %array to two dimensions
    T = zeros(4, 4, 3, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_a(:, :, i);            %Keep the ankle at one point as the MTP changes
            T(:, :, 2, ii) = T_s(:, :, i);            %Subtalar doesn't change
            T(:, :, 3, ii) = T_m(:, :, j);            %Change the MTP while the subtalar stays the same
            ii = ii + 1;
        end
    end
    
    % Creating PAM class module
    %Looking at the flexor digitorus longus, but now for the bipedal robot
    %The robot uses the same physical skeletal geometry, but changes where
    %the PAMs are attached when compared to muscles. Because of this,
    %transformation matrices from the skeletal structure will be reused
    %from above. 

    %For the robot, the wrapping locations also include a global origin
    %point (in this case, t4) and a global insertion point (in this case,
    %a7).

    Axis1 = [3; 1; 3];                           %The axis of interest when calculating the moment arm about each joint

    %FDL Insertion points, found from AttachPoints_Robot. The first column is
    %the insertion point that shares a point with the tibia (t4). The last
    %point is a7, the insertion point for the foot
    FDLLocation = [-0.02, -0.015, 0.044, 0.071, 0.166, -0.002, 0.028, 0.044;
                    -0.085, -0.405, 0.032, 0.018, -0.008, -0.008, -0.007, -0.006;
                    0.002, -0.02, -0.028, -0.026, 0.012, 0.015, 0.022, 0.024];    
    FDLCrossPoints = [3, 3, 6];
    FDLMIF = 310;  %max isometric force
    FDL = PamData('Flexor Digitorum Longus', FDLLocation, FDLCrossPoints, FDLMIF, T, Axis1);
    Axis2 = [3; 1; 3];                           %The axis of interest when calculating the moment arm about each joint

    %FHL t4 to a8
    FHLLocation = [-0.02, -0.019, 0.037, 0.104, 0.173, 0.016, 0.056;
                    -0.085, -.408, 0.028, 0.007, -0.005, -0.006, -0.01;
                    0.002, -0.017, -0.024, -0.026, -0.027, -0.026, -0.018];
    FHLCrossPoints = [3, 3, 6];
    FHLMIF = 322;  %max isometric force
    FHL = PamData('Flexor Hallucis Longus', FHLLocation, FHLCrossPoints, FHLMIF, T, Axis1);

    %EDL, t2 to a5
    EDLLocation = [0.006, 0.029, 0.092, 0.162, 0.005, 0.044;
        -0.062, -0.401, 0.039, 0.006, 0.005, 0.0;
        0.047, 0.007, 0.0, 0.013, 0.015, 0.025];
    EDLCrossPoints = [3, 3, 5];
    EDLMIF = 512;  %max isometric force
    EDL = PamData('Extensor Digitorum Longus', EDLLocation, EDLCrossPoints, EDLMIF, T, Axis1);


    %EHL t2 to a6
    EHLLocation = [0.006, 0.033, 0.097, 0.129, 0.173, 0.03, 0.056;
                    -0.062, -0.398, 0.039, 0.031, 0.014, 0.004, 0.003;
                    0.047, -0.008, -0.021, -0.026, -0.028, -0.024, -0.019];
    EHLViaPoints = [3, 3, 6];
    EHLMIF = 162;  %max isometric force
    EHL = PamData('Extensor Hallucis Longus', EHLLocation, EHLViaPoints, EHLMIF, T, Axis1);

    %Torque Calcs, "R" for robot
    AnkleTorqueR = FDL.Torque(1, :, :) + FHL.Torque(1, :, :) + EDL.Torque(1, :, :) + EHL.Torque(1, :, :);
    SubtalarTorqueR = FDL.Torque(2, :, :) + FHL.Torque(2, :, :) + EDL.Torque(2, :, :) + EHL.Torque(2, :, :);
    MTPTorqueR = FDL.Torque(3, :, :) + FHL.Torque(3, :, :) + EDL.Torque(3, :, :) + EHL.Torque(3, :, :);

    %Create the Mesh of Torques to corespond with the joint angles
    AnkleTorqueMR = zeros(divisions, divisions); SubtalarTorqueMR = zeros(divisions, divisions); MTPTorqueMR = zeros(divisions, divisions);
    for i = 1:divisions
        AnkleTorqueMR(:, i) = AnkleTorqueR(((i-1)*divisions)+1:i*divisions);
        SubtalarTorqueMR(:, i) = SubtalarTorqueR(((i-1)*divisions)+1:i*divisions);
        MTPTorqueMR(:, i) = MTPTorqueR(((i-1)*divisions)+1:i*divisions);
    end
    
    %Including Generic variable name for plotting
    RobotAxis1 = Ankle.Theta;
    RobotAxis1Label = 'Ankle Flexion, Degrees';
    RobotAxis2 = MTP.Theta;
    RobotAxis2Label = 'MTP Flexion, Degrees';
    RobotTorque1 = AnkleTorqueMR;
    RobotTorque2 = SubtalarTorqueMR;
    RobotTorque3 = MTPTorqueMR;
    RobotTitle1 = 'Robot Ankle Torque, Z Axis';
    RobotTitle2 = 'Robot Subtalar Torque, X Axis';
    RobotTitle3 = 'Robot MTP Torque, Z Axis';
    
elseif isequal(ChooseJoint, 'Uni_Hip')
    %Hip Joint x rotation
    Raxis = [1; 0; 0];
    MaxTheta = 20.6*pi/180;
    MinTheta = -30*pi/180;
    HomeH = [-0.0707; -0.0661; 0.0835];                   %Home position from the Tibia
    HipX = JointData('HipX', Raxis, MaxTheta, MinTheta, HomeH, divisions);  
    
    %Hip Joint z rotation
    Raxis = [0; 0; 1];
    MaxTheta = 110*pi/180;
    MinTheta = -15*pi/180;
    HomeH = [-0.0707; -0.0661; 0.0835];                   %Home position from the Tibia
    HipZ = JointData('HipZ', Raxis, MaxTheta, MinTheta, HomeH, divisions);
    
    %Back rotation
    Raxis = [0; 0; 1];
    MaxTheta = 0;
    MinTheta = 0;
    HomeB = [-0.1007; 0.0815; 0];                   %Home position from the Tibia
    Back = JointData('Back', Raxis, MaxTheta, MinTheta, HomeB, divisions);
    
    R = zeros(3, 3, 1, divisions^2);
    T = zeros(4, 4, 1, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            R(:, :, 1, ii) = HipX.RotMat(:, :, j)*HipZ.RotMat(:, :, i);
            T(:, :, 1, ii) = [R(:, :, 1, ii), HomeH; 0 0 0 1];
            ii = ii+1;
        end
    end
    
    %Gluteus Maximus p6 -> f3
    GMLocation = [-0.162, -0.174, -0.069, -0.046, 0.005;
                  0.04, -0.033, -0.042, -0.116, -0.339;
                  0.023, 0.101, 0.067, 0.06, 0.022];
    GMCrossPoints = 3;
    GMMIF = 5047;  
    Axis1 = 10;
    Muscle1a = PamData('Gluteus Maximus', GMLocation, GMCrossPoints, GMMIF, T, Axis1);
    
    Axis1 = 20;
    Muscle1b = PamData('Gluteus Maximus', GMLocation, GMCrossPoints, GMMIF, T, Axis1);
    
    Axis1 = 30;
    Muscle1c = PamData('Gluteus Maximus', GMLocation, GMCrossPoints, GMMIF, T, Axis1);
    
    %Adductor Magnus p9 -> f4
    AMLocation = [-0.163, -0.059, -0.059, 0.015;
                  -0.013, -0.108, -0.108, -0.374;
                  0.005, -0.03, -0.03, -0.025];
    AMCrossPoints = 3;
    AMMIF = 2268;  
    Axis1 = 10;
    Muscle2a = PamData('Adductor Magnus', AMLocation, AMCrossPoints, AMMIF, T, Axis1);
    
    Axis1 = 20;
    Muscle2b = PamData('Adductor Magnus', AMLocation, AMCrossPoints, AMMIF, T, Axis1);
    
    Axis1 = 30;
    Muscle2c = PamData('Adductor Magnus', AMLocation, AMCrossPoints, AMMIF, T, Axis1);
    
    %Iliacus p1 -> f4
    ILocation = [-0.055, -0.024, 0.015;
                0.091, -0.057, -0.374;
                0.085, 0.076, -0.025];
            
    ICrossPoints = 3;
    IMIF = 1073;  
    Axis1 = 10;
    Muscle3a = PamData('Iliacus', ILocation, ICrossPoints, IMIF, T, Axis1);
    
    Axis1 = 20;
    Muscle3b = PamData('Iliacus', ILocation, ICrossPoints, IMIF, T, Axis1);
    
    Axis1 = 30;
    Muscle3c = PamData('Iliacus', ILocation, ICrossPoints, IMIF, T, Axis1);
    
    %Include the back to the transformtion matrix for the Psoas
    T(:, :, 2, :) = T(:, :, 1, :);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = Back.TransformationMat(:, :, j);
            ii = ii+1;
        end
    end
    
    %Psoas b1 -> f4
           b1 = [0.036; 0.007; 0.029]; %origin: psoas
    psoas_wr1 = [-0.024; -0.057; 0.076]; %pelvis
           f4 = [0.015; -0.374; -0.025]; %insertion: add_mag (hip adduction), iliacus, psoas. origin: med_gas
           
    PLocation = [0.036, -0.024, 0.015;
                  0.007, -0.057, -0.374;
                  0.029, 0.076, -0.025];
    PCrossPoints = [2 3];
    PMIF = 1113;
    Axis1 = [10 10];
    Muscle4a = PamData('Psoas', PLocation, PCrossPoints, PMIF, T, Axis1);
    
    Axis1 = [20 20];
    Muscle4b = PamData('Psoas', PLocation, PCrossPoints, PMIF, T, Axis1);
    
    Axis1 = [30 30];
    Muscle4c = PamData('Psoas', PLocation, PCrossPoints, PMIF, T, Axis1);
    
    %Torque Calcs, "R" for robot
    HipXTorque = Muscle1a.Torque(1, :) + Muscle2a.Torque(1, :) + Muscle3a.Torque(1, :) + Muscle4a.Torque(2, :);
    HipYTorque = Muscle1b.Torque(1, :) + Muscle2b.Torque(1, :) + Muscle3b.Torque(1, :) + Muscle4b.Torque(2, :);
    HipZTorque = Muscle1c.Torque(1, :) + Muscle2c.Torque(1, :) + Muscle3c.Torque(1, :) + Muscle4c.Torque(2, :);

    %Create the Mesh of Torques to corespond with the joint angles
    HipXTorqueMR = zeros(divisions, divisions); HipYTorqueMR = zeros(divisions, divisions); HipZTorqueMR = zeros(divisions, divisions);
    for i = 1:divisions
        HipXTorqueMR(:, i) = HipXTorque(((i-1)*divisions)+1:i*divisions);
        HipYTorqueMR(:, i) = HipYTorque(((i-1)*divisions)+1:i*divisions);
        HipZTorqueMR(:, i) = HipZTorque(((i-1)*divisions)+1:i*divisions);
    end
    
    %Including Generic variable name for plotting
    RobotAxis1 = HipZ.Theta;
    RobotAxis1Label = 'Hip Flexion, Degrees';
    RobotAxis2 = HipX.Theta;
    RobotAxis2Label = 'Adduction/Abduction, Degrees';
    RobotTorque1 = HipZTorqueMR;
    RobotTorque2 = HipYTorqueMR;
    RobotTorque3 = HipXTorqueMR;
    RobotTitle1 = 'Hip Torque, uniarticular, Z axis';
    RobotTitle2 = 'Hip Torque, uniarticular, Y axis';
    RobotTitle3 = 'Hip Torque, uniarticular, X axis';

end