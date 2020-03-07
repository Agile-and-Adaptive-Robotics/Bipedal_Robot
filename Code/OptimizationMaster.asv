% Optimization Master
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This master script will call subscripts in order to 
% optimize muscle attachment locations for the bipedal robot. It first 
% grabs the human torque data and initial robot torque data. Once this is 
% complete, it begin perturbing robot muscle locations until it finds a
% solution that matches human torque data. 

clear
close all
clc

%% ------------- User Entered Parameters -------------
%Choose the Joint that you will to observe
%Options are: Back, Bi_Hip, Calves, Foot, Toe, Uni_Hip
ChooseJoint = 'Back';

%Choose the number of divisions for the angles of rotation
divisions = 100;

%Choose the number of iterations for the optimization code
iterations = 2;

%Choose the minimum change for the value of the location
epsilon = 0.01;

%Choose the scaling factor for the cost function, which weights the
%importance of distance from the attachment point to the nearest point on
%the model body
%Scaling Values
GTorque = 0.0001;           %Cost weight for the difference between human and robot torque      
GDiameter40 = 1e4;              %Cost Weight for the diameter of the festo muscle
GDiameter20 = 1e2;
G = 1000;                   %Cost weight for the distance from the attacment point to the model body
GLength = 1000;

%Adjust the axis range for the Torque plots
caxisRange = [-40 150];

%Adjust the axis for the robot model plot
axisLimits = [-1 1 -1 1 -1.25 0.75];

%% ----------------- Setup -------------------------------
%Include relevant folders
addpath('Open_Sim_Bone_Geometry')
addpath('Functions')
addpath('Human_Data')

%% ------------- Humanoid Model --------------
%Runs the humanoid model. Only run if you need to update the data
%run("HumanoidMuscleCalculation.m")

%Loads important data for the human model. Data was previously created and
%then stored, for quicker load time.
%Data saved: Back, Bi_Hip, Calves, Foot, Toe, Uni_Hip
load(strcat('Human_', ChooseJoint, '_Data.mat'));

%% ------------- Robot Model -----------------
%Runs the bipedal model for initial calculations and creating the proper
%classes
run("RobotPAMCalculationOptimization.m")    

%% ------------- Optimization ------------------
%The robot model has now constructed a preliminary model and generated
%torques. We will now begin to move around the attachment points to see if
%we can generate better torques

%Load in points that construct the relevant bone models to the muscles.
%This will be used as a part of our cost function, to evaluate how far away
%our new muscle attachments are. 
if isequal(ChooseJoint, 'Back')
    PointsFile = 'Pelvis_R_Mesh_Points.xlsx';
    Pelvis = xlsread(PointsFile)';
    
    PointsFile = 'Spine_Mesh_Points.xlsx';
    Home = Joint1a.Home;
    Spine = xlsread(PointsFile)'+Home;
    
    PointsFile = 'Sacrum_Mesh_Points.xlsx';
    Sacrum = xlsread(PointsFile)';
end

%Determine how many muscles are included in the algorithm
MuscleNum = size(Muscles, 2);

%Calculate the distance from the attachment point to the nearest point on
%the robot body
%The first step is to get all the attachment points into the same
%reference frame. 

L = Muscle1.Location;
iii = 0;                %variable that remembers which joint we are looking at. Will update to 1 for the first joint, 2 for the second cross over joint, and so on.
for k = 1:size(Muscle1.Location, 2)
    for ii = 1:size(Muscle1.CrossPoints, 2)
        iii = iii+1;
        if k == Muscle1.CrossPoints(ii)
            L(:, k) = L(:, k)+Muscle1.TransformationMat(1:3, 4, iii);
        end
    end
end

%Create a tensor to store all the points that describe the bone mesh
bone{1} = Spine;
bone{2} = Sacrum;
bone{3} = Pelvis;

%Evaluate cost function for the initial set of attachment points
C(1) = 0;
k = 1;

for ii = 1:100
    for iii = 1:100
        C(1) = C(1) + abs(HumanTorque1(ii, iii) - RobotTorque1(ii, iii));
        C(1) = C(1) + abs(HumanTorque2(ii, iii) - RobotTorque2(ii, iii));
        if exist('RobotTorque3', 'var') == 1
            C(1) = C(1) + abs(HumanTorque3(ii, iii) - RobotTorque3(ii, iii));
            if exist('RobotTorque4', 'var') == 1
                C(1) = C(1) + abs(HumanTorque4(ii, iii) - RobotTorque4(ii, iii));
            end
        end
    end
end

%For this part of the cost function, we must sum across all muscles
for i = 1:MuscleNum
    Diameter = Muscles{i}.Diameter;
    MLength = Muscles{i}.MuscleLength;
    %Increase the cost based on length of muscle
    for ii = 1:length(MLength)
        C(k) = C(k) + GLength*MLength(ii);
    end

    %Increase the cost based on the diameter of the muscle
    if Diameter == 40
        C(k) = C(k) + GDiameter40;
    elseif Diameter == 20
        C(k) = C(k) + GDiameter20;
    end
    MLength = [];
end

beginOptimization = 1;          %Flags that optimization has begun for the optimization pam calculations
%Perturb original location by adding epsilon to every cross point
%Will then move to creating an algorith that will add and subtract to
%different axes

%Start the optimization as a for loop to create all the points of a cube.
%later, automate this to be in a while loop that go until a threshold
ep1 = zeros(3, 1);              %Change in the first cross point for the Muscle
ep2 = zeros(3, 1);              %Change in the second cross point for the muscle
neg = [1, -1];
k = 2;                  %index for the cost function. Start at 2, since 1 was the original point
shrinkEp = 0;                   %Flag. Once it reaches the number of muscles, it will scale epsilon down to refine the search

%We begin with the point before the crossing point
for m = 1:MuscleNum
    Cross(m) = Muscles{m}.CrossPoints(1);  %Will need to change this indexing for muscles that have multiple crossing points
    StartingLocation1{m} = Muscles{m}.Location(:, Cross(m) - 1);
    StartingLocation2{m} = Muscles{m}.Location(:, Cross(m));
    LocationTracker1{m} = StartingLocation1{m};
    LocationTracker2{m} = StartingLocation2{m};
end

disp('Beginning Optimization');
for iiii = 1:iterations
    for m = 1:MuscleNum
        for k1 = 1:2
            ep1(1) = epsilon*neg(k1);
            for k2 = 1:2
                ep1(2) = epsilon*neg(k2);
                for k3 = 1:2
                    ep1(3) = epsilon*neg(k3);
                    for k4 = 1:2
                        ep2(1) = epsilon*neg(k4);
                        for k5 = 1:2
                            ep2(2) = epsilon*neg(k5);
                            for k6 = 1:2
                                ep2(3) = epsilon*neg(k6);

                                %Change the Locaiton of the first and
                                %second crossing point
                                Location{m}(:, Cross(m) - 1) = StartingLocation1{m} + ep1;
                                Location{m}(:, Cross(m)) = StartingLocation2{m} + ep2;
                                
                                %Update all of the Location Tracker
                                %variables to identify current
                                %configuration
                                for config = 1:MuscleNum
                                    LocationTracker1{config}(:, k) = Location{config}(:, Cross(config) - 1);
                                    LocationTracker2{config}(:, k) = Location{config}(:, Cross(config));
                                end
                                
                                run('RobotPAMCalculationOptimization.m')

                                %Evaluate cost function for the new set of attachment points
                                C(k) = 0;

                                for ii = 1:100
                                    for iii = 1:100
                                        C(k) = C(k) + abs(HumanTorque1(ii, iii) - RobotTorque1(ii, iii));
                                        C(k) = C(k) + abs(HumanTorque2(ii, iii) - RobotTorque2(ii, iii));
                                        if exist('RobotTorque3', 'var') == 1
                                            C(k) = C(k) + abs(HumanTorque3(ii, iii) - RobotTorque3(ii, iii));
                                            if exist('RobotTorque4', 'var') == 1
                                                C(k) = C(k) + abs(HumanTorque4(ii, iii) - RobotTorque4(ii, iii));
                                            end
                                        end
                                    end
                                end

                                %Increase the cost based on length of muscle
                                for ii = 1:length(Muscles{m}.MuscleLength)
                                    C(k) = C(k) + GLength*Muscles{m}.MuscleLength(ii);
                                end

                                %Increase the cost based on the diameter of the muscle
                                if Muscles{m}.Diameter == 40
                                    C(k) = C(k) + GDiameter40;
                                elseif Muscles{m}.Diameter == 20
                                    C(k) = C(k) + GDiameter20;
                                end

                                disp(['Iteration number ', num2str(k-1), ' out of of ', num2str(2^6*iterations*MuscleNum), '.']);
                                k = k+1;            %Increment k for the next point of the cost function
                            end
                        end
                    end
                end
            end
        end
%         [a, b] = min(C(end-2^6:end));
        [a, b] = min(C);

        %If the previous starting location is still the location with the
        %minimum cost, change the step size of epsilon.
        if StartingLocation1{m} == LocationTracker1{m}(:, b)
            if StartingLocation2{m} == LocationTracker2{m}(:, b)
                shrinkEp = shrinkEp+1;
                if shrinkEp == MuscleNum
                    shrinkEp = 0;
                    epsilon = 0.1*epsilon;
                    if epsilon < 10^-6
                        break
                    end
                end
            end
        end
        StartingLocation1{m} = LocationTracker1{m}(:, b);               %Update the next round of optimization with the location with minimum cost
        StartingLocation2{m} = LocationTracker2{m}(:, b);
    end
end

%%------------------ Plotting ----------------------
% run('OptimizationPlotting.m')

figure
plot(C)
xlabel('Epochs')
ylabel('Cost Value')