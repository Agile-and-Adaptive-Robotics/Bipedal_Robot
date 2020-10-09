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
ChooseJoint = 'Bi_Hip';

%Choose the number of positions for the angles of rotation
positions = 100;

%Choose the number of iterations for the optimization code 
iterations = 1;

    %Note: Time calculation equation (estimates how long this will take)
%     totalTimeSeconds = 3*iterations*2^6*3
%     totalTimeMinutes = totalTimeSeconds/60
%     totalTimeHours = totalTimeMinutes/60

%Choose the minimum change for the value of the location
epsilon = 0.1;
refinementRate = 0.5;       %The gain on epsilon once a local minima has been found
breakCheck = 1e-3;          %Parameter that determines if epsilon has been reduced enough to end the optimization process

%Choose the scaling factor for the cost function, which weights the
%importance of distance from the attachment point to the nearest point on
%the model body
%Scaling Values
% GTorque = 1e-4;           %Cost weight for the difference between human and robot torque      
% GDiameter40 = 1e7;              %Cost Weight for the diameter of the festo muscle
% GDiameter20 = 1e3;
% G = 1000;                   %Cost weight for the distance from the attachment point to the model body
% GLength = 100;

%New approach to setting gains. Going to try to have all of them sum to 1.
GTorque = 1e-4;
GLength = 1e-2;
disG = 500;               %Cost weight for the distance away from the starting point

GDiameter20 = 1e0;
GDiameter40 = 5e0;
G = 1e-5;

%Adjust the axis range for the Torque plots
caxisRange = [-40 150];

%Adjust the axis for the robot model plot
axisLimits = [-1 1 -1 1 -1.25 0.75];

%% ----------------- Setup -------------------------------
%Include relevant folders
addpath('Open_Sim_Bone_Geometry')
addpath('Functions')
addpath('Human_Data')
addpath('Joint_Optimization_Scripts')

%% ------------- Humanoid Model --------------
%Loads important data for the human model. Data was previously created and
%then stored, for quicker load time.
%Data saved: Back, Bi_Hip, Calves, Foot, Toe, Uni_Hip
load(strcat('Human_', ChooseJoint, '_Data.mat'));

%% ------------- Robot Model -----------------
%Runs the bipedal model for initial calculations and creating the proper
%classes
JointOptimizationScript = strcat('RobotPAMCalculationOptimization_', ChooseJoint);
run(JointOptimizationScript)

%Set up the original robot torque if I want to see if the error is better
OriginalRobotTorque1 = RobotTorque1;
OriginalRobotTorque2 = RobotTorque2;
if exist('RobotTorque3', 'var') == 1
    OriginalRobotTorque3 = RobotTorque3;
    if exist('RobotTorque4', 'var') == 1
        OriginalRobotTorque4 = RobotTorque4;
    end
end


%% ------------- Optimization ------------------
%The robot model has now constructed a preliminary model and generated
%torques. We will now begin to move around the attachment points to see if
%we can generate better torques

%Evaluate cost function for the initial set of attachment points
% C = zeros(1, MuscleNum*2^6*iterations);   %Don't do this. Need to look at
% the min later on, and initializing with zeros fucks it up. Maybe init
% with 1's and then multiply by a lot
eC(1) = 0;              %Error component of the cost function
mC(1) = 0;              %Muscle length component of the cost function
dC(1) = 0;              %Muscle diameter component of the cost function
C(1) = 0;               %Cost function value

%Calculate error for original robot and the human model
for ii = 1:100
    for iii = 1:100
        eC(1) = eC(1) + GTorque*abs(HumanTorque1(ii, iii) - RobotTorque1(ii, iii));
        eC(1) = eC(1) + GTorque*abs(HumanTorque2(ii, iii) - RobotTorque2(ii, iii));
        if exist('RobotTorque3', 'var') == 1
            eC(1) = eC(1) + GTorque*abs(HumanTorque3(ii, iii) - RobotTorque3(ii, iii));
            if exist('RobotTorque4', 'var') == 1
                eC(1) = eC(1) + GTorque*abs(HumanTorque4(ii, iii) - RobotTorque4(ii, iii));
            end
        end
    end
end

%Calculate the length and diameter components of the cost function
MuscleNum = size(Muscles, 2);       %Keep track of how many muscles are involved in this joint.
k = 1;                              %Indexing variable for each iteration of the optimization. The original robot sets the starting point
for i = 1:MuscleNum
    Diameter = Muscles{i}.Diameter;
    MLength = Muscles{i}.MuscleLength;
    %Increase the cost based on length of muscle
    for ii = 1:length(MLength)
        mC(k) = mC(k) + GLength*MLength(ii);
    end

    %Increase the cost based on the diameter of the muscle
    if Diameter == 40
        dC(k) = dC(k) + GDiameter40;
    elseif Diameter == 20
        dC(k) = dC(k) + GDiameter20;
    end
   
    MLength = [];
end
C(1) = eC(1) + mC(1) + dC(1);

beginOptimization = 1;          %Flags that optimization has begun for the optimization pam calculations
%Perturb original location by adding epsilon to every cross point
%Will then move to creating an algorith that will add and subtract to
%different axes

%Start the optimization process. Continues to calculate cost function for a
%given amount of iterations or until epsilon becomes a negligible value

ep1 = zeros(3, 1);              %Change in the first cross point for the Muscle
ep2 = zeros(3, 1);              %Change in the second cross point for the muscle
neg = [1, -1];
k = 2;                          %index for the cost function. Start at 2, since 1 was the original point
myBreak = 0;                    %Flag for if epsilon becomes extremely small, in which case the optimization should end
previousBestC = C(1);           %Variable that stores the best value from the previous series of calculations. Compares to newest best ot determine if actions need to be takne
previousBestIteration = 1;      %Variable that stores the iteration number of the lowest cost value

%Initializing a location tracker and starting locations in order to observe
%changes once the optimization process has completed.

LocationTracker1 = cell(MuscleNum, 1);
LocationTracker2 = cell(MuscleNum, 1);
StartingLocation1 = cell(MuscleNum, 1);
StartingLocation2 = cell(MuscleNum, 1);
for m = 1:MuscleNum
    LocationTracker1{m} = zeros(3, MuscleNum*2^6*iterations, size(Muscles{1}.CrossPoints, 2));   %Tracks location of point prior to a joint
    LocationTracker2{m} = zeros(3, MuscleNum*2^6*iterations, size(Muscles{1}.CrossPoints, 2));   %Tracks location of point after a joint
    StartingLocation1{m} = zeros(3, size(Muscles{1}.CrossPoints, 2));           %Updates with the optimization starting location for the attachment point BEFORE a joint
    StartingLocation2{m} = zeros(3, size(Muscles{1}.CrossPoints, 2));           %Updates with the optimization starting location for the attachment point AFTER a joint
end
Cross = zeros(MuscleNum, size(Muscles{1}.CrossPoints, 2));            %Create a variable that tracks all of the crossing points for the muscles

%We begin with the point before the crossing point
for m = 1:MuscleNum
    Cross(m, :) = Muscles{m}.CrossPoints;  
    for i = 1:size(Cross, 2)
        if Cross(m, i) == 0
            i = i + 1;
        end
        StartingLocation1{m}(:, i) = Muscles{m}.Location(:, Cross(m, i) - 1);
        StartingLocation2{m}(:, i) = Muscles{m}.Location(:, Cross(m, i));
        LocationTracker1{m}(:, 1, i) = StartingLocation1{m}(:, i);
        LocationTracker2{m}(:, 1, i) = StartingLocation2{m}(:, i);
    end
end

tic
disp('Beginning Optimization');
for iiii = 1:iterations
    for m = 1:MuscleNum
        %Repeat 8 iterations for a single muscle looking at different
        %locations for the muscle cross points
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

            %Change the Location of the first and
            %second crossing point
            for i = 1:size(Cross, 2)
                if Cross(m, i) == 0
                    i = i+1;
                end
                Location{m}(:, Cross(m, i) - 1) = StartingLocation1{m}(:, i) + ep1;
                Location{m}(:, Cross(m, i)) = StartingLocation2{m}(:, i) + ep2;
            end
            %Update all of the Location Tracker
            %variables to identify current
            %configuration
            for config = 1:MuscleNum
                for i = 1:size(Cross, 2)
                    if Cross(config, i) == 0
                        i = i+1;
                    end
                    LocationTracker1{config}(:, k, i) = Location{config}(:, Cross(config, i) - 1);
                    LocationTracker2{config}(:, k, i) = Location{config}(:, Cross(config, i));
                end
            end

            run(JointOptimizationScript)

            %Evaluate cost function for the new set of attachment points
            eC(k) = 0;
            mC(k) = 0;
            dC(k) = 0;
            disC(k) = 0;
            C(k) = 0;

            for ii = 1:positions
                for iii = 1:positions
                    if(HumanTorque1(ii, iii) > RobotTorque1(ii, iii))
                        eC(k) = eC(k) + GTorque*abs(HumanTorque1(ii, iii) - RobotTorque1(ii, iii));
                    end
                    if(HumanTorque2(ii, iii) > RobotTorque2(ii, iii))
                        eC(k) = eC(k) + GTorque*abs(HumanTorque2(ii, iii) - RobotTorque2(ii, iii));
                    end
                    if exist('RobotTorque3', 'var') == 1
                        if(HumanTorque3(ii, iii) > RobotTorque3(ii, iii))
                            eC(k) = eC(k) + GTorque*abs(HumanTorque3(ii, iii) - RobotTorque3(ii, iii));
                        end
                        if exist('RobotTorque4', 'var') == 1
                            if(HumanTorque4(ii, iii) > RobotTorque4(ii, iii))
                                eC(k) = eC(k) + GTorque*abs(HumanTorque4(ii, iii) - RobotTorque4(ii, iii));
                            end
                        end
                    end
                end
            end

            %Increase the cost based on length of each muscle
            for ii = 1:MuscleNum
                for iii = 1:length(Muscles{ii}.MuscleLength)
                    mC(k) = mC(k) + GLength*Muscles{ii}.MuscleLength(iii);
                end
            end

            %Increase the cost based on the diameter of the muscle
            for ii = 1:MuscleNum
                if Muscles{ii}.Diameter == 40
                    dC(k) = dC(k) + GDiameter40;
                elseif Muscles{ii}.Diameter == 20
                    dC(k) = dC(k) + GDiameter20;
                end
            end

            %Increase the cost based on how far way the
            %new placement is from the original
            %Want to change this to be distance from the bone surface
            for ii = 1:MuscleNum
                for i = 1:size(Cross, 2)
                    if Cross(ii, i) == 0
                        i = i+1;
                    end
                    disC(k) = disC(k) + disG*norm(LocationTracker1{ii}(:, k, i) - LocationTracker1{ii}(:, 1, i))^2;
                    disC(k) = disC(k) + disG*norm(LocationTracker2{ii}(:, k, i) - LocationTracker2{ii}(:, 1, i))^2;
                end
            end

            C(k) = eC(k) + mC(k) + dC(k)+disC(k);

            %Keep track of robot torque for error plots
            %later
            RTorqueTracker1(:, :, k) = RobotTorque1(:, :);
            RTorqueTracker2(:, :, k) = RobotTorque2(:, :);
            if exist('RobotTorque3', 'var') == 1
                RTorqueTracker3(:, :, k) = RobotTorque3(:, :);
                if exist('RobotTorque4', 'var') == 1
                    RTorqueTracker4(:, :, k) = RobotTorque4(:, :);
                end
            end

            disp(['Iteration number ', num2str(k-1), ' out of of ', num2str(2^6*iterations*MuscleNum), '.']);
            
            k = k+1;            %Increment k for the next point of the cost function

        end
        end
        end
        end
        end
        end
        %Once the best position is chosen for one muscles, should update
        %its location to the current best when starting on the other
        %muscles
        [bestMuscleC, bestMuscleIteration] = min(C(end-63:end));                        %Check for the best iteration of the previous cycle for one muscle
        for i = 1:size(Cross, 2)
            if Cross(m, i) == 0
                i = i+1;
            end
            Location{m}(:, Cross(m, i) - 1) = LocationTracker1{m}(:, bestMuscleIteration, i);     %For the following muscles, use the best version of the previous calculations
            Location{m}(:, Cross(m, i)) = LocationTracker2{m}(:, bestMuscleIteration, i);
        end
    end
    [currentBestC, currentBestIteration] = min(C);
    if currentBestIteration == previousBestIteration
        epsilon = epsilon*refinementRate;
        if epsilon < breakCheck
            myBreak = 1;
        end
    end

    for n = 1:MuscleNum
        for i = 1:size(Cross, 2)
            if Cross(n, i) == 0
                i = i+1;
            end
            StartingLocation1{n}(:, i) = LocationTracker1{n}(:, currentBestIteration, i);               %Update the next round of optimization with the location with minimum cost
            StartingLocation2{n}(:, i) = LocationTracker2{n}(:, currentBestIteration, i);
        end
    end
 
    previousBestIteration = currentBestIteration;
    previousBestC = currentBestC;
    
    if myBreak == 1             %Likely should change the loop to be a while loop for the iterations. Can tackle that in another branch soon
        break
    end
end
timeElapsed = toc;

NewRobotTorque1 = RTorqueTracker1(:, :, currentBestIteration);
NewRobotTorque2 = RTorqueTracker2(:, :, currentBestIteration);
if exist('RobotTorque3', 'var') == 1
    NewRobotTorque3 = RTorqueTracker3(:, :, currentBestIteration);
    if exist('RobotTorque4', 'var') == 1
        NewRobotTorque4 = RTorqueTracker4(:, :, currentBestIteration);
    end
end


%%------------------ Plotting ----------------------
% run('OptimizationPlotting.m')

figure
hold on
plot(C, 'k')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')
xlim([0, length(C)])
ylim([0, 800])
set(gca, 'FontSize', 12)
hold off

%Create an average of the iterationst to create viewable epochs
for i = 1:iterations
    startP = (i-1)*MuscleNum*2^6+1;             %Starting Point for summation of the epoch
    endP = i*MuscleNum*2^6;                     %Ending poster for summation of the epoch
    if i*MuscleNum*2^6 < length(C)
        averageC(i) = sum(C(startP:endP));
    else
        averageC(i) = sum(C(startP:end));
    end
end

%Possible GitHub bug. Weird artifacting in my code

% 
% figure
% plot(averageC)
% xlabel('Epochs')
% ylabel('Cost Value')
% <<<<<<< HEAD
% <<<<<<< HEAD
% xlim([0 15])
% =======
% =======
% >>>>>>> 6bc9ba6699a00a8adb68bb672884ee8edc838d8d
% title('Average C')
% 
% figure
% subplot(2, 2, 1)
% plot(eC)
% title('Error Component')
% 
% subplot(2, 2, 2)
% plot(dC)
% title('Diameter Component')
% 
% subplot(2, 2, 3)
% plot(mC)
% title('Muscle Length Component')
% 
% subplot(2, 2, 4)
% plot(disC)
% title('Distance Component')
% 
% %Mean Squared Error for original robot placement
% oMSE = 0; 
% for ii = 1:100
%     for iii = 1:100
%         oMSE = oMSE + (HumanTorque1(ii, iii) - OriginalRobotTorque1(ii, iii))^2;
%         oMSE = oMSE + (HumanTorque2(ii, iii) - OriginalRobotTorque2(ii, iii))^2;
%         if exist('RobotTorque3', 'var') == 1
%             oMSE = oMSE + (HumanTorque3(ii, iii) - OriginalRobotTorque3(ii, iii))^2;
%         end
%         if exist('RobotTorque4', 'var') == 1
%             oMSE = oMSE + (HumanTorque4(ii, iii) - OriginalRobotTorque4(ii, iii))^2;
%         end
%     end
% end
% oMSE = oMSE/divisions^2;
% 
% %Mean squared error for new placements
% MSE = 0; 
% for ii = 1:100
%     for iii = 1:100
%         MSE = MSE + (HumanTorque1(ii, iii) - NewRobotTorque1(ii, iii))^2;
%         MSE = MSE + (HumanTorque2(ii, iii) - NewRobotTorque2(ii, iii))^2;
%         if exist('RobotTorque3', 'var') == 1
%             MSE = MSE + (HumanTorque3(ii, iii) - NewRobotTorque3(ii, iii))^2;
%         end
%         if exist('RobotTorque4', 'var') == 1
%             MSE = MSE + (HumanTorque4(ii, iii) - NewRobotTorque4(ii, iii))^2;
%         end
%     end
% end
% <<<<<<< HEAD
% MSE = MSE/divisions^2;
% >>>>>>> 6bc9ba6699a00a8adb68bb672884ee8edc838d8d
% =======
% MSE = MSE/divisions^2;
% >>>>>>> 6bc9ba6699a00a8adb68bb672884ee8edc838d8d
