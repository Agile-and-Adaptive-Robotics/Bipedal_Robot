%% Human Muscle Data
% This script is a collection of all of the human muscles and ways of
% generating torque information for each of those muscles, based on the
% type of actuation they create

addpath Functions

%% Hip Muscles

% Test
iteration = 100;
R = zeros(3, 3, iteration);
T = zeros(4, 4, iteration);
adducMax = 90*pi/180;
testShiftAxis = [1, 0.5, 0];

for i = 0:iteration-1
    theta = adducMax/iteration*i;
    R(:, :, i+1) = [cos(theta), sin(theta), 0;
                    -sin(theta), cos(theta), 0;
                    0, 0, 1];
    
    T(:, :, i+1) = RpToTrans(R(:, :, i+1), testShiftAxis');
end

Name = 'Test Muscle';
MIF = 20;
Location = [0.1, 0.2, 0;
            0.2, 0.3, 0;
            0.3, 0.4, 0];
CrossPoint = 2;
TestM = MuscleData(Name, Location, CrossPoint, MIF, T); 

% Abduction
iteration = 100;
R = zeros(3, 3, iteration);
T = zeros(4, 4, iteration);
adducMax = 30*pi/180;
pelvisToHip = [-0.0707, -0.0661, 0.0835];

% tibiaToPelvis = [


for i = 0:iteration-1
    theta = adducMax/iteration*i;
    R(:, :, i+1) = [cos(theta), sin(theta), 0;
                    -sin(theta), cos(theta), 0;
                    0, 0, 1];
    
    T(:, :, i+1) = RpToTrans(R(:, :, i+1), pelvisToHip');
end

Name = 'Gluteus Maximus 1';
MIF = 573;
Location = [-0.119, 0.061, 0.07;
            -0.129, 0.001, 0.089;
            -0.046, -0.025, 0.039;
            -0.028, -0.057, 0.047];
CrossPoint = 3;
GlutMax1 = MuscleData(Name, Location, CrossPoint, MIF, T); 

Name = 'Sartorius';
MIF = 156;
Location = [-0.015, -0.001, 0.124;
            -0.003, -0.357, -0.042;
            -0.006, -0.042, -0.04;
            0.006, -0.059, -0.038;
            0.024, -0.084, -0.0235];
CrossPoint = [2, 3];
Sar = MuscleData(Name, Location, CrossPoint, MIF, T);

% Name = 'Gracilis';
% Location = [-0.074, -0.119, 0.028;
%             -0.027, -0.032, -0.038;
%             -0.019, -0.052, -0.036;
%             0.006, -0.084, -0.023];
% MIF = 162;
% OFL = 0.352;
% TSL = 0.126;
% Pen = 0.05235988;









