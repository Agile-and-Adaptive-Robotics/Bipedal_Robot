%% Human Muscle Data
% This script is a collection of all of the human muscles and ways of
% generating torque information for each of those muscles, based on the
% type of actuation they create
clear
clc
close all


addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions

%% ------------- Test Case ----------------
%% Test Monoarticulate Muscle
iteration = 20;

adducMax = -90*pi/180;
adducMin = 0;
theta = linspace(adducMin, adducMax, iteration);

flexMax = 180*pi/180;
extMax = 0;
phi = linspace(extMax, flexMax, iteration);

testShiftAxis = [1, 0.5, 0];

R = zeros(3, 3, iteration);
Rx = zeros(3, 3, iteration);
Ry = zeros(3, 3, iteration);
Rz = zeros(3, 3, iteration);
T = zeros(4, 4, iteration);

pos = 1;
for ii = 1:size(theta, 2)
    for i = 1:size(phi, 2)
        Rz(:, :, i) = [cos(theta(i)), -sin(theta(i)), 0;
                        sin(theta(i)), cos(theta(i)), 0;
                        0, 0, 1];
                    
        Rx(:, :, ii) = [1 0 0;
                        0 cos(phi(ii)), -sin(phi(ii));
                        0 sin(phi(ii)), cos(phi(ii))];
        
        R(:, :, pos) = Rz(:, :, i)*Rx(:, :, ii);            
                    
        T(:, :, pos) = RpToTrans(R(:, :, pos), testShiftAxis');
        pos = pos + 1;
    end
end

Name = 'Test Muscle';
MIF = 20;
Location = [0.1, 0.2, 0;
            0.4, 0.5, 0];
CrossPoint = 2;
TestM = MonoMuscleData(Name, Location, CrossPoint, MIF, T); 

%% Test Biarticulate Muscle with attachment point after the middle joint
clear T
iteration = 2;

adducMax = -90*pi/180;
adducMin = 0;
theta = linspace(adducMin, adducMax, iteration);

testShiftAxis = [1, 0, 0];

R = zeros(3, 3, iteration);
Rx = zeros(3, 3, iteration);
Ry = zeros(3, 3, iteration);
Rz = zeros(3, 3, iteration);
T = zeros(4, 4, iteration);

for i = 1:size(theta, 2)
    Rz(:, :, i) = [cos(theta(i)), -sin(theta(i)), 0;
                    sin(theta(i)), cos(theta(i)), 0;
                    0, 0, 1];

    R(:, :, i) = Rz(:, :, i);            

    T(:, :, i) = RpToTrans(R(:, :, i), testShiftAxis');
end
T(:, :, :, 2) = T(:, :, :, 1);
Name = 'Test Biarticular Muscle';
MIF = 20;
Location = [0.8, 0.7, 0;
            0.5, 0.4, 0];
CrossPoint = [2, 2];
TestBM = BiMuscleData(Name, Location, CrossPoint, MIF, T);

%Verification Calculation
v1 = Location(1, :);
v2 = Location(2, :);
v1p = RowVecTrans(T(:, :, 1, 1)\eye(4), v1);
v1pp = RowVecTrans((T(:, :, 1, 2)*T(:, :, 1, 1))\eye(4), v1);
v2p = RowVecTrans(T(:, :, 1, 2), v2);

TestBM.UnitDirection

unitD1 = (v1p - v2p)/norm(v1p - v2p)
unitD2 = (v1pp - v2)/norm(v1pp - v2)

TestBM.MomentArm

mA1 = v2p - unitD1*dot(unitD1, v2p)
mA2 = v2 - unitD2*dot(unitD2, v2)