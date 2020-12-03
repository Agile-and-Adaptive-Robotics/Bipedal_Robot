%% BPA Data Test
% This script sets up test BPAs to validate the calculations in the BPA
% classes
clear
clc
close all


addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions

%% ------------- Test Case ----------------
%% Test Monoarticulate Muscle
iteration = 2;

adducMax = -90*pi/180;
adducMin = 0;
theta = linspace(adducMin, adducMax, iteration);

flexMax = 180*pi/180;
extMax = 0;
phi = linspace(extMax, flexMax, iteration);

testShiftAxis = [1, 0, 0];

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
Location = [0.5, 0.5, 0;
            0.5, 0.5, 0;
            0.5, -0.6, 0];
CrossPoint = 2;
Dia = 10;
TestBPA = MonoPamData(Name, Location, CrossPoint, Dia, T); 

% %Position 1 verification
% v1 = Location(1, :);
% v2 = Location(2, :);
% v3 = Location(3, :);
% 
% mL1 = norm(v1 - RowVecTrans(T(:, :, 1), v2));
% mL2 = norm(v2 - v3);
% mL = mL1 + mL2;
% 
% d = RowVecTrans(T(:, :, 1)\eye(4), v1) - v2;
% u = d/norm(d);
% 
% TestBPA.MuscleLength(1)
% 
% TestBPA.UnitDirection(1, :)
% 
% TestBPA.MomentArm(1, :)
% 
% %Position 2 verification
% mL1 = norm(v1 - RowVecTrans(T(:, :, 3), v2))
% 
% d = RowVecTrans(T(:, :, 3)\eye(4), v1) - v2
% u = d/norm(d)
% 
% mA = v2 - u*dot(u, v2)
% 
% TestBPA.MuscleLength(3)
% 
% TestBPA.UnitDirection(3, :)
% 
% TestBPA.MomentArm(3, :)


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
Location = [0.5, 0.5, 0;
            0.5, 0.5, 0;
            1, -0.5, 0];
CrossPoint = [2, 2];
Dia = 20;
TestBM = BiPamData(Name, Location, CrossPoint, Dia, T);

%Verification Calculation, home position
v1 = Location(1, :);
v2 = Location(2, :);
v3 = Location(3, :);

mL2 = sqrt((1-0.5)^2 + (-1)^2);

%Verification Calculation, 1st joint rotation
v22 = [1.5, 0.5, 0];
u1 = [-2, -1, 0]/norm([-2, -1, 0]);
u2 = u1;

mA1 = v22 - u1*dot(u1, v22);
mL1 = norm([2, 1, 0])

%Verification Calculation 2nd joint rotation
v22 = [1.5, -0.5, 0];
u1 = [-2, 1, 0]/norm([-2, 1, 0])

u2 = [-1, -2, 0]/norm([-1, -2, 0])

mL1 = norm([-1, -2, 0])



%% Test Biarticulate, with points on the first and second joint
% Name = 'Second Test Biarticular Muscle';
% MIF = 20;
% Location = [0.8, 0.7, 0;
%             0.5, 0.4, 0;
%             0.9, 0.6, 0];
% CrossPoint = [2, 3];
% TestBM2 = BiMuscleData(Name, Location, CrossPoint, MIF, T);     
% 
% %Verification Calculation
% v1 = Location(1, :);
% v2 = Location(2, :);
% v3 = Location(3, :);
% 
% TestBM2.MuscleLength
% 
% mL1 = norm(v1 - RowVecTrans(T(:, :, 2, 1), v2))
% mL2 = norm(v2 - RowVecTrans(T(:, :, 2, 2), v3))
% 
% mL = mL1 + mL2
% 
% v1p = RowVecTrans(T(:, :, 2, 1)\eye(4), v1);
% v2p = RowVecTrans(T(:, :, 2, 2)\eye(4), v2);
% 
% TestBM2.UnitDirection
% 
% unitD1 = (v1p - v2)/norm(v1p - v2)
% unitD2 = (v2p - v3)/norm(v2p - v3)
% 
% TestBM2.MomentArm
% 
% mA1 = v2 - unitD1*dot(unitD1, v2)
% mA2 = v3 - unitD2*dot(unitD2, v3)