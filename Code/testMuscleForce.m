% This script is for playing with the muscle force that can be generated,
% using Hoy, Thelen, and Millard as references

%% Freshen up the workspace
clc
clear
close all

%% Add paths to the muscle and pam calculators
addpath Human_Data
addpath Robot_Data
addpath Functions


%% Joint rotation transformation matrices
positions = 20;
Rx = zeros(3, 3, positions);
Ry = zeros(3, 3, positions);
Rz = zeros(3, 3, positions);
T = zeros(4, 4, positions);

% Adduction/Abduction, X rotation
adducMax = 20*pi/180;
adducMin = -45*pi/180;
theta = linspace(adducMin, adducMax, positions);

% Internal/External rotation, Y rotation
rotationMax = 40*pi/180;
rotationMin = -45*pi/180;
phi = linspace(rotationMin, rotationMax, positions);

% Flexion/Extension, Z rotation
flexMax = 85*pi/180;
flexMin = -25*pi/180;
gamma = linspace(flexMin, flexMax, positions);

pelvisToHip = [-0.0707, -0.0661, 0.0835];

j = 1;
for iii = 1:length(gamma)
    for ii = 1:length(phi)
        for i = 1:length(theta)
            Rx(:, :, i) = [1, 0, 0;
                           0, cos(theta(i)), -sin(theta(i));
                           0, sin(theta(i)), cos(theta(i))];
                        
            Ry(:, :, ii) = [cos(phi(ii)), 0, -sin(phi(ii));
                            0, 1, 0;
                            sin(phi(ii)), 0, cos(phi(ii))];
                
            Rz(:, :, iii) = [cos(gamma(iii)), -sin(gamma(iii)), 0;
                            sin(gamma(iii)), cos(gamma(iii)), 0;
                            0, 0, 1];
                        
            R = Rz(:, :, iii)*Ry(:, :, ii)*Rx(:, :, i);
            T(:, :, j) = RpToTrans(R, pelvisToHip');
            j = j + 1;
        end
    end
end


%% Muscle calculation

Name = 'Adductor Magnus 2';
MIF = 343;
Location = [-0.083, -0.119, 0.031;
            0.005, -0.229, 0.023];
CrossPoint = 2;
Add_Mag2 = MonoMuscleData(Name, Location, CrossPoint, MIF, T);

%% Force Calculation 

gamma = 0.5;                    %Shape factor for force length curve of contractile element
epsilon = 0.6;                  %Passive muscle strain due to maximum isometric force
kPE = 4;                        %Exponential shape facor for passive muscle element

tsL = 0.12;                     %Tendon slack length
alphaP = 0.05235988;            %Optimal pennation angle
ofL = 0.121;                    %Optimal muscle fiber length

mtL = Add_Mag2.MuscleLength(1)                   %Muscle-tendon length

mif = 343;                      %maximum isometric force

%Estimate the nominal muscle length, without the tendon
if mtL == tsL
    mtL0 = (mtL*1.001 - tsL)/ofL;
else
    mtL0 = (mtL - tsL)/ofL;
end

%Set function solver parameters
options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);

mL = fsolve(@fBalance, [mtL0, kPE, epsilon, gamma, alphaP, tsL], options);

function muscleF = fBalance(nmL)
    fPE = (exp(kPE*(nmL - 1)/epsilon - 1)/(exp(kPE) - 1));
    
    fL = exp(-(nmL - 1).^2/gamma);

    cosAlpha = sqrt(1 - (sin(alphaP)/nmL)^2);
    
    fT = 37.5/tsL*(mtL - nmL*cosAlpha - tsL);
    
    muscleF = (fL + fPE)*cosAlpha - fT;
end