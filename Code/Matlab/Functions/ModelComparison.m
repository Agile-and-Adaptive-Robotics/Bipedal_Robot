%% Compare Martens and Sarosi models to ours
% Compare at arbitrary lrest of 150 mm, 300 mm, and 450 mm
clear; clc; close all

lrest = [.15, .3, .45];
%% Load our data
load allData.mat XX YY ZZ xTrain yTrain zTrain xVal yVal zVal
load data20mm_sorted.mat Ax Ay Az xTrain20 yTrain20 zTrain20 xVal20 yVal20 zVal20
load FestoLookup.mat f_10 f20

%% Redimensionalize force
%Calculate max force for resting lengths
Fmax10 = maxBPAforce(lrest,'10',620); %10 mm
Fmax20 = maxBPAforce(lrest,'20',620); %20 mm

%Redimensionalize force
Z10 = Fmax10'*ZZ';
Z20 = Fmax20'*Az';

%Sort it for the curve fitter
Z10_150 = Z10(:,1); %10 mm
Z10_300 = Z10(:,2);
Z10_450 = Z10(:,3);
Z20_150 = Z20(:,1); %20 mm
Z20_300 = Z20(:,2);
Z20_450 = Z20(:,3);

%% Our model predicted values


%% Martens and Boblan coefficients
%Coeff.
% Converted from [MPa, MPA/m, GPa/m^2, MPa/m^3, deg(?)]
c10 = [74.085*10^3 -689.20*10^3 1.837*10^6  -848.79*10^3 -0.1]; %10 mm
c20 = [93.232*10^3 -715.29*10^3 1.5483*10^6 -502.95*10^3 -0.073]; %20 mm

Mar1 = Martens(lrest(1),10,c10,XX,YY);

%% Sarosi coefficients
a10 = [-9.2194029;
a20 = [1 ];

%% Our model
B10_150 = f_10(XX,YY).*maxBPAforce(lrest(1),'10',620);
B10_300 = f_10(XX,YY).*maxBPAforce(lrest(2),'10',620);
B10_450 = f_10(XX,YY).*maxBPAforce(lrest(3),'10',620);

B20_150 = f20(XX,YY).*maxBPAforce(lrest(1),'20',620);
B20_300 = f20(XX,YY).*maxBPAforce(lrest(2),'20',620);
B20_450 = f20(XX,YY).*maxBPAforce(lrest(3),'20',620);


%% Martens function
function [FF1, FF2, FF3] = Martens(rest,d,c,x,y);
%rest = resting length
%d = diameter
%c = coefficient matrix c0 c1 c2 c3 d0
%x = contraction
%y = pressure normalized by 620 kPa
c0 = c(1); 
c1 = c(2);
c2 = c(3);
c3 = c(4);
d0 = c(5);
t0 = 28.6;  %initial fiber angle
tcorr = t0+d0;  %corrected fiber angle
D0 = d; %uninflated diameter
H0 = 0.0018; %Initial membrane thickness
p = 620*y; %Pressure redimensionalized

Lfiber = rest./cosd(tcorr);     %Fiber length, constant
n = rest.*tand(tcorr)./(pi*D0); %Not explained. Number of fibers?
D = ((Lfiber^2-L.^2)^(1/2))/(n*pi);     %Diameter as a function of length

derV = Lfiber^2/(4*pi*n.^2)-L.^2/(2*pi*n^2); %Derivative dV/dL
derD = -(L.*((Lfiber^2-L.^2)^(-1/2)))/(n*pi); %Derivative dD/dL

Eru = c3*L.^3+c2*L.^2+c1*L+c0;  %Modulus of elasticity

sigL = Eru.*(L-L0)/L0;      %Tension in length direction
sigPE = Eru.*(D-D0)/D0;     %Tension in perimeter direction

FL = sigL*H0*pi*D;      %Force to deform membrane in length direction
FPE = sigPE*H0*L*pi;    %Force to deform membrane perimeter

FF = -p.*derV + FPE.*derD - FL;
end

%% Sarosi function
function [FF1, FF2, FF3] = sar(rest,dia,u,x,y)
%rest = resting length;
%dia = diameter
%u = coefficients
%relative contraction
%y = pressure
a = u(1);
b = u(2);
c = u(3);
d = u(4);
e = u(5);
f = u(6);

kmax10 = 0.2; %as seen graphically in https://publicatio.bibl.u-szeged.hu/22381/1/scientific_bulletin_2012.pdf
kmax20 = []; %Solve using fzero?
p = y*620;  %kPa


FF = (a*p+b)*exp(c*k)+d*p*k+e*p+f;

end