clc                 %Command Line Clear
clear               %Clear the workspace of stored variables
close all           %Close all open figures

v1 = [0.1, 0.2, 0];

v2 = [0.4 0.5 0];

v3 = [0.6, -0.7, 0];

T = [1, 0, 0, 1;
     0, 1, 0, 0.5;
     0, 0, 1, 0;
     0, 0, 0, 1];
 
% v = T*x
v = RowVecTrans(T, v2);

mL1 = norm(v1 - v);

mL2 = norm(v2 - v3);
 
mL = mL1 + mL2;

direction = v1 - v;

unitDirection = direction/norm(direction);

mA = v2 - dot(unitDirection, v2)*unitDirection;


theta = 90*pi/180;
nT = [cos(theta) sin(theta) 0 1;
     -sin(theta) cos(theta) 0 0.5;
     0 0 1 0;
     0 0 0 1];
 
 nv = RowVecTrans(nT, v2)`