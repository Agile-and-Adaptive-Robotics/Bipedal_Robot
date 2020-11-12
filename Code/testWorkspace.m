clc                 %Command Line Clear
clear               %Clear the workspace of stored variables
close all           %Close all open figures

addpath Functions

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


%----------- fletcher method
p = RowVecTrans(T\eye(4), v1) - v2;
z = [0 0 1];

mAf = dot(v2,(cross(p, z)/norm(cross(p, z))))

%--------------
theta = 90*pi/180;
nT = [cos(theta) sin(theta) 0 1;
     -sin(theta) cos(theta) 0 0.5;
     0 0 1 0;
     0 0 0 1];
 
 v1p = RowVecTrans(nT\eye(4), v1)
 v2p = RowVecTrans(nT, v2)
 v3p = RowVecTrans(nT, v3)
 
 mAn = v2 - (v1p - v2)/norm(v1p - v2)*(dot((v1p - v2)/norm(v1p - v2), v2))