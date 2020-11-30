clc                 %Command Line Clear
clear               %Clear the workspace of stored variables
close all           %Close all open figures

addpath C:\Users\Connor\Documents\GitHub\Bipedal_Robot\Code\Functions

% v1 = [0.1, 0.2, 0];
% 
% v2 = [0.4 0.5 0];
% 
% v3 = [0.6, -0.7, 0];
% 
% T = [1, 0, 0, 1;
%      0, 1, 0, 0.5;
%      0, 0, 1, 0;
%      0, 0, 0, 1];
%  
% % v = T*x
% v = RowVecTrans(T, v2);
% 
% mL1 = norm(v1 - v);
% 
% mL2 = norm(v2 - v3);
%  
% mL = mL1 + mL2;
% 
% direction = v1 - v;
% 
% unitDirection = direction/norm(direction);
% 
% mA = v2 - dot(unitDirection, v2)*unitDirection;
% 
% 
% %----------- fletcher method
% p = RowVecTrans(T\eye(4), v1) - v2;
% z = [0 0 1];
% 
% mAf = dot(v2,(cross(p, z)/norm(cross(p, z))))
% 
% %--------------
% theta = 90*pi/180;
% nT = [cos(theta) sin(theta) 0 1;
%      -sin(theta) cos(theta) 0 0.5;
%      0 0 1 0;
%      0 0 0 1];
%  
%  v1p = RowVecTrans(nT\eye(4), v1);
%  v2p = RowVecTrans(nT, v2);
%  v3p = RowVecTrans(nT, v3);
%  
%  mAn = v2 - (v1p - v2)/norm(v1p - v2)*(dot((v1p - v2)/norm(v1p - v2), v2));
 
 % Testing biarticulate
 v1 = [0.8, 0.7, 0];
 v2 = [0.5, 0.4, 0];
 
 T1 = [1 0 0 1;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
   
 T2 = [1 0 0 1;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
   
   v2p = RowVecTrans(T2*T1, v2)
   
   theta = -90*pi/180;
   nT2 = [cos(theta) -sin(theta) 0 1;
     sin(theta) cos(theta) 0 0;
     0 0 1 0;
     0 0 0 1];
 
 v2pp = RowVecTrans(T1*nT2, v2)

 v2ppp = RowVecTrans(nT2*nT2, v2)
 
 nT1 = nT2;
 
 v2pppp = RowVecTrans(nT1*T2, v2)