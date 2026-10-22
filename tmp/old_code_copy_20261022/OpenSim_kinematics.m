%%Human Leg Kinematics
%Section 1 Loading prep files, Range of motion (RoM) setup. Choose which
%degrees of freedom (DoFs) to examinate
%Section 2: for loops break DoFs (can also be 1 DoF) into 100 data points
%and calculates T matrix, length, moment arm, force, and torque
%Section 3 is for plots and surface figures
%Section 4 is to run size calculations
clear; clc; close all;

%Units in m, N, rad/s unless otherwise stated.

load AttachPoints.mat  %Attachment points for muscles with stationary attachment points
load MovePoints.mat    %Moving attachment points and knee data
load MuscProp.mat      %Muscle specific properties
%IndexPrep %Matrices for length, force, and moment arm  
k = 180/pi;  %convert radians to degrees
c = pi/180;  %convert degrees to radians

%Range of motion, radians
range_hz = c*[ -20, 140];    %Hip (Wikipedia)
thetamin_hz = range_hz(1);      %20 deg extention
thetamax_hz = range_hz(2);      %140 deg flexion
range_k = c*[-120, 10];    %Knee
thetamin_k = range_k(1);       %120 deg flexion
thetamax_k = range_k(2);       %10 deg extension
thetamin_a = -50*c;       %ankle
thetamax_a = 20*c;
thetamin_s = -35*c;       %subtalar
thetamax_s = 25*c;
thetamin_m = -30*c;       %mtp
thetamax_m = 80*c;
thetamin_bx = -20*c;       %Lumbar bending  (I.A. Kapandji, Phsiology of the Joints. Churchill Livingstone, 1986)
thetamax_bx = 20*c;
thetamin_by = -5*c;       %Lumbar rotation  (Kapandji 1986)
thetamax_by = 5*c;
thetamin_bz = -60*c;       %Back flexion-/extension+  (Kapandji 1986)
thetamax_bz = 20*c;

%Indexing prep, see IndexPrep.m
s = 100;
t = 1; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)
%%
%%T matrix, length, moment arm, force, and torque calculations
for j=1:t
 for i = 1:s
%Degree of freedom, position or range
theta_h_z = 0; %thetamin_hz + (thetamax_hz-thetamin_hz)*(j-1)/s; %Hip, right, extension-/flexion+
thetamin_hy = c*(-32.5 - (theta_h_z/(range_hz(2)-range_hz(1)))*20); %Max External Rotation (linear interpolation of Wikipedia RoM)
thetamax_hy = c*(40);  %Max Internal Rotation (Wikipedia)
theta_h_y = 0; %Hip, right, external-/internal+ rotation
thetamin_hx = c*(-53.75 - (theta_h_z/(range_hz(2)-range_hz(1)))*30); %Max Abduction (linear interpolation of Wikipedia RoM)
thetamax_hx = c*(28.75 - (theta_h_z/(range_hz(2)-range_hz(1)))*10); %Max Adduction (linear interpolation of Wikipedia RoM)
theta_h_x = 0; %thetamin_hx + (thetamax_hx-thetamin_hx)*(i-1)/s; %Hip, right, adduction-/abduction+
theta_k = thetamin_k + (thetamax_k-thetamin_k)*(i-1)/s; %Knee, right, flexion-/extension+ 
theta_b_x = thetamin_bx + (thetamax_bx-thetamin_bx)*(i-1)/s; %Lumbar bending left-/right+
theta_b_y = 0; %thetamin_by + (thetamax_by-thetamin_by)*(j-1)/s; %Lumbar rotation cw-/ccw+
theta_b_z = thetamin_bz + (thetamax_bz-thetamin_bz)*(j-1)/s; %Back flexion-/extension+
theta_a = 0; %thetamin_a + (thetamax_a-thetamin_a)*(j-1)/s; %Ankle, right, plantar flexion-/dorsiflexion+
theta_s = 0; %thetamin_s + (thetamax_s-thetamin_s)*(i-1)/s; %Subtalar, right, inversion-/eversion+ rotation
theta_m = 0; %thetamin_m + (thetamax_m-thetamin_m)*(i-1)/s; %Metatarsophalangeal (mtp), right, flexion-/extension+

%Index for plotting (for 2 dof, j index to x axis, i index to y axis)
%theta_hZ(j,1) = theta_h_z*k; %Convert from radians to degrees
%theta_hY(i,1) = theta_h_y*k;
%theta_hX(i,1) = theta_h_x*k;
%theta_K(i,1) = theta_k*k;
theta_bX(i,1) = theta_b_x*k;
%theta_bY(j,1) = theta_b_y*k;
theta_bZ(j,1) = theta_b_z*k;
%theta_A(j,1) = theta_a*k;
%theta_S(i,1) = theta_s*k;
%theta_M(i,1) = theta_m*k;

%Rotation matrices
%Hip
R_h_x = [1, 0, 0; 0, cos(theta_h_x), -sin(theta_h_x); 0, sin(theta_h_x), cos(theta_h_x)];
R_h_y = [cos(theta_h_y), 0, sin(theta_h_y); 0, 1, 0; -sin(theta_h_y), 0, cos(theta_h_y)];
R_h_z = [cos(theta_h_z), -sin(theta_h_z), 0; sin(theta_h_z), cos(theta_h_z), 0; 0, 0, 1]; 
R_h = R_h_x*R_h_y*R_h_z;
%Knee
R_k = [cos(theta_k), -sin(theta_k), 0; sin(theta_k), cos(theta_k), 0; 0, 0, 1];
%Back
R_b_x = [1, 0, 0; 0, cos(theta_b_x), -sin(theta_b_x); 0, sin(theta_b_x), cos(theta_b_x)];
R_b_y = [cos(theta_b_y), 0, sin(theta_b_y); 0, 1, 0; -sin(theta_b_y), 0, cos(theta_b_y)];
R_b_z = [cos(theta_b_z), -sin(theta_b_z), 0; sin(theta_b_z), cos(theta_b_z), 0; 0, 0, 1];
R_b = R_b_x*R_b_y*R_b_z;
%Ankle
omg_a = [-0.10501355; -0.17402245; 0.97912632]; %Axis of rotation
R_a = [cos(theta_a)+(omg_a(1))^2*(1-cos(theta_a)), omg_a(1)*omg_a(2)*(1-cos(theta_a))-omg_a(3)*sin(theta_a), omg_a(1)*omg_a(3)*(1-cos(theta_a))+omg_a(2)*sin(theta_a); ...
         omg_a(1)*omg_a(2)*(1-cos(theta_a))+omg_a(3)*sin(theta_a), cos(theta_a)+(omg_a(2))^2*(1-cos(theta_a)), omg_a(2)*omg_a(3)*(1-cos(theta_a))-omg_a(1)*sin(theta_a); ...
         omg_a(1)*omg_a(3)*(1-cos(theta_a))-omg_a(2)*sin(theta_a), omg_a(2)*omg_a(3)*(1-cos(theta_a))+omg_a(1)*sin(theta_a), cos(theta_a)+(omg_a(3))^2*(1-cos(theta_a))]; %Ref. Modern Robotics page 74
%Subtalar
omg_s = [0.7871796; 0.60474746; -0.12094949]; %Axis of rotation
R_s = [cos(theta_s)+(omg_s(1))^2*(1-cos(theta_s)), omg_s(1)*omg_s(2)*(1-cos(theta_s))-omg_s(3)*sin(theta_s), omg_s(1)*omg_s(3)*(1-cos(theta_s))+omg_s(2)*sin(theta_s); ...
         omg_s(1)*omg_s(2)*(1-cos(theta_s))+omg_s(3)*sin(theta_s), cos(theta_s)+(omg_s(2))^2*(1-cos(theta_s)), omg_s(2)*omg_s(3)*(1-cos(theta_s))-omg_s(1)*sin(theta_s); ...
         omg_s(1)*omg_s(3)*(1-cos(theta_s))-omg_s(2)*sin(theta_s), omg_s(2)*omg_s(3)*(1-cos(theta_s))+omg_s(1)*sin(theta_s), cos(theta_s)+(omg_s(3))^2*(1-cos(theta_s))];
%Mtp
omg_m = [-0.5809544; 0; 0.81393611]; %Axis of rotation
R_m = [cos(theta_m)+(omg_m(1))^2*(1-cos(theta_m)), omg_m(1)*omg_m(2)*(1-cos(theta_m))-omg_m(3)*sin(theta_m), omg_m(1)*omg_m(3)*(1-cos(theta_m))+omg_m(2)*sin(theta_m); ...
         omg_m(1)*omg_m(2)*(1-cos(theta_m))+omg_m(3)*sin(theta_m), cos(theta_m)+(omg_m(2))^2*(1-cos(theta_m)), omg_m(2)*omg_m(3)*(1-cos(theta_m))-omg_m(1)*sin(theta_m); ...
         omg_m(1)*omg_m(3)*(1-cos(theta_m))-omg_m(2)*sin(theta_m), omg_m(2)*omg_m(3)*(1-cos(theta_m))+omg_m(1)*sin(theta_m), cos(theta_m)+(omg_m(3))^2*(1-cos(theta_m))];

%Joint location in home position
p_h = [-0.0707; -0.0661; 0.0835]; %Hip, right
p_k = [fcn1(theta_k); fcn2(theta_k); 0]; %Knee, function of knee angle (see MovePoints)
p_b = [-0.1007; 0.0815; 0]; %Back
p_a = [0; -0.43; 0]; %Ankle, right
p_s = [-0.04877; -0.04195; 0.00792]; %Subtalar, right
p_m = [0.1788; -0.002; 0.00108]; %mtp, right

%Transformation matrices
T_p = [1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1]; %Pelvis
T_h = RpToTrans(R_h,p_h); %Hip to Pelvis, right
T_k = RpToTrans(R_k,p_k); %Knee to Hip, right
T_b = RpToTrans(R_b,p_b); %Back to Pelvis
T_a = RpToTrans(R_a,p_a); %Ankle to Knee, right
T_s = RpToTrans(R_s,p_s); %Subtalar to Ankle, right
T_m = RpToTrans(R_m,p_m); %Mtp to Subtalar, right

%Muscle length, moment arm, force, and torque
%Gluteus Medius, 1
% glut_med1_length(i,j) = norm(VecTrans(T_h,glut_med1_i) - glut_med1_o); %Total length
% glut_med1_maX(i,j) = CrossProd(glut_med1_i, (VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i)/norm(VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i),1); %Moment arm, x axis
% glut_med1_maY(i,j) = CrossProd(glut_med1_i, (VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i)/norm(VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i),2); %Moment arm, y axis
% glut_med1_maZ(i,j) = CrossProd(glut_med1_i, (VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i)/norm(VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i),3); %Moment arm, z axis
% glut_med1_f(i,j) = forz(glut_med1_length(i,j), glut_med1_mif, glut_med1_ofl, glut_med1_tsl, glut_med1_pa); %Force
% glut_med1_tqX(i,j) = glut_med1_f(i,j)*glut_med1_maX(i,j); %Torque, x axis
% glut_med1_tqY(i,j) = glut_med1_f(i,j)*glut_med1_maY(i,j); %Torque, y axis
% glut_med1_tqZ(i,j) = glut_med1_f(i,j)*glut_med1_maZ(i,j); %Torque, z axis
%Gluteus Medius, 2
% glut_med2_length(i,j) = norm(VecTrans(T_h,glut_med2_i) - glut_med2_o); %Total length
% glut_med2_maX(i,j) = CrossProd(glut_med2_i, (VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i)/norm(VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i),1); %Moment arm, x axis
% glut_med2_maY(i,j) = CrossProd(glut_med2_i, (VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i)/norm(VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i),2); %Moment arm, y axis
% glut_med2_maZ(i,j) = CrossProd(glut_med2_i, (VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i)/norm(VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i),3); %Moment arm, z axis
% glut_med2_f(i,j) = forz(glut_med2_length(i,j), glut_med2_mif, glut_med2_ofl, glut_med2_tsl, glut_med2_pa); %Force
% glut_med2_tqX(i,j) = glut_med2_f(i,j)*glut_med2_maX(i,j); %Torque, x axis
% glut_med2_tqY(i,j) = glut_med2_f(i,j)*glut_med2_maY(i,j); %Torque, y axis
% glut_med2_tqZ(i,j) = glut_med2_f(i,j)*glut_med2_maZ(i,j); %Torque, z axis
%Gluteus Medius, 3
% glut_med3_length(i,j) = norm(VecTrans(T_h,glut_med3_i) - glut_med3_o); %Total length
% glut_med3_maX(i,j) = CrossProd(glut_med3_i, (VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i)/norm(VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i),1); %Moment arm, x axis
% glut_med3_maY(i,j) = CrossProd(glut_med3_i, (VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i)/norm(VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i),2); %Moment arm, y axis
% glut_med3_maZ(i,j) = CrossProd(glut_med3_i, (VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i)/norm(VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i),3); %Moment arm, z axis
% glut_med3_f(i,j) = forz(glut_med3_length(i,j), glut_med3_mif, glut_med3_ofl, glut_med3_tsl, glut_med3_pa); %Force
% glut_med3_tqX(i,j) = glut_med3_f(i,j)*glut_med3_maX(i,j); %Torque, x axis
% glut_med3_tqY(i,j) = glut_med3_f(i,j)*glut_med3_maY(i,j); %Torque, y axis
% glut_med3_tqZ(i,j) = glut_med3_f(i,j)*glut_med3_maZ(i,j); %Torque, z axis
%Gluteus Minimus, 1
% glut_min1_length(i,j) = norm(VecTrans(T_h,glut_min1_i) - glut_min1_o); %Total length
% glut_min1_maX(i,j) = CrossProd(glut_min1_i, (VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i)/norm(VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i),1); %Moment arm, x axis
% glut_min1_maY(i,j) = CrossProd(glut_min1_i, (VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i)/norm(VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i),2); %Moment arm, y axis
% glut_min1_maZ(i,j) = CrossProd(glut_min1_i, (VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i)/norm(VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i),3); %Moment arm, z axis
% glut_min1_f(i,j) = forz(glut_min1_length(i,j), glut_min1_mif, glut_min1_ofl, glut_min1_tsl, glut_min1_pa); %Force
% glut_min1_tqX(i,j) = glut_min1_f(i,j)*glut_min1_maX(i,j); %Torque, x axis
% glut_min1_tqY(i,j) = glut_min1_f(i,j)*glut_min1_maY(i,j); %Torque, y axis
% glut_min1_tqZ(i,j) = glut_min1_f(i,j)*glut_min1_maZ(i,j); %Torque, z axis
%Gluteus Minimus, 2
% glut_min2_length(i,j) = norm(VecTrans(T_h,glut_min2_i) - glut_min2_o); %Total length
% glut_min2_maX(i,j) = CrossProd(glut_min2_i, (VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i)/norm(VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i),1); %Moment arm, x axis
% glut_min2_maY(i,j) = CrossProd(glut_min2_i, (VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i)/norm(VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i),2); %Moment arm, y axis
% glut_min2_maZ(i,j) = CrossProd(glut_min2_i, (VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i)/norm(VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i),3); %Moment arm, z axis
% glut_min2_f(i,j) = forz(glut_min2_length(i,j), glut_min2_mif, glut_min2_ofl, glut_min2_tsl, glut_min2_pa); %Force
% glut_min2_tqX(i,j) = glut_min2_f(i,j)*glut_min2_maX(i,j); %Torque, x axis
% glut_min2_tqY(i,j) = glut_min2_f(i,j)*glut_min2_maY(i,j); %Torque, y axis
% glut_min2_tqZ(i,j) = glut_min2_f(i,j)*glut_min2_maZ(i,j); %Torque, z axis
%Gluteus Minimus, 3
% glut_min3_length(i,j) = norm(VecTrans(T_h,glut_min3_i) - glut_min3_o); %Total length
% glut_min3_maX(i,j) = CrossProd(glut_min3_i, (VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i)/norm(VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i),1); %Moment arm, x axis
% glut_min3_maY(i,j) = CrossProd(glut_min3_i, (VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i)/norm(VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i),2); %Moment arm, y axis
% glut_min3_maZ(i,j) = CrossProd(glut_min3_i, (VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i)/norm(VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i),3); %Moment arm, z axis
% glut_min3_f(i,j) = forz(glut_min3_length(i,j), glut_min3_mif, glut_min3_ofl, glut_min3_tsl, glut_min3_pa); %Force
% glut_min3_tqX(i,j) = glut_min3_f(i,j)*glut_min3_maX(i,j); %Torque, x axis
% glut_min3_tqY(i,j) = glut_min3_f(i,j)*glut_min3_maY(i,j); %Torque, y axis
% glut_min3_tqZ(i,j) = glut_min3_f(i,j)*glut_min3_maZ(i,j); %Torque, z axis
%Gluteus Maximus, 1
% glut_max1_length(i,j) = norm(glut_max1_wr1-glut_max1_o)+norm(VecTrans(T_h,glut_max1_wr2) - glut_max1_wr1)+norm(glut_max1_i-glut_max1_wr2); %Total length
% glut_max1_maX(i,j) = CrossProd(glut_max1_wr2, (VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2)/norm(VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2),1); %Moment arm, x axis
% glut_max1_maY(i,j) = CrossProd(glut_max1_wr2, (VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2)/norm(VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2),2); %Moment arm, y axis
% glut_max1_maZ(i,j) = CrossProd(glut_max1_wr2, (VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2)/norm(VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2),3); %Moment arm, z axis
% glut_max1_f(i,j) = forz(glut_max1_length(i,j), glut_max1_mif, glut_max1_ofl, glut_max1_tsl, glut_max1_pa); %Force
% glut_max1_tqX(i,j) = glut_max1_f(i,j)*glut_max1_maX(i,j); %Torque, x axis
% glut_max1_tqY(i,j) = glut_max1_f(i,j)*glut_max1_maY(i,j); %Torque, y axis
% glut_max1_tqZ(i,j) = glut_max1_f(i,j)*glut_max1_maZ(i,j); %Torque, z axis
%Gluteus Maximus, 2
% glut_max2_length(i,j) = norm(glut_max2_wr1-glut_max2_o)+norm(VecTrans(T_h,glut_max2_wr2) - glut_max2_wr1)+norm(glut_max2_i-glut_max2_wr2); %Total length
% glut_max2_maX(i,j) = CrossProd(glut_max2_wr2, (VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2)/norm(VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2),1); %Moment arm, x axis
% glut_max2_maY(i,j) = CrossProd(glut_max2_wr2, (VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2)/norm(VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2),2); %Moment arm, y axis
% glut_max2_maZ(i,j) = CrossProd(glut_max2_wr2, (VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2)/norm(VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2),3); %Moment arm, z axis
% glut_max2_f(i,j) = forz(glut_max2_length(i,j), glut_max2_mif, glut_max2_ofl, glut_max2_tsl, glut_max2_pa); %Force
% glut_max2_tqX(i,j) = glut_max2_f(i,j)*glut_max2_maX(i,j); %Torque, x axis
% glut_max2_tqY(i,j) = glut_max2_f(i,j)*glut_max2_maY(i,j); %Torque, y axis
% glut_max2_tqZ(i,j) = glut_max2_f(i,j)*glut_max2_maZ(i,j); %Torque, z axis
%Gluteus Maximus, 3
% glut_max3_length(i,j) = norm(glut_max3_wr1-glut_max3_o)+norm(VecTrans(T_h,glut_max3_wr2) - glut_max3_wr1)+norm(glut_max3_i-glut_max3_wr2); %Total length
% glut_max3_maX(i,j) = CrossProd(glut_max3_wr2, (VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2)/norm(VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2),1); %Moment arm, x axis
% glut_max3_maY(i,j) = CrossProd(glut_max3_wr2, (VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2)/norm(VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2),2); %Moment arm, y axis
% glut_max3_maZ(i,j) = CrossProd(glut_max3_wr2, (VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2)/norm(VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2),3); %Moment arm, z axis
% glut_max3_f(i,j) = forz(glut_max3_length(i,j), glut_max3_mif, glut_max3_ofl, glut_max3_tsl, glut_max3_pa); %Force
% glut_max3_tqX(i,j) = glut_max3_f(i,j)*glut_max3_maX(i,j); %Torque, x axis
% glut_max3_tqY(i,j) = glut_max3_f(i,j)*glut_max3_maY(i,j); %Torque, y axis
% glut_max3_tqZ(i,j) = glut_max3_f(i,j)*glut_max3_maZ(i,j); %Torque, z axis
%Semimembranosus
% semimem_length(i,j) = norm(VecTrans(T_h*T_k,semimem_i)-semimem_o); %Total length
% semimem_maHX(i,j) = CrossProd(VecTrans(T_k,semimem_i), (VecTrans(T_h\eye(4),semimem_o)-VecTrans(T_k,semimem_i))/norm(VecTrans(T_h\eye(4),semimem_o)-VecTrans(T_k,semimem_i)),1); %Moment arm, Hip, x axis
% semimem_maHY(i,j) = CrossProd(VecTrans(T_k,semimem_i), (VecTrans(T_h\eye(4),semimem_o)-VecTrans(T_k,semimem_i))/norm(VecTrans(T_h\eye(4),semimem_o)-VecTrans(T_k,semimem_i)),2); %Moment arm, Hip, y axis
% semimem_maHZ(i,j) = CrossProd(VecTrans(T_k,semimem_i), (VecTrans(T_h\eye(4),semimem_o)-VecTrans(T_k,semimem_i))/norm(VecTrans(T_h\eye(4),semimem_o)-VecTrans(T_k,semimem_i)),3); %Moment arm, Hip, z axis
% semimem_maKZ(i,j) = CrossProd(semimem_i, (VecTrans((T_h*T_k)\eye(4),semimem_o)-semimem_i)/norm(VecTrans((T_h*T_k)\eye(4),semimem_o)-semimem_i),3); %Moment arm, Knee, z axis
% semimem_f(i,j) = forz(semimem_length(i,j), semimem_mif, semimem_ofl, semimem_tsl, semimem_pa); %Force
% semimem_tqHX(i,j) = semimem_f(i,j)*semimem_maHX(i,j); %Torque, Hip, x axis
% semimem_tqHY(i,j) = semimem_f(i,j)*semimem_maHY(i,j); %Torque, Hip, y axis
% semimem_tqHZ(i,j) = semimem_f(i,j)*semimem_maHZ(i,j); %Torque, Hip, z axis
% semimem_tqKZ(i,j) = semimem_f(i,j)*semimem_maKZ(i,j); %Torque, Knee, z axis
%Semitendinosus
% semiten_length(i,j) = norm(VecTrans(T_h*T_k,semiten_wr2)-semiten_o)+norm(semiten_wr3-semiten_wr2)+norm(semiten_wr3-semiten_i); %Total length
% semiten_maHX(i,j) = CrossProd(VecTrans(T_k,semiten_wr2), (VecTrans(T_h\eye(4),semiten_o)-VecTrans(T_k,semiten_wr2))/norm(VecTrans(T_h\eye(4),semiten_o)-VecTrans(T_k,semiten_wr2)),1); %Moment arm, Hip, x axis
% semiten_maHY(i,j) = CrossProd(VecTrans(T_k,semiten_wr2), (VecTrans(T_h\eye(4),semiten_o)-VecTrans(T_k,semiten_wr2))/norm(VecTrans(T_h\eye(4),semiten_o)-VecTrans(T_k,semiten_wr2)),2); %Moment arm, Hip, y axis
% semiten_maHZ(i,j) = CrossProd(VecTrans(T_k,semiten_wr2), (VecTrans(T_h\eye(4),semiten_o)-VecTrans(T_k,semiten_wr2))/norm(VecTrans(T_h\eye(4),semiten_o)-VecTrans(T_k,semiten_wr2)),3); %Moment arm, Hip, z axis
% semiten_maKZ(i,j) = CrossProd(semiten_wr2, (VecTrans((T_h*T_k)\eye(4),semiten_o)-semiten_wr2)/norm(VecTrans((T_h*T_k)\eye(4),semiten_o)-semiten_wr2),3); %Moment arm, Knee, z axis
% semiten_f(i,j) = forz(semiten_length(i,j), semiten_mif, semiten_ofl, semiten_tsl, semiten_pa); %Force
% semiten_tqHX(i,j) = semiten_f(i,j)*semiten_maHX(i,j); %Torque, Hip, x axis
% semiten_tqHY(i,j) = semiten_f(i,j)*semiten_maHY(i,j); %Torque, Hip, y axis
% semiten_tqHZ(i,j) = semiten_f(i,j)*semiten_maHZ(i,j); %Torque, Hip, z axis
% semiten_tqKZ(i,j) = semiten_f(i,j)*semiten_maKZ(i,j); %Torque, Knee, z axis
%Bicep Femoris, Long Head
% bifemlh_length(i,j) = norm(VecTrans(T_h*T_k,bifemlh_wr1)-bifemlh_o)+norm(bifemlh_i-bifemlh_wr1); %Total length
% bifemlh_maHX(i,j) = CrossProd(VecTrans(T_k,bifemlh_wr1), (VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1))/norm(VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1)),1); %Moment arm, Hip, x axis
% bifemlh_maHY(i,j) = CrossProd(VecTrans(T_k,bifemlh_wr1), (VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1))/norm(VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1)),2); %Moment arm, Hip, y axis
% bifemlh_maHZ(i,j) = CrossProd(VecTrans(T_k,bifemlh_wr1), (VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1))/norm(VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1)),3); %Moment arm, Hip, z axis
% bifemlh_maKZ(i,j) = CrossProd(bifemlh_wr1, (VecTrans((T_h*T_k)\eye(4),bifemlh_o)-bifemlh_wr1)/norm(VecTrans((T_h*T_k)\eye(4),bifemlh_o)-bifemlh_wr1),3); %Moment arm, Knee, z axis
% bifemlh_f(i,j) = forz(bifemlh_length(i,j), bifemlh_mif, bifemlh_ofl, bifemlh_tsl, bifemlh_pa); %Force
% bifemlh_tqHX(i,j) = bifemlh_f(i,j)*semiten_maHX(i,j); %Torque, Hip, x axis
% bifemlh_tqHY(i,j) = bifemlh_f(i,j)*bifemlh_maHY(i,j); %Torque, Hip, y axis
% bifemlh_tqHZ(i,j) = bifemlh_f(i,j)*bifemlh_maHZ(i,j); %Torque, Hip, z axis
% bifemlh_tqKZ(i,j) = bifemlh_f(i,j)*bifemlh_maKZ(i,j); %Torque, Knee, z axis
%Bicep Femoris, Short Head
% bifemsh_length(i,j) = norm(VecTrans(T_k,bifemsh_wr1)-bifemsh_o)+norm(bifemsh_i-bifemsh_wr1); %Total length, right
% bifemsh_maZ(i,j) = CrossProd(bifemsh_wr1, (VecTrans(T_k\eye(4),bifemsh_o)-bifemsh_wr1)/norm(VecTrans(T_k\eye(4),bifemsh_o)-bifemsh_wr1),3); %Moment arm, Knee, z axis
% bifemsh_f(i,j) = forz(bifemsh_length(i,j), bifemsh_mif, bifemsh_ofl, bifemsh_tsl, bifemsh_pa); %Force
% bifemsh_tqZ(i,j) = bifemsh_f(i,j)*bifemsh_maZ(i,j); %Torque, Knee, z axis
%Sartorius
% sar_length(i,j) = norm(VecTrans(T_h,sar_wr1)-sar_o)+norm(VecTrans(T_k,sar_wr2)-sar_wr1)+norm(sar_wr3-sar_wr2)+norm(sar_wr3-sar_i); %Total length, right
% sar_maHX(i,j) = CrossProd(sar_wr1, (VecTrans(T_h\eye(4),sar_o)-sar_wr1)/norm(VecTrans(T_h\eye(4),sar_o)-sar_wr1),1); %Moment arm, Hip, x axis
% sar_maHY(i,j) = CrossProd(sar_wr1, (VecTrans(T_h\eye(4),sar_o)-sar_wr1)/norm(VecTrans(T_h\eye(4),sar_o)-sar_wr1),2); %Moment arm, Hip, y axis
% sar_maHZ(i,j) = CrossProd(sar_wr1, (VecTrans(T_h\eye(4),sar_o)-sar_wr1)/norm(VecTrans(T_h\eye(4),sar_o)-sar_wr1),3); %Moment arm, Hip, z axis
% sar_maKZ(i,j) = CrossProd(sar_wr2, (VecTrans(T_k\eye(4),sar_wr1)-sar_wr2)/norm(VecTrans(T_k\eye(4),sar_wr1)-sar_wr2),3); %Moment arm, Knee, z axis
% sar_f(i,j) = forz(sar_length(i,j), sar_mif, sar_ofl, sar_tsl, sar_pa); %Force
% sar_tqHX(i,j) = sar_f(i,j)*sar_maHX(i,j); %Torque, Hip, x axis
% sar_tqHY(i,j) = sar_f(i,j)*sar_maHY(i,j); %Torque, Hip, y axis
% sar_tqHZ(i,j) = sar_f(i,j)*sar_maHZ(i,j); %Torque, Hip, z axis
% sar_tqKZ(i,j) = sar_f(i,j)*sar_maKZ(i,j); %Torque, Knee, z axis
%Adductor Longus
% add_long_length(i,j) = norm(VecTrans(T_h,add_long_i) - add_long_o); %Total length, right
% add_long_maX(i,j) = CrossProd(add_long_i, (VecTrans(T_h\eye(4),add_long_o)-add_long_i)/norm(VecTrans(T_h\eye(4),add_long_o)-add_long_i),1); %Moment arm, x axis
% add_long_maY(i,j) = CrossProd(add_long_i, (VecTrans(T_h\eye(4),add_long_o)-add_long_i)/norm(VecTrans(T_h\eye(4),add_long_o)-add_long_i),2); %Moment arm, y axis
% add_long_maZ(i,j) = CrossProd(add_long_i, (VecTrans(T_h\eye(4),add_long_o)-add_long_i)/norm(VecTrans(T_h\eye(4),add_long_o)-add_long_i),3); %Moment arm, z axis
% add_long_f(i,j) = forz(add_long_length(i,j), add_long_mif, add_long_ofl, add_long_tsl, add_long_pa); %Force
% add_long_tqX(i,j) = add_long_f(i,j)*add_long_maX(i,j); %Torque, x axis
% add_long_tqY(i,j) = add_long_f(i,j)*add_long_maY(i,j); %Torque, y axis
% add_long_tqZ(i,j) = add_long_f(i,j)*add_long_maZ(i,j); %Torque, z axis
%Adductor Brevis
% add_brev_length(i,j) = norm(VecTrans(T_h,add_brev_i) - add_brev_o); %Total length
% add_brev_maX(i,j) = CrossProd(add_brev_i, (VecTrans(T_h\eye(4),add_brev_o)-add_brev_i)/norm(VecTrans(T_h\eye(4),add_brev_o)-add_brev_i),1); %Moment arm, x axis
% add_brev_maY(i,j) = CrossProd(add_brev_i, (VecTrans(T_h\eye(4),add_brev_o)-add_brev_i)/norm(VecTrans(T_h\eye(4),add_brev_o)-add_brev_i),2); %Moment arm, y axis
% add_brev_maZ(i,j) = CrossProd(add_brev_i, (VecTrans(T_h\eye(4),add_brev_o)-add_brev_i)/norm(VecTrans(T_h\eye(4),add_brev_o)-add_brev_i),3); %Moment arm, z axis
% add_brev_f(i,j) = forz(add_brev_length(i,j), add_brev_mif, add_brev_ofl, add_brev_tsl, add_brev_pa); %Force
% add_brev_tqX(i,j) = add_brev_f(i,j)*add_brev_maX(i,j); %Torque, x axis
% add_brev_tqY(i,j) = add_brev_f(i,j)*add_brev_maY(i,j); %Torque, y axis
% add_brev_tqZ(i,j) = add_brev_f(i,j)*add_brev_maZ(i,j); %Torque, z axis
%Adductor Magnus, 1
% add_mag1_length(i,j) = norm(VecTrans(T_h,add_mag1_i) - add_mag1_o); %Total length
% add_mag1_maX(i,j) = CrossProd(add_mag1_i, (VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i)/norm(VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i),1); %Moment arm, x axis
% add_mag1_maY(i,j) = CrossProd(add_mag1_i, (VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i)/norm(VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i),2); %Moment arm, y axis
% add_mag1_maZ(i,j) = CrossProd(add_mag1_i, (VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i)/norm(VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i),3); %Moment arm, z axis
% add_mag1_f(i,j) = forz(add_mag1_length(i,j), add_mag1_mif, add_mag1_ofl, add_mag1_tsl, add_mag1_pa); %Force
% add_mag1_tqX(i,j) = add_mag1_f(i,j)*add_mag1_maX(i,j); %Torque, x axis
% add_mag1_tqY(i,j) = add_mag1_f(i,j)*add_mag1_maY(i,j); %Torque, y axis
% add_mag1_tqZ(i,j) = add_mag1_f(i,j)*add_mag1_maZ(i,j); %Torque, z axis
%Adductor Magnus, 2
% add_mag2_length(i,j) = norm(VecTrans(T_h,add_mag2_i) - add_mag2_o); %Total length
% add_mag2_maX(i,j) = CrossProd(add_mag2_i, (VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i)/norm(VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i),1); %Moment arm, x axis
% add_mag2_maY(i,j) = CrossProd(add_mag2_i, (VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i)/norm(VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i),2); %Moment arm, y axis
% add_mag2_maZ(i,j) = CrossProd(add_mag2_i, (VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i)/norm(VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i),3); %Moment arm, z axis
% add_mag2_f(i,j) = forz(add_mag2_length(i,j), add_mag2_mif, add_mag2_ofl, add_mag2_tsl, add_mag2_pa); %Force
% add_mag2_tqX(i,j) = add_mag2_f(i,j)*add_mag2_maX(i,j); %Torque, x axis
% add_mag2_tqY(i,j) = add_mag2_f(i,j)*add_mag2_maY(i,j); %Torque, y axis
% add_mag2_tqZ(i,j) = add_mag2_f(i,j)*add_mag2_maZ(i,j); %Torque, z axis
%Adductor Magnus, 3
% add_mag3_length(i,j) = norm(VecTrans(T_h,add_mag3_i) - add_mag3_o); %Total length
% add_mag3_maX(i,j) = CrossProd(add_mag3_i, (VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i)/norm(VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i),1); %Moment arm, x axis
% add_mag3_maY(i,j) = CrossProd(add_mag3_i, (VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i)/norm(VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i),2); %Moment arm, y axis
% add_mag3_maZ(i,j) = CrossProd(add_mag3_i, (VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i)/norm(VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i),3); %Moment arm, z axis
% add_mag3_f(i,j) = forz(add_mag3_length(i,j), add_mag3_mif, add_mag3_ofl, add_mag3_tsl, add_mag3_pa); %Force
% add_mag3_tqX(i,j) = add_mag3_f(i,j)*add_mag3_maX(i,j); %Torque, x axis
% add_mag3_tqY(i,j) = add_mag3_f(i,j)*add_mag3_maY(i,j); %Torque, y axis
% add_mag3_tqZ(i,j) = add_mag3_f(i,j)*add_mag3_maZ(i,j); %Torque, z axis
%Tensor Fasciae Latae
% tfl_length(i,j) = norm(VecTrans(T_h,tfl_wr1)-tfl_o)+norm(tfl_wr2 - tfl_wr1)+norm(tfl_wr2-VecTrans(T_k,tfl_i)); %Total length, right
% tfl_maHX(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),tfl_o)-tfl_wr1)/norm(VecTrans(T_h\eye(4),tfl_o)-tfl_wr1),1); %Moment arm, Hip, x axis
% tfl_maHY(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),tfl_o)-tfl_wr1)/norm(VecTrans(T_h\eye(4),tfl_o)-tfl_wr1),2); %Moment arm, Hip, y axis
% tfl_maHZ(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),tfl_o)-tfl_wr1)/norm(VecTrans(T_h\eye(4),tfl_o)-tfl_wr1),3); %Moment arm, Hip, z axis
% tfl_maKZ(i,j) = CrossProd(tfl_i, (VecTrans(T_k\eye(4),tfl_wr2)-tfl_i)/norm(VecTrans(T_k\eye(4),tfl_wr2)-tfl_i),3); %Moment arm, Knee, z axis
% tfl_f(i,j) = forz(tfl_length(i,j), tfl_mif, tfl_ofl, tfl_tsl, tfl_pa); %Force
% tfl_tqHX(i,j) = tfl_f(i,j)*tfl_maHX(i,j); %Torque, Hip, x axis
% tfl_tqHY(i,j) = tfl_f(i,j)*tfl_maHY(i,j); %Torque, Hip, y axis
% tfl_tqHZ(i,j) = tfl_f(i,j)*tfl_maHZ(i,j); %Torque, Hip, z axis
% tfl_tqKZ(i,j) = tfl_f(i,j)*tfl_maKZ(i,j); %Torque, Knee, z axis
%Pectineus
% pect_length(i,j) = norm(VecTrans(T_h,pect_i) - pect_o); %Total length
% pect_maX(i,j) = CrossProd(pect_i, (VecTrans(T_h\eye(4),pect_o)-pect_i)/norm(VecTrans(T_h\eye(4),pect_o)-pect_i),1); %Moment arm, x axis
% pect_maY(i,j) = CrossProd(pect_i, (VecTrans(T_h\eye(4),pect_o)-pect_i)/norm(VecTrans(T_h\eye(4),pect_o)-pect_i),2); %Moment arm, y axis
% pect_maZ(i,j) = CrossProd(pect_i, (VecTrans(T_h\eye(4),pect_o)-pect_i)/norm(VecTrans(T_h\eye(4),pect_o)-pect_i),3); %Moment arm, z axis
% pect_f(i,j) = forz(pect_length(i,j), pect_mif, pect_ofl, pect_tsl, pect_pa); %Force
% pect_tqX(i,j) = pect_f(i,j)*pect_maX(i,j); %Torque, x axis
% pect_tqY(i,j) = pect_f(i,j)*pect_maX(i,j); %Torque, y axis
% pect_tqZ(i,j) = pect_f(i,j)*pect_maX(i,j); %Torque, z axis
%Gracilis
% grac_length(i,j) = norm(VecTrans(T_h*T_k,grac_wr2)-grac_o)+norm(grac_i-grac_wr2); %Total length, right
% grac_maHX(i,j) = CrossProd(VecTrans(T_k,grac_wr2), (VecTrans(T_h\eye(4),grac_o)-VecTrans(T_k,grac_wr2))/norm(VecTrans(T_h\eye(4),grac_o)-VecTrans(T_k,grac_wr2)),1); %Moment arm, Hip, x axis
% grac_maHY(i,j) = CrossProd(VecTrans(T_k,grac_wr2), (VecTrans(T_h\eye(4),grac_o)-VecTrans(T_k,grac_wr2))/norm(VecTrans(T_h\eye(4),grac_o)-VecTrans(T_k,grac_wr2)),2); %Moment arm, Hip, y axis
% grac_maHZ(i,j) = CrossProd(VecTrans(T_k,grac_wr2), (VecTrans(T_h\eye(4),grac_o)-VecTrans(T_k,grac_wr2))/norm(VecTrans(T_h\eye(4),grac_o)-VecTrans(T_k,grac_wr2)),3); %Moment arm, Hip, z axis
% grac_maKZ(i,j) = CrossProd(grac_wr2, (VecTrans((T_h*T_k)\eye(4),grac_o)-grac_wr2)/norm(VecTrans((T_h*T_k)\eye(4),grac_o)-grac_wr2),3); %Moment arm, Knee, z axis
% grac_f(i,j) = forz(grac_length(i,j), grac_mif, grac_ofl, grac_tsl, grac_pa); %Force
% grac_tqHX(i,j) = grac_f(i,j)*grac_maHX(i,j); %Torque, Hip, x axis
% grac_tqHY(i,j) = grac_f(i,j)*grac_maHY(i,j); %Torque, Hip, y axis
% grac_tqHZ(i,j) = grac_f(i,j)*grac_maHZ(i,j); %Torque, Hip, z axis
% grac_tqKZ(i,j) = grac_f(i,j)*grac_maKZ(i,j); %Torque, Knee, z axis
%Iliacus
% iliacus_length(i,j) = norm(iliacus_wr1-iliacus_o)+norm(VecTrans(T_h,iliacus_wr3)-iliacus_wr1)+norm(iliacus_wr3-iliacus_i); %Total length
% iliacus_maX(i,j) = CrossProd(iliacus_wr3, (VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3),1); %Moment arm, x axis
% iliacus_maY(i,j) = CrossProd(iliacus_wr3, (VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3),2); %Moment arm, y axis
% iliacus_maZ(i,j) = CrossProd(iliacus_wr3, (VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3),3); %Moment arm, z axis
% iliacus_f(i,j) = forz(iliacus_length(i,j), iliacus_mif, iliacus_ofl, iliacus_tsl, iliacus_pa); %Force
% iliacus_tqX(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, x axis
% iliacus_tqY(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, y axis
% iliacus_tqZ(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, z axis
%Psoas
% psoas_length(i,j) = norm(psoas_wr1-psoas_o)+norm(VecTrans(T_h,psoas_wr3)-psoas_wr1)+norm(psoas_wr3-psoas_i); %Total length, right
% psoas_maHX(i,j) = CrossProd(psoas_wr3, (VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3)/norm(VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3),1); %Moment arm, x axis
% psoas_maHY(i,j) = CrossProd(psoas_wr3, (VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3)/norm(VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3),2); %Moment arm, y axis
% psoas_maHZ(i,j) = CrossProd(psoas_wr3, (VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3)/norm(VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3),3); %Moment arm, z axis
% psoas_f(i,j) = forz(psoas_length(i,j), psoas_mif, psoas_ofl, psoas_tsl, psoas_pa); %Force
% psoas_tqHX(i,j) = psoas_f(i,j)*psoas_maHX(i,j); %Torque, x axis
% psoas_tqHY(i,j) = psoas_f(i,j)*psoas_maHY(i,j); %Torque, y axis
% psoas_tqHZ(i,j) = psoas_f(i,j)*psoas_maHZ(i,j); %Torque, z axis
%Quadriceps Femoris
% quad_fem_length(i,j) = norm(VecTrans(T_h,quad_fem_i) - quad_fem_o); %Total length
% quad_fem_maX(i,j) = CrossProd(quad_fem_i, (VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i)/norm(VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i),1); %Moment arm, x axis
% quad_fem_maY(i,j) = CrossProd(quad_fem_i, (VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i)/norm(VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i),2); %Moment arm, y axis
% quad_fem_maZ(i,j) = CrossProd(quad_fem_i, (VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i)/norm(VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i),3); %Moment arm, z axis
% quad_fem_f(i,j) = forz(quad_fem_length(i,j), quad_fem_mif, quad_fem_ofl, quad_fem_tsl, quad_fem_pa); %Force
% quad_fem_tqX(i,j) = quad_fem_f(i,j)*quad_fem_maX(i,j); %Torque, x axis
% quad_fem_tqY(i,j) = quad_fem_f(i,j)*quad_fem_maY(i,j); %Torque, x axis
% quad_fem_tqZ(i,j) = quad_fem_f(i,j)*quad_fem_maZ(i,j); %Torque, x axis
%Gemelli
% gem_length(i,j) = norm(VecTrans(T_h,gem_i) - gem_o); %Total length
% gem_maX(i,j) = CrossProd(gem_i, (VecTrans(T_h\eye(4),gem_o)-gem_i)/norm(VecTrans(T_h\eye(4),gem_o)-gem_i),1); %Moment arm, x axis
% gem_maY(i,j) = CrossProd(gem_i, (VecTrans(T_h\eye(4),gem_o)-gem_i)/norm(VecTrans(T_h\eye(4),gem_o)-gem_i),2); %Moment arm, y axis
% gem_maZ(i,j) = CrossProd(gem_i, (VecTrans(T_h\eye(4),gem_o)-gem_i)/norm(VecTrans(T_h\eye(4),gem_o)-gem_i),3); %Moment arm, z axis
% gem_f(i,j) = forz(gem_length(i,j), gem_mif, gem_ofl, gem_tsl, gem_pa); %Force
% gem_tqX(i,j) = gem_f(i,j)*gem_maX(i,j); %Torque, x axis
% gem_tqY(i,j) = gem_f(i,j)*gem_maX(i,j); %Torque, y axis
% gem_tqZ(i,j) = gem_f(i,j)*gem_maX(i,j); %Torque, z axis
%Piriformis
% peri_length(i,j) = norm(VecTrans(T_h,peri_i)-peri_wr1)+norm(peri_o-peri_wr1); %Total length
% peri_maX(i,j) = CrossProd(peri_i, (VecTrans(T_h\eye(4),peri_wr1)-peri_i)/norm(VecTrans(T_h\eye(4),peri_wr1)-peri_i),1); %Moment arm, x axis
% peri_maY(i,j) = CrossProd(peri_i, (VecTrans(T_h\eye(4),peri_wr1)-peri_i)/norm(VecTrans(T_h\eye(4),peri_wr1)-peri_i),2); %Moment arm, y axis
% peri_maZ(i,j) = CrossProd(peri_i, (VecTrans(T_h\eye(4),peri_wr1)-peri_i)/norm(VecTrans(T_h\eye(4),peri_wr1)-peri_i),3); %Moment arm, z axis
% peri_f(i,j) = forz(peri_length(i,j), peri_mif, peri_ofl, peri_tsl, peri_pa); %Force
% peri_tqX(i,j) = peri_f(i,j)*peri_maX(i,j); %Torque, x axis
% peri_tqY(i,j) = peri_f(i,j)*peri_maY(i,j); %Torque, y axis
% peri_tqZ(i,j) = peri_f(i,j)*peri_maZ(i,j); %Torque, z axis
%Rectus Femoris
% rect_fem_o = [-0.03; -0.031; 0.097]; %Pelvis, right, origin
% rect_fem_wr1 = [0.058; -0.339; -0.003]; %Femur
% rect_fem_i = [fcn3(theta_k); fcn4(theta_k); 0.001]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% rect_fem_length(i,j) = norm(VecTrans(T_h,rect_fem_wr1)-rect_fem_o)+norm(VecTrans(T_k,rect_fem_i)-rect_fem_wr1); %Total length
% rect_fem_maHX(i,j) = CrossProd(VecTrans(T_k,rect_fem_i), (VecTrans(T_h\eye(4),rect_fem_o)-VecTrans(T_k,rect_fem_i))/norm(VecTrans(T_h\eye(4),rect_fem_o)-VecTrans(T_k,rect_fem_i)),1); %Moment arm, Hip, x axis
% rect_fem_maHY(i,j) = CrossProd(VecTrans(T_k,rect_fem_i), (VecTrans(T_h\eye(4),rect_fem_o)-VecTrans(T_k,rect_fem_i))/norm(VecTrans(T_h\eye(4),rect_fem_o)-VecTrans(T_k,rect_fem_i)),2); %Moment arm, Hip, y axis
% rect_fem_maHZ(i,j) = CrossProd(VecTrans(T_k,rect_fem_i), (VecTrans(T_h\eye(4),rect_fem_o)-VecTrans(T_k,rect_fem_i))/norm(VecTrans(T_h\eye(4),rect_fem_o)-VecTrans(T_k,rect_fem_i)),3); %Moment arm, Hip, z axis
% rect_fem_maKZ(i,j) = CrossProd(rect_fem_i, (VecTrans((T_h*T_k)\eye(4),rect_fem_o)-rect_fem_i)/norm(VecTrans((T_h*T_k)\eye(4),rect_fem_o)-rect_fem_i),3); %Moment arm, Knee, z axis
% rect_fem_f(i,j) = forz(rect_fem_length(i,j), rect_fem_mif, rect_fem_ofl, rect_fem_tsl, rect_fem_pa); %Force
% rect_fem_tqHX(i,j) = rect_fem_f(i,j)*rect_fem_maHX(i,j); %Torque, Hip, x axis
% rect_fem_tqHY(i,j) = rect_fem_f(i,j)*rect_fem_maHY(i,j); %Torque, Hip, y axis
% rect_fem_tqHZ(i,j) = rect_fem_f(i,j)*rect_fem_maHZ(i,j); %Torque, Hip, z axis
% rect_fem_tqKZ(i,j) = rect_fem_f(i,j)*rect_fem_maKZ(i,j); %Torque, Knee, z axis
%Vastus Medialis
% vas_med_o = [0.014; -0.21; 0.019]; %Femur, right, origin
% vas_med_wr1 = [0.036; -0.277; 0.001]; %Femur, right, wrapping point 1
% vas_med_wr2 = [0.054; -0.332; -0.006]; %Femur, right, wrapping point 2
% vas_med_i = [fcn5(theta_k); fcn6(theta_k); -0.015]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% vas_med_length(i,j) = norm(vas_med_wr1-vas_med_o)+norm(vas_med_wr1-vas_med_wr2)+norm(vas_med_wr2-VecTrans(T_k,vas_med_i)); %Total length
% vas_med_maZ(i,j) = CrossProd(vas_med_i, (VecTrans(T_k\eye(4),vas_med_wr1)-vas_med_i)/norm(VecTrans(T_k\eye(4),vas_med_wr1)-vas_med_i),3); %Moment arm, z axis
% vas_med_f(i,j) = forz(vas_med_length(i,j), vas_med_mif, vas_med_ofl, vas_med_tsl, vas_med_pa); %Force
% vas_med_tqZ(i,j) = vas_med_f(i,j)*vas_med_maZ(i,j); %Torque, z axis
%Vastus Intermedius
% vas_int_o = [0.029; -0.192; 0.031]; %Femur, right, origin
% vas_int_wr1 = [0.034; -0.208; 0.028]; %Femur, right, wrapping point 1
% vas_int_wr2 = [0.054; -0.332; 0.013]; %Femur, right, wrapping point 2
% vas_int_i = [fcn7(theta_k); fcn8(theta_k); 0.002]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% vas_int_length(i,j) = norm(vas_int_wr1-vas_int_o)+norm(vas_int_wr1-vas_int_wr2)+norm(vas_int_wr2-VecTrans(T_k,vas_int_i)); %Total length, right
% vas_int_maZ(i,j) = CrossProd(vas_int_i, (VecTrans(T_k\eye(4),vas_int_wr1)-vas_int_i)/norm(VecTrans(T_k\eye(4),vas_int_wr1)-vas_int_i),3); %Moment arm, z axis
% vas_int_f(i,j) = forz(vas_int_length(i,j), vas_int_mif, vas_int_ofl, vas_int_tsl, vas_int_pa); %Force
% vas_int_tqZ(i,j) = vas_int_f(i,j)*vas_int_maZ(i,j); %Torque, z axis
%Vastus Lateralis
% vas_lat_o = [0.005; -0.185; 0.035]; %Femur, right, origin
% vas_lat_wr1 = [0.027; -0.259; 0.041]; %Femur, right, wrapping point 1
% vas_lat_wr2 = [0.054; -0.332; 0.029]; %Femur, right, wrapping point 2
% vas_lat_i = [fcn9(theta_k); fcn10(theta_k); 0.016]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% vas_lat_length(i,j) = norm(vas_lat_wr1-vas_lat_o)+norm(vas_lat_wr1-vas_lat_wr2)+norm(vas_lat_wr2-(VecTrans(T_k,vas_lat_i))); %Total length
% vas_lat_maZ(i,j) = CrossProd(vas_lat_i, (VecTrans(T_k\eye(4),vas_lat_wr1)-vas_lat_i)/norm(VecTrans(T_k\eye(4),vas_lat_wr1)-vas_lat_i),3); %Moment arm, z axis
% vas_lat_f(i,j) = forz(vas_lat_length(i,j), vas_lat_mif, vas_lat_ofl, vas_lat_tsl, vas_lat_pa); %Force
% vas_lat_tqZ(i,j) = vas_lat_f(i,j)*vas_lat_maZ(i,j); %Torque, z axis
% %Medial Gastroocnemius
% med_gas_length(i,j) = norm(VecTrans(T_k*T_a*T_s,med_gas_i)-med_gas_o); %Total length
% med_gas_maKZ(i,j) = CrossProd(VecTrans(T_a*T_s,med_gas_i), (VecTrans(T_k\eye(4),med_gas_o)-VecTrans(T_a*T_s,med_gas_i))/norm(VecTrans(T_k\eye(4),med_gas_o)-VecTrans(T_a*T_s,med_gas_i)),3); %Moment arm, Knee, z axis
% med_gas_maAZ(i,j) = CrossProd(VecTrans(T_s,med_gas_i), (VecTrans((T_k*T_a)\eye(4),med_gas_o)-VecTrans(T_s,med_gas_i))/norm(VecTrans((T_k*T_a)\eye(4),med_gas_o)-VecTrans(T_s,med_gas_i)),3,R_a); %Moment arm, Ankle, z' axis
% med_gas_maSX(i,j) = CrossProd(med_gas_i, (VecTrans((T_k*T_a*T_s)\eye(4),med_gas_o)-med_gas_i)/norm(VecTrans((T_k*T_a*T_s)\eye(4),med_gas_o)-med_gas_i),1,R_s); %Moment arm, Subtalar, x' axis
% med_gas_f(i,j) = forz(med_gas_length(i,j), med_gas_mif, med_gas_ofl, med_gas_tsl, med_gas_pa); %Force
% med_gas_tqKZ(i,j) = med_gas_f(i,j)*med_gas_maKZ(i,j); %Torque, Knee, z axis
% med_gas_tqAZ(i,j) = med_gas_f(i,j)*med_gas_maAZ(i,j); %Torque, Ankle, z' axis
% med_gas_tqSX(i,j) = med_gas_f(i,j)*med_gas_maSX(i,j); %Torque, Subtalar, x' axis
%Lateral Gastrocenemius
% lat_gas_length(i,j) = norm(VecTrans(T_k*T_a*T_s,lat_gas_i)-lat_gas_o); %Total length
% lat_gas_maKZ(i,j) = CrossProd(VecTrans(T_a*T_s,lat_gas_i), (VecTrans(T_k\eye(4),lat_gas_o)-VecTrans(T_a*T_s,lat_gas_i))/norm(VecTrans(T_k\eye(4),lat_gas_o)-VecTrans(T_a*T_s,lat_gas_i)),3); %Moment arm, Knee, z axis
% lat_gas_maAZ(i,j) = CrossProd(VecTrans(T_s,lat_gas_i), (VecTrans((T_k*T_a)\eye(4),lat_gas_o)-VecTrans(T_s,lat_gas_i))/norm(VecTrans((T_k*T_a)\eye(4),lat_gas_o)-VecTrans(T_s,lat_gas_i)),3,R_a); %Moment arm, Ankle, z' axis
% lat_gas_maSX(i,j) = CrossProd(lat_gas_i, (VecTrans((T_k*T_a*T_s)\eye(4),lat_gas_o)-lat_gas_i)/norm(VecTrans((T_k*T_a*T_s)\eye(4),lat_gas_o)-lat_gas_i),1,R_s); %Moment arm, Subtalar, x' axis
% lat_gas_f(i,j) = forz(lat_gas_length(i,j), lat_gas_mif, lat_gas_ofl, lat_gas_tsl, lat_gas_pa); %Force
% lat_gas_tqKZ(i,j) = lat_gas_f(i,j)*lat_gas_maKZ(i,j); %Torque, Knee, z axis
% lat_gas_tqAZ(i,j) = lat_gas_f(i,j)*lat_gas_maAZ(i,j); %Torque, Ankle, z' axis
% lat_gas_tqSX(i,j) = lat_gas_f(i,j)*lat_gas_maSX(i,j); %Torque, Subtalar, x' axis
%Soleus
% soleus_length(i,j) = norm(VecTrans(T_a*T_s,soleus_i)-soleus_o); %Total length
% soleus_maAZ(i,j) = CrossProd(VecTrans(T_s,soleus_i), (VecTrans(T_a\eye(4),soleus_o)-VecTrans(T_s,soleus_i))/norm(VecTrans(T_a\eye(4),soleus_o)-VecTrans(T_s,soleus_i)),3,R_a); %Moment arm, Ankle, z' axis
% soleus_maSX(i,j) = CrossProd(soleus_i, (VecTrans((T_a*T_s)\eye(4),soleus_o)-soleus_i)/norm(VecTrans((T_a*T_s)\eye(4),soleus_o)-soleus_i),1,R_s); %Moment arm, Subtalar, x' axis
% soleus_f(i,j) = forz(soleus_length(i,j), soleus_mif, soleus_ofl, lat_gas_tsl, lat_gas_pa); %Force
% soleus_tqAZ(i,j) = soleus_f(i,j)*soleus_maAZ(i,j); %Torque, Ankle, z' axis
% soleus_tqSX(i,j) = soleus_f(i,j)*soleus_maSX(i,j); %Torque, Subtalar, x' axis
%Tibialis Posterior
% tib_post_length(i,j) = norm(tib_post_wr1-tib_post_o)+norm(VecTrans(T_a*T_s,tib_post_wr2)-tib_post_wr1)+norm(tib_post_i-tib_post_wr2); %Total length
% tib_post_maAZ(i,j) = CrossProd(VecTrans(T_s,tib_post_wr2), (VecTrans(T_a\eye(4),tib_post_wr1)-VecTrans(T_s,tib_post_wr2))/norm(VecTrans(T_a\eye(4),tib_post_wr1)-VecTrans(T_s,tib_post_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% tib_post_maSX(i,j) = CrossProd(tib_post_wr2, (VecTrans((T_a*T_s)\eye(4),tib_post_wr1)-tib_post_wr2)/norm(VecTrans((T_a*T_s)\eye(4),tib_post_wr1)-tib_post_wr2),1,R_s); %Moment arm, subtalar, x' axis
% tib_post_f(i,j) = forz(tib_post_length(i,j), tib_post_mif, tib_post_ofl, tib_post_tsl, tib_post_pa); %Force
% tib_post_tqAZ(i,j) = tib_post_f(i,j)*tib_post_maAZ(i,j); %Torque, Ankle, z' axis
% tib_post_tqSX(i,j) = tib_post_f(i,j)*tib_post_maSX(i,j); %Torque, Subtalar, x' axis
%Flexor Digitorum Longus
% flex_dig_length(i,j) = norm(flex_dig_o-flex_dig_wr1)+norm(VecTrans(T_a*T_s,flex_dig_wr2)-flex_dig_wr1)+norm(flex_dig_wr2-flex_dig_wr3)+norm(flex_dig_wr3-flex_dig_wr4)+norm(VecTrans(T_m,flex_dig_wr5)-flex_dig_wr4)+norm(flex_dig_wr5-flex_dig_wr6)+norm(flex_dig_wr6-flex_dig_i); %Total length
% flex_dig_maAZ(i,j) = CrossProd(VecTrans(T_s,flex_dig_wr2), (VecTrans(T_a\eye(4),flex_dig_wr1)-VecTrans(T_s,flex_dig_wr2))/norm(VecTrans(T_a\eye(4),flex_dig_wr1)-VecTrans(T_s,flex_dig_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% flex_dig_maSX(i,j) = CrossProd(flex_dig_wr2, (VecTrans((T_a*T_s)\eye(4),flex_dig_wr1)-flex_dig_wr2)/norm(VecTrans((T_a*T_s)\eye(4),flex_dig_wr1)-flex_dig_wr2),1,R_s); %Moment arm, subtalar, x' axis
% flex_dig_maMZ(i,j) = CrossProd(flex_dig_wr5, (VecTrans(T_m\eye(4),flex_dig_wr4)-flex_dig_wr5)/norm(VecTrans(T_m\eye(4),flex_dig_wr4)-flex_dig_wr5),3,R_m); %Moment arm, MTP, z' axis
% flex_dig_f(i,j) = forz(flex_dig_length(i,j), flex_dig_mif, flex_dig_ofl, flex_dig_tsl, flex_dig_pa); %Force
% flex_dig_tqAZ(i,j) = flex_dig_f(i,j)*flex_dig_maAZ(i,j); %Torque, Ankle, z' axis
% flex_dig_tqSX(i,j) = flex_dig_f(i,j)*flex_dig_maSX(i,j); %Torque, subtalar, x' axis
% flex_dig_tqMZ(i,j) = flex_dig_f(i,j)*flex_dig_maMZ(i,j); %Torque, MTP, z' axis
%Flexor Hallucis Longus
% flex_hal_length(i,j) = norm(flex_hal_o-flex_hal_wr1)+norm(VecTrans(T_a*T_s,flex_hal_wr2)-flex_hal_wr1)+norm(flex_hal_wr2-flex_hal_wr3)+norm(flex_hal_wr3-flex_hal_wr4)+norm(VecTrans(T_m,flex_hal_wr5)-flex_hal_wr4)+norm(flex_hal_wr5-flex_hal_i); %Total length
% flex_hal_maAZ(i,j) = CrossProd(VecTrans(T_s,flex_hal_wr2), (VecTrans(T_a\eye(4),flex_hal_wr1)-VecTrans(T_s,flex_hal_wr2))/norm(VecTrans(T_a\eye(4),flex_hal_wr1)-VecTrans(T_s,flex_hal_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% flex_hal_maSX(i,j) = CrossProd(flex_hal_wr2, (VecTrans((T_a*T_s)\eye(4),flex_hal_wr1)-flex_hal_wr2)/norm(VecTrans((T_a*T_s)\eye(4),flex_hal_wr1)-flex_hal_wr2),1,R_s); %Moment arm, subtalar, x' axis
% flex_hal_maMZ(i,j) = CrossProd(flex_hal_wr5, (VecTrans(T_m\eye(4),flex_hal_wr4)-flex_hal_wr5)/norm(VecTrans(T_m\eye(4),flex_hal_wr4)-flex_hal_wr5),3,R_m); %Moment arm, MTP, z' axis
% flex_hal_f(i,j) = forz(flex_hal_length(i,j), flex_hal_mif, flex_hal_ofl, flex_hal_tsl, flex_hal_pa); %Force
% flex_hal_tqAZ(i,j) = flex_hal_f(i,j)*flex_hal_maAZ(i,j); %Torque, Ankle, z' axis
% flex_hal_tqSX(i,j) = flex_hal_f(i,j)*flex_hal_maSX(i,j); %Torque, subtalar, x' axis
% flex_hal_tqMZ(i,j) = flex_hal_f(i,j)*flex_hal_maMZ(i,j); %Torque, MTP, z' axis
%Tibialis Anterior
% tib_ant_length(i,j) = norm(tib_ant_wr1-tib_ant_o)+norm(VecTrans(T_a*T_s,tib_ant_i)-tib_ant_wr1); %Total length
% tib_ant_maAZ(i,j) = CrossProd(VecTrans(T_s,tib_ant_i), (VecTrans(T_a\eye(4),tib_ant_wr1)-VecTrans(T_s,tib_ant_i))/norm(VecTrans(T_a\eye(4),tib_ant_wr1)-VecTrans(T_s,tib_ant_i)),3,R_a); %Moment arm, Ankle, z' axis
% tib_ant_maSX(i,j) = CrossProd(tib_ant_i, (VecTrans((T_a*T_s)\eye(4),tib_ant_wr1)-tib_ant_i)/norm(VecTrans((T_a*T_s)\eye(4),tib_ant_wr1)-tib_ant_i),1,R_s); %Moment arm, subtalar, x' axis
% tib_ant_f(i,j) = forz(tib_ant_length(i,j), tib_ant_mif, tib_ant_ofl, tib_ant_tsl, tib_ant_pa); %Force
% tib_ant_tqAZ(i,j) = tib_ant_f(i,j)*tib_ant_maAZ(i,j); %Torque, Ankle, z' axis
% tib_ant_tqSX(i,j) = tib_ant_f(i,j)*tib_ant_maSX(i,j); %Torque, Subtalar, x' axis
%Peroneus Brevis
% per_brev_wr2 = [-0.014; -0.429; 0.029]; %Tibia, right, wrapping point 2
% per_brev_wr3 = [0.047; 0.027; 0.023]; %Calcn, right, wrapping point 3
% per_brev_length(i,j) = norm(per_brev_o-per_brev_wr1)+norm(per_brev_wr1-per_brev_wr2)+norm(VecTrans(T_a*T_s,per_brev_wr3)-per_brev_wr2)+norm(per_brev_wr3-per_brev_i); %Total length
% per_brev_maAZ(i,j) = CrossProd(VecTrans(T_s,per_brev_wr3), (VecTrans(T_a\eye(4),per_brev_wr2)-VecTrans(T_s,per_brev_wr3))/norm(VecTrans(T_a\eye(4),per_brev_wr2)-VecTrans(T_s,per_brev_wr3)),3,R_a); %Moment arm, Ankle, z' axis
% per_brev_maSX(i,j) = CrossProd(per_brev_wr3, (VecTrans((T_a*T_s)\eye(4),per_brev_wr2)-per_brev_wr3)/norm(VecTrans((T_a*T_s)\eye(4),per_brev_wr2)-per_brev_wr3),1,R_s); %Moment arm, subtalar, x' axis
% per_brev_f(i,j) = forz(per_brev_length(i,j), per_brev_mif, per_brev_ofl, per_brev_tsl, per_brev_pa); %Force
% per_brev_tqAZ(i,j) = per_brev_f(i,j)*per_brev_maAZ(i,j); %Torque, Ankle, z' axis
% per_brev_tqSX(i,j) = per_brev_f(i,j)*per_brev_maSX(i,j); %Torque, Subtalar, x' axis
%Peroneus Longus
% per_long_length(i,j) = norm(per_long_o-per_long_wr1)+norm(per_long_wr1-per_long_wr2)+norm(VecTrans(T_a*T_s,per_long_wr3)-per_long_wr2)+norm(per_long_wr3-per_long_wr4)+norm(per_long_wr5-per_long_wr4)+norm(per_long_wr5-per_long_i); %Total length
% per_long_maAZ(i,j) = CrossProd(VecTrans(T_s,per_long_wr3), (VecTrans(T_a\eye(4),per_long_wr2)-VecTrans(T_s,per_long_wr3))/norm(VecTrans(T_a\eye(4),per_long_wr2)-VecTrans(T_s,per_long_wr3)),3,R_a); %Moment arm, Ankle, z' axis
% per_long_maSX(i,j) = CrossProd(per_long_wr3, (VecTrans((T_a*T_s)\eye(4),per_long_wr2)-per_long_wr3)/norm(VecTrans((T_a*T_s)\eye(4),per_long_wr2)-per_long_wr3),1,R_s); %Moment arm, subtalar, x' axis
% per_long_f(i,j) = forz(per_long_length(i,j), per_long_mif, per_long_ofl, per_long_tsl, per_long_pa); %Force
% per_long_tqAZ(i,j) = per_long_f(i,j)*per_long_maAZ(i,j); %Torque, Ankle, z' axis
% per_long_tqSX(i,j) = per_long_f(i,j)*per_long_maSX(i,j); %Torque, Subtalar, x' axis
%Peroneus Tertius
% per_tert_length(i,j) = norm(per_tert_wr1-per_tert_o)+norm(VecTrans(T_a*T_s,per_tert_i)-per_tert_wr1) ; %Total length
% per_tert_maAZ(i,j) = CrossProd(VecTrans(T_s,per_tert_i), (VecTrans(T_a\eye(4),per_tert_wr1)-VecTrans(T_s,per_tert_i))/norm(VecTrans(T_a\eye(4),per_tert_wr1)-VecTrans(T_s,per_tert_i)),3,R_a); %Moment arm, Ankle, z' axis
% per_tert_maSX(i,j) = CrossProd(per_tert_i, (VecTrans((T_a*T_s)\eye(4),per_tert_wr1)-per_tert_i)/norm(VecTrans((T_a*T_s)\eye(4),per_tert_wr1)-per_tert_i),1,R_s); %Moment arm, subtalar, x' axis
% per_tert_f(i,j) = forz(per_tert_length(i,j), per_tert_mif, per_tert_ofl, per_tert_tsl, per_tert_pa); %Force
% per_tert_tqAZ(i,j) = per_tert_f(i,j)*per_tert_maAZ(i,j); %Torque, Ankle, z' axis
% per_tert_tqSX(i,j) = per_tert_f(i,j)*per_tert_maSX(i,j); %Torque, Subtalar, x' axis
%Extensor Digitorum Longus
% ext_dig_length(i,j) = norm(ext_dig_o-ext_dig_wr1)+norm(VecTrans(T_a*T_s,ext_dig_wr2)-ext_dig_wr1)+norm(ext_dig_wr2-ext_dig_wr3)+norm(VecTrans(T_m,ext_dig_wr4)-ext_dig_wr3)+norm(ext_dig_wr4-ext_dig_i); %Total length
% ext_dig_maAZ(i,j) = CrossProd(VecTrans(T_s,ext_dig_wr2), (VecTrans(T_a\eye(4),ext_dig_wr1)-VecTrans(T_s,ext_dig_wr2))/norm(VecTrans(T_a\eye(4),ext_dig_wr1)-VecTrans(T_s,ext_dig_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% ext_dig_maSX(i,j) = CrossProd(ext_dig_wr2, (VecTrans((T_a*T_s)\eye(4),ext_dig_wr1)-ext_dig_wr2)/norm(VecTrans((T_a*T_s)\eye(4),ext_dig_wr1)-ext_dig_wr2),1,R_s); %Moment arm, subtalar, x' axis
% ext_dig_maMZ(i,j) = CrossProd(ext_dig_wr4, (VecTrans(T_m\eye(4),ext_dig_wr3)-ext_dig_wr4)/norm(VecTrans(T_m\eye(4),ext_dig_wr3)-ext_dig_wr4),3,R_m); %Moment arm, MTP, z' axis
% ext_dig_f(i,j) = forz(ext_dig_length(i,j), ext_dig_mif, ext_dig_ofl, ext_dig_tsl, ext_dig_pa); %Force
% ext_dig_tqAZ(i,j) = ext_dig_f(i,j)*ext_dig_maAZ(i,j); %Torque, Ankle, z' axis
% ext_dig_tqSX(i,j) = ext_dig_f(i,j)*ext_dig_maSX(i,j); %Torque, subtalar, x' axis
% ext_dig_tqMZ(i,j) = ext_dig_f(i,j)*ext_dig_maMZ(i,j); %Torque, MTP, z' axis
%Extensor Hallucis Longus
% ext_hal_length(i,j) = norm(ext_hal_o-ext_hal_wr1)+norm(VecTrans(T_a*T_s,ext_hal_wr2)-ext_hal_wr1)+norm(ext_hal_wr2-ext_hal_wr3)+norm(ext_hal_wr3-ext_hal_wr4)+norm(VecTrans(T_m,ext_hal_wr5)-ext_hal_wr4)+norm(ext_hal_wr5-ext_hal_i); %Total length
% ext_hal_maAZ(i,j) = CrossProd(VecTrans(T_s,ext_hal_wr2), (VecTrans(T_a\eye(4),ext_hal_wr1)-VecTrans(T_s,ext_hal_wr2))/norm(VecTrans(T_a\eye(4),ext_hal_wr1)-VecTrans(T_s,ext_hal_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% ext_hal_maSX(i,j) = CrossProd(ext_hal_wr2, (VecTrans((T_a*T_s)\eye(4),ext_hal_wr1)-ext_hal_wr2)/norm(VecTrans((T_a*T_s)\eye(4),ext_hal_wr1)-ext_hal_wr2),1,R_s); %Moment arm, subtalar, x' axis
% ext_hal_maMZ(i,j) = CrossProd(ext_hal_wr5, (VecTrans(T_m\eye(4),ext_hal_wr4)-ext_hal_wr5)/norm(VecTrans(T_m\eye(4),ext_hal_wr4)-ext_hal_wr5),3,R_m); %Moment arm, MTP, z' axis
% ext_hal_f(i,j) = forz(ext_hal_length(i,j), ext_hal_mif, ext_hal_ofl, ext_hal_tsl, ext_hal_pa); %Force
% ext_hal_tqAZ(i,j) = ext_hal_f(i,j)*ext_hal_maAZ(i,j); %Torque, Ankle, z' axis
% ext_hal_tqSX(i,j) = ext_hal_f(i,j)*ext_hal_maSX(i,j); %Torque, subtalar, x' axis
% ext_hal_tqMZ(i,j) = ext_hal_f(i,j)*ext_hal_maMZ(i,j); %Torque, MTP, z' axis
%Erector Spinae
ercspn_length(i,j) = norm(VecTrans(T_b,ercspn_i) - ercspn_o); %Total length, right
ercspn_maX(i,j) = CrossProd(ercspn_i, (VecTrans(T_b\eye(4),ercspn_o)-ercspn_i)/norm(VecTrans(T_b\eye(4),ercspn_o)-ercspn_i),1); %Moment arm, x axis
ercspn_maY(i,j) = CrossProd(ercspn_i, (VecTrans(T_b\eye(4),ercspn_o)-ercspn_i)/norm(VecTrans(T_b\eye(4),ercspn_o)-ercspn_i),2); %Moment arm, y axis
ercspn_maZ(i,j) = CrossProd(ercspn_i, (VecTrans(T_b\eye(4),ercspn_o)-ercspn_i)/norm(VecTrans(T_b\eye(4),ercspn_o)-ercspn_i),3); %Moment arm, z axis
ercspn_f(i,j) = forz(ercspn_length(i,j), ercspn_mif, ercspn_ofl, ercspn_tsl, ercspn_pa); %Force
ercspn_tqX(i,j) = ercspn_f(i,j)*ercspn_maX(i,j); %Torque, x axis
ercspn_tqY(i,j) = ercspn_f(i,j)*ercspn_maY(i,j); %Torque, y axis
ercspn_tqZ(i,j) = ercspn_f(i,j)*ercspn_maZ(i,j); %Torque, z axis
%Internal Oblique
% intobl_length(i,j) = norm(VecTrans(T_b,intobl_i) - intobl_o); %Total length, right
% intobl_maX(i,j) = CrossProd(intobl_i, (VecTrans(T_b\eye(4),intobl_o)-intobl_i)/norm(VecTrans(T_b\eye(4),intobl_o)-intobl_i),1); %Moment arm, x axis
% intobl_maY(i,j) = CrossProd(intobl_i, (VecTrans(T_b\eye(4),intobl_o)-intobl_i)/norm(VecTrans(T_b\eye(4),intobl_o)-intobl_i),2); %Moment arm, y axis
% intobl_maZ(i,j) = CrossProd(intobl_i, (VecTrans(T_b\eye(4),intobl_o)-intobl_i)/norm(VecTrans(T_b\eye(4),intobl_o)-intobl_i),3); %Moment arm, z axis
% intobl_f(i,j) = forz(intobl_length(i,j), intobl_mif, intobl_ofl, intobl_tsl, intobl_pa); %Force
% intobl_tqX(i,j) = intobl_f(i,j)*intobl_maX(i,j); %Torque, x axis
% intobl_tqY(i,j) = intobl_f(i,j)*intobl_maY(i,j); %Torque, y axis
% intobl_tqZ(i,j) = intobl_f(i,j)*intobl_maZ(i,j); %Torque, z axis
%External Oblique extobl
% extobl_length(i,j) = norm(VecTrans(T_b,extobl_i) - extobl_o); %Total length, right
% extobl_maX(i,j) = CrossProd(extobl_i, (VecTrans(T_b\eye(4),extobl_o)-extobl_i)/norm(VecTrans(T_b\eye(4),extobl_o)-extobl_i),1); %Moment arm, x axis
% extobl_maY(i,j) = CrossProd(extobl_i, (VecTrans(T_b\eye(4),extobl_o)-extobl_i)/norm(VecTrans(T_b\eye(4),extobl_o)-extobl_i),2); %Moment arm, y axis
% extobl_maZ(i,j) = CrossProd(extobl_i, (VecTrans(T_b\eye(4),extobl_o)-extobl_i)/norm(VecTrans(T_b\eye(4),extobl_o)-extobl_i),3); %Moment arm, z axis
% extobl_f(i,j) = forz(extobl_length(i,j), extobl_mif, extobl_ofl, extobl_tsl, extobl_pa); %Force
% extobl_tqX(i,j) = extobl_f(i,j)*extobl_maX(i,j); %Torque, x axis
% extobl_tqY(i,j) = extobl_f(i,j)*extobl_maY(i,j); %Torque, y axis
% extobl_tqZ(i,j) = extobl_f(i,j)*extobl_maZ(i,j); %Torque, z axis
end
end
%%
%%Plots and Surface figures
% figure
% plot(theta_K, rect_fem_length)
% title('Rect Fem length for knee flexion')
% ylabel('Length, m')
% xlabel('Flexion, Degrees')

% uni_hip_tqZ = glut_min1_tqZ+glut_min2_tqZ+glut_min3_tqZ+glut_med1_tqZ+glut_med2_tqZ+glut_med3_tqZ+glut_max1_tqZ+glut_max2_tqZ+glut_max3_tqZ+add_long_tqZ+add_brev_tqZ+add_mag1_tqZ+add_mag2_tqZ+add_mag3_tqZ+pect_tqZ+iliacus_tqZ+psoas_tqHZ+quad_fem_tqZ+gem_tqZ+peri_tqZ;
% uni_hip_tqY = glut_min1_tqY+glut_min2_tqY+glut_min3_tqY+glut_med1_tqY+glut_med2_tqY+glut_med3_tqY+glut_max1_tqY+glut_max2_tqY+glut_max3_tqY+add_long_tqY+add_brev_tqY+add_mag1_tqY+add_mag2_tqY+add_mag3_tqY+pect_tqY+iliacus_tqY+psoas_tqHY+quad_fem_tqY+gem_tqY+peri_tqY;
% uni_hip_tqX = glut_min1_tqX+glut_min2_tqX+glut_min3_tqX+glut_med1_tqX+glut_med2_tqX+glut_med3_tqX+glut_max1_tqX+glut_max2_tqX+glut_max3_tqX+add_long_tqX+add_brev_tqX+add_mag1_tqX+add_mag2_tqX+add_mag3_tqX+pect_tqX+iliacus_tqX+psoas_tqHX+quad_fem_tqX+gem_tqX+peri_tqX;

% bi_hip_tqHX = semimem_tqHX+semiten_tqHX+bifemlh_tqHX+sar_tqHX+tfl_tqHX+rect_fem_tqHX+grac_tqHX;
% bi_hip_tqHY = semimem_tqHY+semiten_tqHY+bifemlh_tqHY+sar_tqHY+tfl_tqHY+rect_fem_tqHY+grac_tqHY;
% bi_hip_tqHZ = semimem_tqHZ+semiten_tqHZ+bifemlh_tqHZ+sar_tqHZ+tfl_tqHZ+rect_fem_tqHZ+grac_tqHZ;
% bi_knee_tqKZ = semimem_tqKZ+semiten_tqKZ+bifemlh_tqKZ+sar_tqKZ+tfl_tqKZ+rect_fem_tqKZ+grac_tqKZ;
% uni_knee_tqZ = bifemsh_tqZ+vas_int_tqZ+vas_med_tqZ+vas_lat_tqZ;
% vas_tqZ = vas_int_tqZ+vas_med_tqZ+vas_lat_tqZ;

% calves_tqKZ = med_gas_tqKZ + lat_gas_tqKZ;
% calves_tqAZ = med_gas_tqAZ + lat_gas_tqAZ + soleus_tqAZ;
% calves_tqSX = med_gas_tqSX + lat_gas_tqSX + soleus_tqSX;

% toes_tqAZ = flex_dig_tqAZ + flex_hal_tqAZ +ext_dig_tqAZ + ext_hal_tqAZ;
% toes_tqSX = flex_dig_tqSX + flex_hal_tqSX +ext_dig_tqSX + ext_hal_tqSX;
% toes_tqMZ = flex_dig_tqMZ + flex_hal_tqMZ +ext_dig_tqMZ + ext_hal_tqMZ;

% foot_tqAZ = tib_post_tqAZ + tib_ant_tqAZ + per_brev_tqAZ + per_long_tqAZ + per_tert_tqAZ;
% foot_tqSX = tib_post_tqSX + tib_ant_tqSX + per_brev_tqSX + per_long_tqSX + per_tert_tqSX;

figure
surf(theta_bZ, theta_bX, ercspn_tqZ, 'EdgeColor','none')
title('Back Torque, ercspn, Z axis')
ylabel('Lumbar Bending, Degrees')
xlabel('Back Flexion, Degrees')
zlabel('Torque, N*m')
colorbar

figure
surf(theta_bZ, theta_bX, ercspn_tqY, 'EdgeColor','none')
title('Back Torque, ercspn, Y axis')
ylabel('Lumbar Bending, Degrees')
xlabel('Back Flexion, Degrees')
zlabel('Torque, N*m')
colorbar

figure
surf(theta_bZ, theta_bX, ercspn_tqX, 'EdgeColor','none')
title('Back Torque, ercspn, X axis')
ylabel('Lumbar Bending, Degrees')
xlabel('Back Flexion, Degrees')
zlabel('Torque, N*m')
colorbar
% 
%%
%%Size Calculations
% [glut_med1_dia, glut_med1_size]= Size(glut_med1_length, glut_med1_mif)
% [glut_med2_dia, glut_med2_size]= Size(glut_med2_length, glut_med2_mif)
% [glut_med3_dia, glut_med3_size]= Size(glut_med3_length, glut_med3_mif)
% [glut_min1_dia, glut_min1_size]= Size(glut_min1_length, glut_min1_mif)
% [glut_min2_dia, glut_min2_size]= Size(glut_min2_length, glut_min2_mif)
% [glut_min3_dia, glut_min3_size]= Size(glut_min3_length, glut_min3_mif)
% [glut_max1_dia, glut_max1_size]= Size(glut_max1_length, glut_max1_mif)
% [glut_max2_dia, glut_max2_size]= Size(glut_max2_length, glut_max2_mif)
% [glut_max3_dia, glut_max3_size]= Size(glut_max3_length, glut_max3_mif)
% [semimem_dia, semimem_size]= Size(semimem_length, semimem_mif)
% [semiten_dia, semiten_size]= Size(semiten_length, semiten_mif)
% [bifemlh_dia, bifemlh_size]= Size(bifemlh_length, bifemlh_mif)
% [bifemsh_dia, bifemsh_size]= Size(bifemsh_length, bifemsh_mif)
% [sar_dia, sar_size]= Size(sar_length, sar_mif)
% [add_long_dia, add_long_size]= Size(add_long_length, add_long_mif)
% [add_brev_dia, add_brev_size]= Size(add_brev_length, add_brev_mif)
% [add_mag1_dia, add_mag1_size]= Size(add_mag1_length, add_mag1_mif)
% [add_mag2_dia, add_mag2_size]= Size(add_mag2_length, add_mag2_mif)
% [add_mag3_dia, add_mag3_size]= Size(add_mag3_length, add_mag3_mif)
% [tfl_dia, tfl_size]= Size(tfl_length, tfl_mif)
% [pect_dia, pect_size]= Size(pect_length, pect_mif)
% [grac_dia, grac_size]= Size(grac_length, grac_mif)
% [iliacus_dia, iliacus_size]= Size(iliacus_length, iliacus_mif)
% [psoas_dia, psoas_size]= Size(psoas_length, psoas_mif)
% [quad_fem_dia, quad_fem_size]= Size(quad_fem_length, quad_fem_mif)
% [gem_dia, gem_size]= Size(gem_length, gem_mif)
% [peri_dia, peri_size]= Size(peri_length, peri_mif)
% [rect_fem_dia, rect_fem_size]= Size(rect_fem_length, rect_fem_mif)
% [vas_med_dia, vas_med_size]= Size(vas_med_length, vas_med_mif)
% [vas_int_dia, vas_int_size]= Size(vas_int_length, vas_int_mif)
% [vas_lat_dia, vas_lat_size]= Size(vas_lat_length, vas_lat_mif)
% [med_gas_dia, med_gas_size]= Size(med_gas_length, med_gas_mif)
% [lat_gas_dia, lat_gas_size]= Size(lat_gas_length, lat_gas_mif)
% [soleus_dia, soleus_size]= Size(soleus_length, soleus_mif)
% [tib_post_dia, tib_post_size]= Size(tib_post_length, tib_post_mif)
% [flex_dig_dia, flex_dig_size]= Size(flex_dig_length, flex_dig_mif)
% [flex_hal_dia, flex_hal_size]= Size(flex_hal_length, flex_hal_mif)
% [tib_ant_dia, tib_ant_size]= Size(tib_ant_length, tib_ant_mif)
% [per_brev_dia, per_brev_size]= Size(per_brev_length, per_brev_mif)
% [per_long_dia, per_long_size]= Size(per_long_length, per_long_mif)
% [per_tert_dia, per_tert_size]= Size(per_tert_length, per_tert_mif)
% [ext_dig_dia, ext_dig_size]= Size(ext_dig_length, ext_dig_mif)
% [ext_hal_dia, ext_hal_size]= Size(ext_hal_length, ext_hal_mif)
% [ercspn_dia, ercspn_size]= Size(ercspn_length, ercspn_mif)
% [intobl_dia, intobl_size]= Size(intobl_length, intobl_mif)
% [extobl_dia, extobl_size]= Size(extobl_length, extobl_mif)

