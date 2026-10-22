%%Initialize
clear; clc; %close all;
% Robotic Legs Kinematics
%Units in m, N, rad/s unless otherwise stated.

load AttachPoints_Robot.mat  %Attachment points for muscles with stationary attachment points
load MovePoints.mat %Moving attachment points and knee data
load PAM_Prop.mat      %Muscle specific properties
IndexPrep %Matrices for length, force, and moment arm  
k = 180/pi;  %convert radians to degrees
c = pi/180;  %convert degrees to radians

%Range of motion, radians
range_hz = c*[ -15, 110];   %Hip (Wikipedia)
thetamin_hz = range_hz(1);  %15 deg extention
thetamax_hz = range_hz(2);  %110 deg flexion
range_k = c*[-120, 10];     %Knee
thetamin_k = range_k(1);    %120 deg flexion
thetamax_k = range_k(2);    %10 deg extension
thetamin_a = -50*c;         %ankle
thetamax_a = 20*c;
thetamin_s = -25*c;         %subtalar
thetamax_s = 35*c;
thetamin_m = -30*c;         %mtp
thetamax_m = 80*c;
thetamin_bx = -20*c;        %Lumbar bending  (I.A. Kapandji, Phsiology of the Joints. Churchill Livingstone, 1986)
thetamax_bx = 20*c;
thetamin_by = -5*c;         %Lumbar rotation  (Kapandji 1986)
thetamax_by = 5*c;
thetamin_bz = -60*c;        %Back flexion-/extension+  (Kapandji 1986)
thetamax_bz = 20*c;


%Indexing prep, see IndexPrep.m
s = 100;
t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)

%% T matrix, length, and moment arm calculations
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
theta_b_x = 0; %thetamin_bx + (thetamax_bx-thetamin_bx)*(i-1)/s; %Lumbar bending left-/right+
theta_b_y = 0; %thetamin_by + (thetamax_by-thetamin_by)*(j-1)/s; %Lumbar rotation cw-/ccw+
theta_b_z = 0; %thetamin_bz + (thetamax_bz-thetamin_bz)*(j-1)/s; %Back flexion-/extension+
theta_a = 0; %Ankle, right, plantar flexion-/dorsiflexion+
theta_s = 0; %Subtalar, right, inversion-/eversion+ rotation
theta_m = 0; %Metatarsophalangeal (mtp), right, flexion-/extension+

%Index for plotting (for 2 dof, j index to x axis, i index to y axis)
theta_hZ(j,1) = theta_h_z*k; %Convert from radians to degrees
%theta_hY(i,j) = theta_h_y*k;
% theta_hX(i,j) = theta_h_x*k;
theta_K(i,1) = theta_k*k;
%theta_bX(i,1) = theta_b_x*k;
%theta_bY(j,1) = theta_b_y*k;
%theta_bZ(j,1) = theta_b_z*k;
%theta_A(j,1) = theta_a*k;
%theta_S(i,1) = theta_s*k;
%theta_M(i,1) = theta_m*k;

%Rotation matrices
%Hip
R_h_x = [1, 0, 0; 0, cos(theta_h_x), -sin(theta_h_x); 0, sin(theta_h_x), cos(theta_h_x)];
R_h_y = [cos(theta_h_y), 0, sin(theta_h_y); 0, 1, 0; -sin(theta_h_y), 0, cos(theta_h_y)];
R_h_z = [cos(theta_h_z), -sin(theta_h_z), 0; sin(theta_h_z), cos(theta_h_z), 0; 0, 0, 1]; 
R_h = R_h_x*R_h_y*R_h_z;
%Knee (Patella is the same)
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
p_pa = [0.25; 0.3; 0]; %Patella (Guess)
p_b = [-0.1007; 0.0815; 0]; %Back
p_a = [0; -0.43; 0]; %Ankle, right
p_s = [-0.04877; -0.04195; 0.00792]; %Subtalar, right
p_m = [0.1788; -0.002; 0.00108]; %mtp, right

%Transformation matrices
T_p = [1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1]; %Pelvis
T_h = RpToTrans(R_h,p_h); %Hip to Pelvis, right
T_k = RpToTrans(R_k,p_k); %Knee to Hip, right
T_pa = RpToTrans(R_k,p_pa); %Patella to knee, right
T_b = RpToTrans(R_b,p_b); %Back to Pelvis
T_a = RpToTrans(R_a,p_a); %Ankle to Knee, right
T_s = RpToTrans(R_s,p_s); %Subtalar to Ankle, right
T_m = RpToTrans(R_m,p_m); %Mtp to Subtalar, right

%Muscle length and moment arm
%Gluteus Maximus
% glut_max_length(i,j) = norm(glut_max_wr1 - p6) + norm(VecTrans(T_h,glut_max_wr2) - glut_max_wr1) + norm(glut_max_wr3 - glut_max_wr2) + norm(f3 - glut_max_wr3); %Total length
% glut_max_maX(i,j) = CrossProd(glut_max_wr2, (VecTrans(T_h\eye(4),glut_max_wr1)-glut_max_wr2)/norm(VecTrans(T_h\eye(4),glut_max_wr1)-glut_max_wr2),1); %Moment arm, x axis
% glut_max_maY(i,j) = CrossProd(glut_max_wr2, (VecTrans(T_h\eye(4),glut_max_wr1)-glut_max_wr2)/norm(VecTrans(T_h\eye(4),glut_max_wr1)-glut_max_wr2),2); %Moment arm, y axis
% glut_max_maZ(i,j) = CrossProd(glut_max_wr2, (VecTrans(T_h\eye(4),glut_max_wr1)-glut_max_wr2)/norm(VecTrans(T_h\eye(4),glut_max_wr1)-glut_max_wr2),3); %Moment arm, z axis

%Semimembranosus
% semimem_length(i,j) = norm(VecTrans(T_h,semimem_wr1)-p5)+norm(semimem_wr2-semimem_wr1)+norm(VecTrans(T_h,twr1)-semimem_wr1)+norm(t1-twr1)+norm(t6-t1); %Total length
% semimem_maHX(i,j) = CrossProd(semimem_wr1, (VecTrans(T_h\eye(4),p5)-semimem_wr1)/norm(VecTrans(T_h\eye(4),p5)-semimem_wr1),1); %Moment arm, Hip, x axis
% semimem_maHY(i,j) = CrossProd(semimem_wr1, (VecTrans(T_h\eye(4),p5)-semimem_wr1)/norm(VecTrans(T_h\eye(4),p5)-semimem_wr1),2); %Moment arm, Hip, y axis
% semimem_maHZ(i,j) = CrossProd(semimem_wr1, (VecTrans(T_h\eye(4),p5)-semimem_wr1)/norm(VecTrans(T_h\eye(4),p5)-semimem_wr1),3); %Moment arm, Hip, z axis
% semimem_maKZ(i,j) = CrossProd(twr1, (VecTrans(T_k\eye(4),semimem_wr2)-twr1)/norm(VecTrans(T_k\eye(4),semimem_wr2)-twr1),3); %Moment arm, Knee, z axis

%Bicep Femoris, Long Head
% bifemlh_length(i,j) = norm(VecTrans(T_h*T_k,twr2)-p5)+norm(twr3-twr2)+norm(twr3-t5); %Total length
% bifemlh_maHX(i,j) = CrossProd(VecTrans(T_k,twr2), (VecTrans(T_h\eye(4),p5)-VecTrans(T_k,twr2))/norm(VecTrans(T_h\eye(4),p5)-VecTrans(T_k,twr2)),1); %Moment arm, Hip, x axis
% bifemlh_maHY(i,j) = CrossProd(VecTrans(T_k,twr2), (VecTrans(T_h\eye(4),p5)-VecTrans(T_k,twr2))/norm(VecTrans(T_h\eye(4),p5)-VecTrans(T_k,twr2)),2); %Moment arm, Hip, y axis
% bifemlh_maHZ(i,j) = CrossProd(VecTrans(T_k,twr2), (VecTrans(T_h\eye(4),p5)-VecTrans(T_k,twr2))/norm(VecTrans(T_h\eye(4),p5)-VecTrans(T_k,twr2)),3); %Moment arm, Hip, z axis
% bifemlh_maKZ(i,j) = CrossProd(twr2, (VecTrans((T_h*T_k)\eye(4),p5)-twr2)/norm(VecTrans((T_h*T_k)\eye(4),p5)-twr2),3); %Moment arm, Knee, z axis

%Bicep Femoris, Short Head
bifemsh_length(i,j) = norm(VecTrans(T_k,twr2)-f2)+norm(twr3-twr2); %Total length, right
bifemsh_maZ(i,j) = CrossProd(twr2, (VecTrans(T_k\eye(4),f2)-twr2)/norm(VecTrans(T_k\eye(4),f2)-twr2),3); %Moment arm, Knee, z axis

%Sartorius
% sar_length(i,j) = norm(VecTrans(T_h,sar_wr1)-p3)+norm(sar_wr1-sar_wr2)+norm(sar_wr3-sar_wr2)+norm(VecTrans(T_k,twr1)-sar_wr3)+norm(t1-twr1)+norm(t6-t1); %Total length, right
% sar_maHX(i,j) = CrossProd(sar_wr1, (VecTrans(T_h\eye(4),p3)-sar_wr1)/norm(VecTrans(T_h\eye(4),p3)-sar_wr1),1); %Moment arm, Hip, x axis
% sar_maHY(i,j) = CrossProd(sar_wr1, (VecTrans(T_h\eye(4),p3)-sar_wr1)/norm(VecTrans(T_h\eye(4),p3)-sar_wr1),2); %Moment arm, Hip, y axis
% sar_maHZ(i,j) = CrossProd(sar_wr1, (VecTrans(T_h\eye(4),p3)-sar_wr1)/norm(VecTrans(T_h\eye(4),p3)-sar_wr1),3); %Moment arm, Hip, z axis
% sar_maKZ(i,j) = CrossProd(twr1, (VecTrans(T_k\eye(4),sar_wr3)-twr1)/norm(VecTrans(T_k\eye(4),sar_wr3)-twr1),3); %Moment arm, Knee, z axis

%Adductor Magnus
% add_mag_length(i,j) = norm(p9 - add_mag_wr1) + norm(VecTrans(T_h,add_mag_wr2) - add_mag_wr1)+norm(add_mag_wr2-f4); %Total length
% add_mag_maX(i,j) = CrossProd(add_mag_wr2, (VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2)/norm(VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2),1); %Moment arm, x axis
% add_mag_maY(i,j) = CrossProd(add_mag_wr2, (VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2)/norm(VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2),2); %Moment arm, y axis
% add_mag_maZ(i,j) = CrossProd(add_mag_wr2, (VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2)/norm(VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2),3); %Moment arm, z axis

%Tensor Fasciae Latae
% tfl_length(i,j) = norm(VecTrans(T_h,tfl_wr1)-p3)+norm(tfl_wr2 - tfl_wr1)+norm(tfl_wr2-VecTrans(T_k,t2))+norm(t2-t5); %Total length, right
% tfl_maHX(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),p3)-tfl_wr1)/norm(VecTrans(T_h\eye(4),p3)-tfl_wr1),1); %Moment arm, Hip, x axis
% tfl_maHY(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),p3)-tfl_wr1)/norm(VecTrans(T_h\eye(4),p3)-tfl_wr1),2); %Moment arm, Hip, y axis
% tfl_maHZ(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),p3)-tfl_wr1)/norm(VecTrans(T_h\eye(4),p3)-tfl_wr1),3); %Moment arm, Hip, z axis
% tfl_maKZ(i,j) = CrossProd(t2, (VecTrans(T_k\eye(4),tfl_wr2)-t2)/norm(VecTrans(T_k\eye(4),tfl_wr2)-t2),3); %Moment arm, Knee, z axis

%Gracilis
% grac_length(i,j) = norm(VecTrans(T_h,grac_wr1)-p8)+norm(grac_wr2-grac_wr1)+norm(VecTrans(T_h,twr1)-grac_wr1)+norm(t1-twr1)+norm(t6-t1); %Total length
% grac_maHX(i,j) = CrossProd(grac_wr1, (VecTrans(T_h\eye(4),p8)-grac_wr1)/norm(VecTrans(T_h\eye(4),p8)-grac_wr1),1); %Moment arm, Hip, x axis
% grac_maHY(i,j) = CrossProd(grac_wr1, (VecTrans(T_h\eye(4),p8)-grac_wr1)/norm(VecTrans(T_h\eye(4),p8)-grac_wr1),2); %Moment arm, Hip, y axis
% grac_maHZ(i,j) = CrossProd(grac_wr1, (VecTrans(T_h\eye(4),p8)-grac_wr1)/norm(VecTrans(T_h\eye(4),p8)-grac_wr1),3); %Moment arm, Hip, z axis
% grac_maKZ(i,j) = CrossProd(twr1, (VecTrans(T_k\eye(4),grac_wr2)-twr1)/norm(VecTrans(T_k\eye(4),grac_wr2)-twr1),3); %Moment arm, Knee, z axis

%Iliacus
% iliacus_length(i,j) = norm(iliacus_wr1-p1)+norm(VecTrans(T_h,f4)-iliacus_wr1); %Total length
% iliacus_maX(i,j) = CrossProd(f4, (VecTrans(T_h\eye(4),iliacus_wr1)-f4)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-f4),1); %Moment arm, x axis
% iliacus_maY(i,j) = CrossProd(f4, (VecTrans(T_h\eye(4),iliacus_wr1)-f4)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-f4),2); %Moment arm, y axis
% iliacus_maZ(i,j) = CrossProd(f4, (VecTrans(T_h\eye(4),iliacus_wr1)-f4)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-f4),3); %Moment arm, z axis

%Psoas
% psoas_length(i,j) = norm(VecTrans(T_b,b1)-psoas_wr1)+norm(VecTrans(T_h,f4)-psoas_wr1); %Total length, right
% psoas_maHX(i,j) = CrossProd(f4, (VecTrans(T_h\eye(4),psoas_wr1)-f4)/norm(VecTrans(T_h\eye(4),psoas_wr1)-f4),1); %Moment arm, Hip, x axis
% psoas_maHY(i,j) = CrossProd(f4, (VecTrans(T_h\eye(4),psoas_wr1)-f4)/norm(VecTrans(T_h\eye(4),psoas_wr1)-f4),2); %Moment arm, Hip, y axis
% psoas_maHZ(i,j) = CrossProd(f4, (VecTrans(T_h\eye(4),psoas_wr1)-f4)/norm(VecTrans(T_h\eye(4),psoas_wr1)-f4),3); %Moment arm, Hip, z axis
% psoas_maBX(i,j) = CrossProd(b1, (VecTrans(T_b\eye(4),psoas_wr1)-b1)/norm(VecTrans(T_b\eye(4),psoas_wr1)-b1),1); %Moment arm, Back, x axis
% psoas_maBY(i,j) = CrossProd(b1, (VecTrans(T_b\eye(4),psoas_wr1)-b1)/norm(VecTrans(T_b\eye(4),psoas_wr1)-b1),2); %Moment arm, Back, y axis
% psoas_maBZ(i,j) = CrossProd(b1, (VecTrans(T_b\eye(4),psoas_wr1)-b1)/norm(VecTrans(T_b\eye(4),psoas_wr1)-b1),3); %Moment arm, Back, z axis

%Rectus Femoris
% pa1 = [fcn3(theta_k); fcn4(theta_k); 0.002]; %Tibia, from MovePoints, based on OpenSim data
% rect_fem_length(i,j) = norm(VecTrans(T_h,rect_fem_wr1)-p10)+norm(rect_fem_wr2-rect_fem_wr1)+norm(rect_fem_wr2-rect_fem_wr3)+norm(VecTrans(T_k,pa1)-rect_fem_wr3); %Total length
% rect_fem_maHX(i,j) = CrossProd(rect_fem_wr1, (VecTrans(T_h\eye(4),p10)-rect_fem_wr1)/norm(VecTrans(T_h\eye(4),p10)-rect_fem_wr1),1); %Moment arm, Hip, x axis
% rect_fem_maHY(i,j) = CrossProd(rect_fem_wr1, (VecTrans(T_h\eye(4),p10)-rect_fem_wr1)/norm(VecTrans(T_h\eye(4),p10)-rect_fem_wr1),2); %Moment arm, Hip, y axis
% rect_fem_maHZ(i,j) = CrossProd(rect_fem_wr1, (VecTrans(T_h\eye(4),p10)-rect_fem_wr1)/norm(VecTrans(T_h\eye(4),p10)-rect_fem_wr1),3); %Moment arm, Hip, z axis
% rect_fem_maKZ(i,j) = CrossProd(pa1, (VecTrans(T_k\eye(4),rect_fem_wr3)-pa1)/norm(VecTrans(T_k\eye(4),rect_fem_wr3)-pa1),3); %Moment arm, Knee, z axis

%Vastus Intermedius
pa2 = [fcn7(theta_k); fcn8(theta_k); 0.002]; %From MovePoints, based on OpenSim data
vas_int_length(i,j) = norm(vas_int_wr1-f1)+norm(vas_int_wr2-vas_int_wr1)+norm(vas_int_wr2-VecTrans(T_k,pa2)); %Total length, right
vas_int_maZ(i,j) = CrossProd(pa2, (VecTrans(T_k\eye(4),vas_int_wr2)-pa2)/norm(VecTrans(T_k\eye(4),vas_int_wr2)-pa2),3); %Moment arm, z axis

%Medial Gastroocnemius
% med_gas_length(i,j) = norm(VecTrans(T_k*T_a*T_s,a1)-f4); %Total length
% med_gas_maKZ(i,j) = CrossProd(VecTrans(T_a*T_s,a1), (VecTrans(T_k\eye(4),f4)-VecTrans(T_a*T_s,a1))/norm(VecTrans(T_k\eye(4),f4)-VecTrans(T_a*T_s,a1)),3); %Moment arm, Knee, z axis
% med_gas_maAZ(i,j) = CrossProd(VecTrans(T_s,a1), (VecTrans((T_k*T_a)\eye(4),f4)-VecTrans(T_s,a1))/norm(VecTrans((T_k*T_a)\eye(4),f4)-VecTrans(T_s,a1)),3,R_a); %Moment arm, Ankle, z' axis
% med_gas_maSX(i,j) = CrossProd(a1, (VecTrans((T_k*T_a*T_s)\eye(4),f4)-a1)/norm(VecTrans((T_k*T_a*T_s)\eye(4),f4)-a1),1,R_s); %Moment arm, Subtalar, x' axis

%Lateral Gastrocenemius
% lat_gas_length(i,j) = norm(VecTrans(T_k*T_a*T_s,a1)-f5); %Total length
% lat_gas_maKZ(i,j) = CrossProd(VecTrans(T_a*T_s,a1), (VecTrans(T_k\eye(4),f5)-VecTrans(T_a*T_s,a1))/norm(VecTrans(T_k\eye(4),f5)-VecTrans(T_a*T_s,a1)),3); %Moment arm, Knee, z axis
% lat_gas_maAZ(i,j) = CrossProd(VecTrans(T_s,a1), (VecTrans((T_k*T_a)\eye(4),f5)-VecTrans(T_s,a1))/norm(VecTrans((T_k*T_a)\eye(4),f5)-VecTrans(T_s,a1)),3,R_a); %Moment arm, Ankle, z' axis
% lat_gas_maSX(i,j) = CrossProd(a1, (VecTrans((T_k*T_a*T_s)\eye(4),f5)-a1)/norm(VecTrans((T_k*T_a*T_s)\eye(4),f5)-a1),1,R_s); %Moment arm, Subtalar, x' axis

%Soleus
% soleus_length(i,j) = norm(VecTrans(T_a*T_s,a1)-t1); %Total length
% soleus_maAZ(i,j) = CrossProd(VecTrans(T_s,a1), (VecTrans(T_a\eye(4),t1)-VecTrans(T_s,a1))/norm(VecTrans(T_a\eye(4),t1)-VecTrans(T_s,a1)),3,R_a); %Moment arm, Ankle, z' axis
% soleus_maSX(i,j) = CrossProd(a1, (VecTrans((T_a*T_s)\eye(4),t1)-a1)/norm(VecTrans((T_a*T_s)\eye(4),t1)-a1),1,R_s); %Moment arm, Subtalar, x' axis

%Tibialis Posterior
% tib_post_length(i,j) = norm(tib_post_wr1-t4)+norm(VecTrans(T_a*T_s,tib_post_wr2)-tib_post_wr1)+norm(a4-tib_post_wr2); %Total length
% tib_post_maAZ(i,j) = CrossProd(VecTrans(T_s,tib_post_wr2), (VecTrans(T_a\eye(4),tib_post_wr1)-VecTrans(T_s,tib_post_wr2))/norm(VecTrans(T_a\eye(4),tib_post_wr1)-VecTrans(T_s,tib_post_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% tib_post_maSX(i,j) = CrossProd(tib_post_wr2, (VecTrans((T_a*T_s)\eye(4),tib_post_wr1)-tib_post_wr2)/norm(VecTrans((T_a*T_s)\eye(4),tib_post_wr1)-tib_post_wr2),1,R_s); %Moment arm, subtalar, x' axis

%Flexor Digitorum Longus
% flex_dig_length(i,j) = norm(t4-flex_dig_wr1)+norm(VecTrans(T_a*T_s,flex_dig_wr2)-flex_dig_wr1)+norm(flex_dig_wr2-flex_dig_wr3)+norm(flex_dig_wr3-flex_dig_wr4)+norm(VecTrans(T_m,flex_dig_wr5)-flex_dig_wr4)+norm(flex_dig_wr5-flex_dig_wr6)+norm(flex_dig_wr6-a7); %Total length
% flex_dig_maAZ(i,j) = CrossProd(VecTrans(T_s,flex_dig_wr2), (VecTrans(T_a\eye(4),flex_dig_wr1)-VecTrans(T_s,flex_dig_wr2))/norm(VecTrans(T_a\eye(4),flex_dig_wr1)-VecTrans(T_s,flex_dig_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% flex_dig_maSX(i,j) = CrossProd(flex_dig_wr2, (VecTrans((T_a*T_s)\eye(4),flex_dig_wr1)-flex_dig_wr2)/norm(VecTrans((T_a*T_s)\eye(4),flex_dig_wr1)-flex_dig_wr2),1,R_s); %Moment arm, subtalar, x' axis
% flex_dig_maMZ(i,j) = CrossProd(flex_dig_wr5, (VecTrans(T_m\eye(4),flex_dig_wr4)-flex_dig_wr5)/norm(VecTrans(T_m\eye(4),flex_dig_wr4)-flex_dig_wr5),3,R_m); %Moment arm, MTP, z' axis

%Flexor Hallucis Longus
% flex_hal_length(i,j) = norm(t4-flex_hal_wr1)+norm(VecTrans(T_a*T_s,flex_hal_wr2)-flex_hal_wr1)+norm(flex_hal_wr2-flex_hal_wr3)+norm(flex_hal_wr3-flex_hal_wr4)+norm(VecTrans(T_m,flex_hal_wr5)-flex_hal_wr4)+norm(flex_hal_wr5-a8); %Total length
% flex_hal_maAZ(i,j) = CrossProd(VecTrans(T_s,flex_hal_wr2), (VecTrans(T_a\eye(4),flex_hal_wr1)-VecTrans(T_s,flex_hal_wr2))/norm(VecTrans(T_a\eye(4),flex_hal_wr1)-VecTrans(T_s,flex_hal_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% flex_hal_maSX(i,j) = CrossProd(flex_hal_wr2, (VecTrans((T_a*T_s)\eye(4),flex_hal_wr1)-flex_hal_wr2)/norm(VecTrans((T_a*T_s)\eye(4),flex_hal_wr1)-flex_hal_wr2),1,R_s); %Moment arm, subtalar, x' axis
% flex_hal_maMZ(i,j) = CrossProd(flex_hal_wr5, (VecTrans(T_m\eye(4),flex_hal_wr4)-flex_hal_wr5)/norm(VecTrans(T_m\eye(4),flex_hal_wr4)-flex_hal_wr5),3,R_m); %Moment arm, MTP, z' axis

%Tibialis Anterior
% tib_ant_length(i,j) = norm(tib_ant_wr1-t2)+norm(VecTrans(T_a*T_s,a3)-tib_ant_wr1); %Total length
% tib_ant_maAZ(i,j) = CrossProd(VecTrans(T_s,a3), (VecTrans(T_a\eye(4),tib_ant_wr1)-VecTrans(T_s,a3))/norm(VecTrans(T_a\eye(4),tib_ant_wr1)-VecTrans(T_s,a3)),3,R_a); %Moment arm, Ankle, z' axis
% tib_ant_maSX(i,j) = CrossProd(a3, (VecTrans((T_a*T_s)\eye(4),tib_ant_wr1)-a3)/norm(VecTrans((T_a*T_s)\eye(4),tib_ant_wr1)-a3),1,R_s); %Moment arm, subtalar, x' axis

%Peroneus Brevis
% per_brev_wr2 = [-0.014; -0.429; 0.029]; %Tibia, right, wrapping point 2
% per_brev_wr3 = [0.047; 0.027; 0.023]; %Calcn, right, wrapping point 3
% per_brev_length(i,j) = norm(t2-per_brev_wr1)+norm(per_brev_wr1-per_brev_wr2)+norm(VecTrans(T_a*T_s,per_brev_wr3)-per_brev_wr2)+norm(per_brev_wr3-a2); %Total length
% per_brev_maAZ(i,j) = CrossProd(VecTrans(T_s,per_brev_wr3), (VecTrans(T_a\eye(4),per_brev_wr2)-VecTrans(T_s,per_brev_wr3))/norm(VecTrans(T_a\eye(4),per_brev_wr2)-VecTrans(T_s,per_brev_wr3)),3,R_a); %Moment arm, Ankle, z' axis
% per_brev_maSX(i,j) = CrossProd(per_brev_wr3, (VecTrans((T_a*T_s)\eye(4),per_brev_wr2)-per_brev_wr3)/norm(VecTrans((T_a*T_s)\eye(4),per_brev_wr2)-per_brev_wr3),1,R_s); %Moment arm, subtalar, x' axis

%Peroneus Longus
% per_long_length(i,j) = norm(t2-per_long_wr1)+norm(per_long_wr1-per_long_wr2)+norm(VecTrans(T_a*T_s,per_long_wr3)-per_long_wr2)+norm(per_long_wr3-per_long_wr4)+norm(per_long_wr5-per_long_wr4)+norm(per_long_wr5-a3); %Total length
% per_long_maAZ(i,j) = CrossProd(VecTrans(T_s,per_long_wr3), (VecTrans(T_a\eye(4),per_long_wr2)-VecTrans(T_s,per_long_wr3))/norm(VecTrans(T_a\eye(4),per_long_wr2)-VecTrans(T_s,per_long_wr3)),3,R_a); %Moment arm, Ankle, z' axis
% per_long_maSX(i,j) = CrossProd(per_long_wr3, (VecTrans((T_a*T_s)\eye(4),per_long_wr2)-per_long_wr3)/norm(VecTrans((T_a*T_s)\eye(4),per_long_wr2)-per_long_wr3),1,R_s); %Moment arm, subtalar, x' axis

%Peroneus Tertius
% per_tert_length(i,j) = norm(per_tert_wr1-t2)+norm(VecTrans(T_a*T_s,a2)-per_tert_wr1) ; %Total length
% per_tert_maAZ(i,j) = CrossProd(VecTrans(T_s,a2), (VecTrans(T_a\eye(4),per_tert_wr1)-VecTrans(T_s,a2))/norm(VecTrans(T_a\eye(4),per_tert_wr1)-VecTrans(T_s,a2)),3,R_a); %Moment arm, Ankle, z' axis
% per_tert_maSX(i,j) = CrossProd(a2, (VecTrans((T_a*T_s)\eye(4),per_tert_wr1)-a2)/norm(VecTrans((T_a*T_s)\eye(4),per_tert_wr1)-a2),1,R_s); %Moment arm, subtalar, x' axis

%Extensor Digitorum Longus
% ext_dig_length(i,j) = norm(t2-ext_dig_wr1)+norm(VecTrans(T_a*T_s,ext_dig_wr2)-ext_dig_wr1)+norm(ext_dig_wr2-ext_dig_wr3)+norm(VecTrans(T_m,ext_dig_wr4)-ext_dig_wr3)+norm(ext_dig_wr4-a5); %Total length
% ext_dig_maAZ(i,j) = CrossProd(VecTrans(T_s,ext_dig_wr2), (VecTrans(T_a\eye(4),ext_dig_wr1)-VecTrans(T_s,ext_dig_wr2))/norm(VecTrans(T_a\eye(4),ext_dig_wr1)-VecTrans(T_s,ext_dig_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% ext_dig_maSX(i,j) = CrossProd(ext_dig_wr2, (VecTrans((T_a*T_s)\eye(4),ext_dig_wr1)-ext_dig_wr2)/norm(VecTrans((T_a*T_s)\eye(4),ext_dig_wr1)-ext_dig_wr2),1,R_s); %Moment arm, subtalar, x' axis
% ext_dig_maMZ(i,j) = CrossProd(ext_dig_wr4, (VecTrans(T_m\eye(4),ext_dig_wr3)-ext_dig_wr4)/norm(VecTrans(T_m\eye(4),ext_dig_wr3)-ext_dig_wr4),3,R_m); %Moment arm, MTP, z' axis

%Extensor Hallucis Longus
% ext_hal_length(i,j) = norm(t2-ext_hal_wr1)+norm(VecTrans(T_a*T_s,ext_hal_wr2)-ext_hal_wr1)+norm(ext_hal_wr2-ext_hal_wr3)+norm(ext_hal_wr3-ext_hal_wr4)+norm(VecTrans(T_m,ext_hal_wr5)-ext_hal_wr4)+norm(ext_hal_wr5-a6); %Total length
% ext_hal_maAZ(i,j) = CrossProd(VecTrans(T_s,ext_hal_wr2), (VecTrans(T_a\eye(4),ext_hal_wr1)-VecTrans(T_s,ext_hal_wr2))/norm(VecTrans(T_a\eye(4),ext_hal_wr1)-VecTrans(T_s,ext_hal_wr2)),3,R_a); %Moment arm, Ankle, z' axis
% ext_hal_maSX(i,j) = CrossProd(ext_hal_wr2, (VecTrans((T_a*T_s)\eye(4),ext_hal_wr1)-ext_hal_wr2)/norm(VecTrans((T_a*T_s)\eye(4),ext_hal_wr1)-ext_hal_wr2),1,R_s); %Moment arm, subtalar, x' axis
% ext_hal_maMZ(i,j) = CrossProd(ext_hal_wr5, (VecTrans(T_m\eye(4),ext_hal_wr4)-ext_hal_wr5)/norm(VecTrans(T_m\eye(4),ext_hal_wr4)-ext_hal_wr5),3,R_m); %Moment arm, MTP, z' axis

%Erector Spinae
% ercspn_length(i,j) = norm(VecTrans(T_b,b2) - ercspn_wr1)+norm(ercspn_wr1 - p7); %Total length, right
% ercspn_maX(i,j) = CrossProd(b2, (VecTrans(T_b\eye(4),ercspn_wr1)-b2)/norm(VecTrans(T_b\eye(4),ercspn_wr1)-b2),1); %Moment arm, x axis
% ercspn_maY(i,j) = CrossProd(b2, (VecTrans(T_b\eye(4),ercspn_wr1)-b2)/norm(VecTrans(T_b\eye(4),ercspn_wr1)-b2),2); %Moment arm, y axis
% ercspn_maZ(i,j) = CrossProd(b2, (VecTrans(T_b\eye(4),ercspn_wr1)-b2)/norm(VecTrans(T_b\eye(4),ercspn_wr1)-b2),3); %Moment arm, z axis

%Internal Oblique
% intobl_length(i,j) = norm(VecTrans(T_b,b3) - p2); %Total length, right
% intobl_maX(i,j) = CrossProd(b3, (VecTrans(T_b\eye(4),p2)-b3)/norm(VecTrans(T_b\eye(4),p2)-b3),1); %Moment arm, x axis
% intobl_maY(i,j) = CrossProd(b3, (VecTrans(T_b\eye(4),p2)-b3)/norm(VecTrans(T_b\eye(4),p2)-b3),2); %Moment arm, y axis
% intobl_maZ(i,j) = CrossProd(b3, (VecTrans(T_b\eye(4),p2)-b3)/norm(VecTrans(T_b\eye(4),p2)-b3),3); %Moment arm, z axis

%External Oblique extobl
% extobl_length(i,j) = norm(VecTrans(T_b,b4) - p4); %Total length, right
% extobl_maX(i,j) = CrossProd(b4, (VecTrans(T_b\eye(4),p4)-b4)/norm(VecTrans(T_b\eye(4),p4)-b4),1); %Moment arm, x axis
% extobl_maY(i,j) = CrossProd(b4, (VecTrans(T_b\eye(4),p4)-b4)/norm(VecTrans(T_b\eye(4),p4)-b4),2); %Moment arm, y axis
% extobl_maZ(i,j) = CrossProd(b4, (VecTrans(T_b\eye(4),p4)-b4)/norm(VecTrans(T_b\eye(4),p4)-b4),3); %Moment arm, z axis
end
end
%%
%%Size function calculations
% [glut_max_dia, glut_max_rest, glut_max_long]= Size(glut_max_length, glut_max_mif)
% [semimem_dia, semimem_rest, semimem_long]= Size(semimem_length, semimem_mif)
% [bifemlh_dia, bifemlh_rest, bifemlh_long]= Size(bifemlh_length, bifemlh_mif)
% [bifemsh_dia, bifemsh_rest, bifemsh_long]= Size(bifemsh_length, bifemsh_mif)
% [sar_dia, sar_rest, sar_long]= Size(sar_length, sar_mif)
% [add_mag_dia, add_mag_rest, add_mag_long]= Size(add_mag_length, add_mag_mif)
% [tfl_dia, tfl_rest, tfl_long]= Size(tfl_length, tfl_mif)
% [grac_dia, grac_rest, grac_long]= Size(grac_length, grac_mif)
% [iliacus_dia, iliacus_rest, iliacus_long]= Size(iliacus_length, iliacus_mif)
% [psoas_dia, psoas_rest, psoas_long]= Size(psoas_length, psoas_mif)
% [rect_fem_dia, rect_fem_rest, rect_fem_long]= Size(rect_fem_length, rect_fem_mif)
% [vas_int_dia, vas_int_rest, vas_int_long]= Size(vas_int_length, vas_int_mif)
% [med_gas_dia, med_gas_rest, med_gas_long]= Size(med_gas_length, med_gas_mif)
% [lat_gas_dia, lat_gas_rest, lat_gas_long]= Size(lat_gas_length, lat_gas_mif)
% [soleus_dia, soleus_rest, soleus_long]= Size(soleus_length, soleus_mif)
% [tib_post_dia, tib_post_rest, tib_post_long]= Size(tib_post_length, tib_post_mif)
% [flex_dig_dia, flex_dig_rest, flex_dig_long]= Size(flex_dig_length, flex_dig_mif)
% [flex_hal_dia, flex_hal_rest, flex_hal_long]= Size(flex_hal_length, flex_hal_mif)
% [tib_ant_dia, tib_ant_rest, tib_ant_long]= Size(tib_ant_length, tib_ant_mif)
% [per_brev_dia, per_brev_rest, per_brev_long]= Size(per_brev_length, per_brev_mif)
% [per_long_dia, per_long_rest, per_long_long]= Size(per_long_length, per_long_mif)
% [per_tert_dia, per_tert_rest, per_tert_long]= Size(per_tert_length, per_tert_mif)
% [ext_dig_dia, ext_dig_rest, ext_dig_long]= Size(ext_dig_length, ext_dig_mif)
% [ext_hal_dia, ext_hal_rest, ext_hal_long]= Size(ext_hal_length, ext_hal_mif)
% [ercspn_dia, ercspn_rest, ercspn_long]= Size(ercspn_length, ercspn_mif)
% [intobl_dia, intobl_rest, intobl_long]= Size(intobl_length, intobl_mif)
% [extobl_dia, extobl_rest, extobl_long]= Size(extobl_length, extobl_mif)
%%
for j=1:t
 for i = 1:s

%%Muscle force, and torque
%Gluteus Maximus
% glut_max_f(i,j) = festo(glut_max_length(i,j), glut_max_rest, glut_max_dia, glut_max_long); %Force
% glut_max_tqX(i,j) = glut_max_f(i,j)*glut_max_maX(i,j); %Torque, x axis
% glut_max_tqY(i,j) = glut_max_f(i,j)*glut_max_maY(i,j); %Torque, y axis
% glut_max_tqZ(i,j) = glut_max_f(i,j)*glut_max_maZ(i,j); %Torque, z axis

%Semimembranosus
% semimem_f(i,j) = festo(semimem_length(i,j), semimem_rest, semimem_dia, semimem_long); %Force
% semimem_tqHX(i,j) = semimem_f(i,j)*semimem_maHX(i,j); %Torque, Hip, x axis
% semimem_tqHY(i,j) = semimem_f(i,j)*semimem_maHY(i,j); %Torque, Hip, y axis
% semimem_tqHZ(i,j) = semimem_f(i,j)*semimem_maHZ(i,j); %Torque, Hip, z axis
% semimem_tqKZ(i,j) = semimem_f(i,j)*semimem_maKZ(i,j); %Torque, Knee, z axis

%Bicep Femoris, Long Head
% bifemlh_f(i,j) = festo(bifemlh_length(i,j), bifemlh_rest, bifemlh_dia, bifemlh_long); %Force
% bifemlh_tqHX(i,j) = bifemlh_f(i,j)*bifemlh_maHX(i,j); %Torque, Hip, x axis
% bifemlh_tqHY(i,j) = bifemlh_f(i,j)*bifemlh_maHY(i,j); %Torque, Hip, y axis
% bifemlh_tqHZ(i,j) = bifemlh_f(i,j)*bifemlh_maHZ(i,j); %Torque, Hip, z axis
% bifemlh_tqKZ(i,j) = bifemlh_f(i,j)*bifemlh_maKZ(i,j); %Torque, Knee, z axis

%Bicep Femoris, Short Head
bifemsh_f(i,j) = festo(bifemsh_length(i,j), bifemsh_rest, bifemsh_dia, bifemsh_long); %Force
bifemsh_tqZ(i,j) = bifemsh_f(i,j)*bifemsh_maZ(i,j); %Torque, Knee, z axis

%Sartorius
% sar_f(i,j) = festo(sar_length(i,j), sar_rest, sar_dia, sar_long); %Force
% sar_tqHX(i,j) = sar_f(i,j)*sar_maHX(i,j); %Torque, Hip, x axis
% sar_tqHY(i,j) = sar_f(i,j)*sar_maHY(i,j); %Torque, Hip, y axis
% sar_tqHZ(i,j) = sar_f(i,j)*sar_maHZ(i,j); %Torque, Hip, z axis
% sar_tqKZ(i,j) = sar_f(i,j)*sar_maKZ(i,j); %Torque, Knee, z axis

%Adductor Magnus
% add_mag_length(i,j) = norm(p9 - add_mag_wr1) + norm(VecTrans(T_h,add_mag_wr2) - add_mag_wr1)+norm(add_mag_wr2-f4); %Total length
% add_mag_maX(i,j) = CrossProd(add_mag_wr2, (VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2)/norm(VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2),1); %Moment arm, x axis
% add_mag_maY(i,j) = CrossProd(add_mag_wr2, (VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2)/norm(VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2),2); %Moment arm, y axis
% add_mag_maZ(i,j) = CrossProd(add_mag_wr2, (VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2)/norm(VecTrans(T_h\eye(4),add_mag_wr1)-add_mag_wr2),3); %Moment arm, z axis

%Tensor Fasciae Latae
% tfl_length(i,j) = norm(VecTrans(T_h,tfl_wr1)-p3)+norm(tfl_wr2 - tfl_wr1)+norm(tfl_wr2-VecTrans(T_k,t2))+norm(t2-t5); %Total length, right
% tfl_maHX(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),p3)-tfl_wr1)/norm(VecTrans(T_h\eye(4),p3)-tfl_wr1),1); %Moment arm, Hip, x axis
% tfl_maHY(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),p3)-tfl_wr1)/norm(VecTrans(T_h\eye(4),p3)-tfl_wr1),2); %Moment arm, Hip, y axis
% tfl_maHZ(i,j) = CrossProd(tfl_wr1, (VecTrans(T_h\eye(4),p3)-tfl_wr1)/norm(VecTrans(T_h\eye(4),p3)-tfl_wr1),3); %Moment arm, Hip, z axis
% tfl_maKZ(i,j) = CrossProd(t2, (VecTrans(T_k\eye(4),tfl_wr2)-t2)/norm(VecTrans(T_k\eye(4),tfl_wr2)-t2),3); %Moment arm, Knee, z axis

%Gracilis
% grac_f(i,j) = festo(grac_length(i,j), grac_rest, grac_dia, grac_long); %Force
% grac_tqHX(i,j) = grac_f(i,j)*grac_maHX(i,j); %Torque, Hip, x axis
% grac_tqHY(i,j) = grac_f(i,j)*grac_maHY(i,j); %Torque, Hip, y axis
% grac_tqHZ(i,j) = grac_f(i,j)*grac_maHZ(i,j); %Torque, Hip, z axis
% grac_tqKZ(i,j) = grac_f(i,j)*grac_maKZ(i,j); %Torque, Knee, z axis

%Iliacus
% iliacus_f(i,j) = festo(iliacus_length(i,j), iliacus_rest, iliacus_dia, iliacus_long); %Force
% iliacus_tqX(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, x axis
% iliacus_tqY(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, y axis
% iliacus_tqZ(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, z axis

%Psoas
% psoas_f(i,j) = festo(psoas_length(i,j), psoas_rest, psoas_dia, psoas_long); %Force
% psoas_tqHX(i,j) = psoas_f(i,j)*psoas_maHX(i,j); %Torque, Hip, x axis
% psoas_tqHY(i,j) = psoas_f(i,j)*psoas_maHY(i,j); %Torque, Hip, y axis
% psoas_tqHZ(i,j) = psoas_f(i,j)*psoas_maHZ(i,j); %Torque, Hip, z axis
% psoas_tqBX(i,j) = psoas_f(i,j)*psoas_maBX(i,j); %Torque, Back, x axis
% psoas_tqBY(i,j) = psoas_f(i,j)*psoas_maBY(i,j); %Torque, Back, y axis
% psoas_tqBZ(i,j) = psoas_f(i,j)*psoas_maBZ(i,j); %Torque, Back, z axis

%Rectus Femoris
% rect_fem_f(i,j) = festo(rect_fem_length(i,j), rect_fem_rest, rect_fem_dia, rect_fem_long); %Force
% rect_fem_tqHX(i,j) = rect_fem_f(i,j)*rect_fem_maHX(i,j); %Torque, Hip, x axis
% rect_fem_tqHY(i,j) = rect_fem_f(i,j)*rect_fem_maHY(i,j); %Torque, Hip, y axis
% rect_fem_tqHZ(i,j) = rect_fem_f(i,j)*rect_fem_maHZ(i,j); %Torque, Hip, z axis
% rect_fem_tqKZ(i,j) = rect_fem_f(i,j)*rect_fem_maKZ(i,j); %Torque, Knee, z axis

%Vastus Intermedius
vas_int_f(i,j) = festo(vas_int_length(i,j), vas_int_rest, vas_int_dia, vas_int_long); %Force
vas_int_tqZ(i,j) = vas_int_f(i,j)*vas_int_maZ(i,j); %Torque, z axis

%Medial Gastroocnemius
% med_gas_f(i,j) = festo(med_gas_length(i,j), med_gas_rest, med_gas_dia, med_gas_long); %Force
% med_gas_tqKZ(i,j) = med_gas_f(i,j)*med_gas_maKZ(i,j); %Torque, Knee, z axis
% med_gas_tqAZ(i,j) = med_gas_f(i,j)*med_gas_maAZ(i,j); %Torque, Ankle, z' axis
% med_gas_tqSX(i,j) = med_gas_f(i,j)*med_gas_maSX(i,j); %Torque, Subtalar, x' axis

%Lateral Gastrocenemius
% lat_gas_f(i,j) = festo(lat_gas_length(i,j), lat_gas_rest, lat_gas_dia, lat_gas_long); %Force
% lat_gas_tqKZ(i,j) = lat_gas_f(i,j)*lat_gas_maKZ(i,j); %Torque, Knee, z axis
% lat_gas_tqAZ(i,j) = lat_gas_f(i,j)*lat_gas_maAZ(i,j); %Torque, Ankle, z' axis
% lat_gas_tqSX(i,j) = lat_gas_f(i,j)*lat_gas_maSX(i,j); %Torque, Subtalar, x' axis

%Soleus
% soleus_f(i,j) = festo(soleus_length(i,j), soleus_rest, soleus_dia, soleus_long); %Force
% soleus_tqAZ(i,j) = soleus_f(i,j)*soleus_maAZ(i,j); %Torque, Ankle, z' axis
% soleus_tqSX(i,j) = soleus_f(i,j)*soleus_maSX(i,j); %Torque, Subtalar, x' axis

%Tibialis Posterior
% tib_post_f(i,j) = festo(tib_post_length(i,j), tib_post_rest, tib_post_dia, tib_post_long); %Force
% tib_post_tqAZ(i,j) = tib_post_f(i,j)*tib_post_maAZ(i,j); %Torque, Ankle, z' axis
% tib_post_tqSX(i,j) = tib_post_f(i,j)*tib_post_maSX(i,j); %Torque, Subtalar, x' axis

%Flexor Digitorum Longus
% flex_dig_f(i,j) = festo(flex_dig_length(i,j), flex_dig_rest, flex_dig_dia, flex_dig_long); %Force
% flex_dig_tqAZ(i,j) = flex_dig_f(i,j)*flex_dig_maAZ(i,j); %Torque, Ankle, z' axis
% flex_dig_tqSX(i,j) = flex_dig_f(i,j)*flex_dig_maSX(i,j); %Torque, subtalar, x' axis
% flex_dig_tqMZ(i,j) = flex_dig_f(i,j)*flex_dig_maMZ(i,j); %Torque, MTP, z' axis

%Flexor Hallucis Longus
% flex_hal_f(i,j) = festo(flex_hal_length(i,j), flex_hal_rest, flex_hal_dia, flex_hal_long); %Force
% flex_hal_tqAZ(i,j) = flex_hal_f(i,j)*flex_hal_maAZ(i,j); %Torque, Ankle, z' axis
% flex_hal_tqSX(i,j) = flex_hal_f(i,j)*flex_hal_maSX(i,j); %Torque, subtalar, x' axis
% flex_hal_tqMZ(i,j) = flex_hal_f(i,j)*flex_hal_maMZ(i,j); %Torque, MTP, z' axis

%Tibialis Anterior
% tib_ant_f(i,j) = festo(tib_ant_length(i,j), tib_ant_rest, tib_ant_dia, tib_ant_long); %Force
% tib_ant_tqAZ(i,j) = tib_ant_f(i,j)*tib_ant_maAZ(i,j); %Torque, Ankle, z' axis
% tib_ant_tqSX(i,j) = tib_ant_f(i,j)*tib_ant_maSX(i,j); %Torque, Subtalar, x' axis

%Peroneus Brevis
% per_brev_f(i,j) = festo(per_brev_length(i,j), per_brev_rest, per_brev_dia, per_brev_long); %Force
% per_brev_tqAZ(i,j) = per_brev_f(i,j)*per_brev_maAZ(i,j); %Torque, Ankle, z' axis
% per_brev_tqSX(i,j) = per_brev_f(i,j)*per_brev_maSX(i,j); %Torque, Subtalar, x' axis

%Peroneus Longus
% per_long_f(i,j) = festo(per_long_length(i,j), per_long_rest, per_long_dia, per_long_long); %Force
% per_long_tqAZ(i,j) = per_long_f(i,j)*per_long_maAZ(i,j); %Torque, Ankle, z' axis
% per_long_tqSX(i,j) = per_long_f(i,j)*per_long_maSX(i,j); %Torque, Subtalar, x' axis

%Peroneus Tertius
% per_tert_f(i,j) = festo(per_tert_length(i,j), per_tert_rest, per_tert_dia, per_tert_long); %Force
% per_tert_tqAZ(i,j) = per_tert_f(i,j)*per_tert_maAZ(i,j); %Torque, Ankle, z' axis
% per_tert_tqSX(i,j) = per_tert_f(i,j)*per_tert_maSX(i,j); %Torque, Subtalar, x' axis

%Extensor Digitorum Longus
% ext_dig_f(i,j) = festo(ext_dig_length(i,j), ext_dig_rest, ext_dig_dia, ext_dig_long); %Force
% ext_dig_tqAZ(i,j) = ext_dig_f(i,j)*ext_dig_maAZ(i,j); %Torque, Ankle, z' axis
% ext_dig_tqSX(i,j) = ext_dig_f(i,j)*ext_dig_maSX(i,j); %Torque, subtalar, x' axis
% ext_dig_tqMZ(i,j) = ext_dig_f(i,j)*ext_dig_maMZ(i,j); %Torque, MTP, z' axis

%Extensor Hallucis Longus
% ext_hal_f(i,j) = festo(ext_hal_length(i,j), ext_hal_rest, ext_hal_dia, ext_hal_long); %Force
% ext_hal_tqAZ(i,j) = ext_hal_f(i,j)*ext_hal_maAZ(i,j); %Torque, Ankle, z' axis
% ext_hal_tqSX(i,j) = ext_hal_f(i,j)*ext_hal_maSX(i,j); %Torque, subtalar, x' axis
% ext_hal_tqMZ(i,j) = ext_hal_f(i,j)*ext_hal_maMZ(i,j); %Torque, MTP, z' axis

%Erector Spinae
% ercspn_f(i,j) = festo(ercspn_length(i,j), ercspn_rest, ercspn_dia, ercspn_long); %Force
% ercspn_tqX(i,j) = ercspn_f(i,j)*ercspn_maX(i,j); %Torque, x axis
% ercspn_tqY(i,j) = ercspn_f(i,j)*ercspn_maY(i,j); %Torque, y axis
% ercspn_tqZ(i,j) = ercspn_f(i,j)*ercspn_maZ(i,j); %Torque, z axis

%Internal Oblique
% intobl_f(i,j) = festo(intobl_length(i,j), intobl_rest, intobl_dia, intobl_long); %Force
% intobl_tqX(i,j) = intobl_f(i,j)*intobl_maX(i,j); %Torque, x axis
% intobl_tqY(i,j) = intobl_f(i,j)*intobl_maY(i,j); %Torque, y axis
% intobl_tqZ(i,j) = intobl_f(i,j)*intobl_maZ(i,j); %Torque, z axis

%External Oblique extobl
% extobl_f(i,j) = festo(extobl_length(i,j), extobl_rest, extobl_dia, extobl_long); %Force
% extobl_tqX(i,j) = extobl_f(i,j)*extobl_maX(i,j); %Torque, x axis
% extobl_tqY(i,j) = extobl_f(i,j)*extobl_maY(i,j); %Torque, y axis
% extobl_tqZ(i,j) = extobl_f(i,j)*extobl_maZ(i,j); %Torque, z axis
end
end
%%
%%Plots and Surface Figures
% figure
% plot(theta_K, rect_fem_tqKZ)
% title('rect_fem torque for knee flexion')
% ylabel('Torque, N-m')
% xlabel('Extension, Degrees')

% uni_hip_tqZ = glut_max_tqZ+add_mag_tqZ+iliacus_tqZ+psoas_tqHZ;
% uni_hip_tqY = glut_max_tqY+add_mag_tqY+iliacus_tqY+psoas_tqHY;
% uni_hip_tqX = glut_max_tqX+add_mag_tqX+iliacus_tqX+psoas_tqHX;

% bi_hip_tqHX = semimem_tqHX+bifemlh_tqHX+sar_tqHX+tfl_tqHX+rect_fem_tqHX+grac_tqHX;
% bi_hip_tqHY = semimem_tqHY+bifemlh_tqHY+sar_tqHY+tfl_tqHY+rect_fem_tqHY+grac_tqHY;
% bi_hip_tqHZ = semimem_tqHZ+bifemlh_tqHZ+sar_tqHZ+tfl_tqHZ+rect_fem_tqHZ+grac_tqHZ;
% bi_knee_tqKZ = semimem_tqKZ+bifemlh_tqKZ+sar_tqKZ+tfl_tqKZ+rect_fem_tqKZ+grac_tqKZ;
% uni_knee_tqZ = bifemsh_tqZ+vas_int_tqZ;

% calves_tqKZ = med_gas_tqKZ + lat_gas_tqKZ;
% calves_tqAZ = med_gas_tqAZ + lat_gas_tqAZ + soleus_tqAZ;
% calves_tqSX = med_gas_tqSX + lat_gas_tqSX + soleus_tqSX;

% toes_tqAZ = flex_dig_tqAZ + flex_hal_tqAZ +ext_dig_tqAZ + ext_hal_tqAZ;
% toes_tqSX = flex_dig_tqSX + flex_hal_tqSX +ext_dig_tqSX + ext_hal_tqSX;
% toes_tqMZ = flex_dig_tqMZ + flex_hal_tqMZ +ext_dig_tqMZ + ext_hal_tqMZ;

% foot_tqAZ = tib_post_tqAZ + tib_ant_tqAZ + per_brev_tqAZ + per_long_tqAZ + per_tert_tqAZ;
% foot_tqSX = tib_post_tqSX + tib_ant_tqSX + per_brev_tqSX + per_long_tqSX + per_tert_tqSX;

figure
surf(theta_hZ, theta_K, bifemsh_tqZ, 'EdgeColor','none')
title('Biceps femoris length')
ylabel('Knee Flexion, Degrees')
xlabel('Hip Flexion, Degrees')
zlabel('Length, m')
colorbar
colormap(GetBichromaticColormap())

figure
surf(theta_hZ, theta_K, vas_int_length, 'EdgeColor','none')
title('Vastus Int length')
ylabel('Knee Flexion, Degrees')
xlabel('Hip Flexion, Degrees')
zlabel('Length, m')
colorbar
colormap(GetBichromaticColormap())

% figure
% surf(theta_hZ, theta_K, bifemlh_tqHZ, 'EdgeColor','none')
% title('Hip Torque, bifemlh, Z axis')
% ylabel('Knee Flexion, Degrees')
% xlabel('Hip Flexion, Degrees')
% zlabel('Torque, N*m')
% colorbar
% colormap(GetBichromaticColormap())
% 
% figure
% surf(theta_hZ, theta_K, bifemlh_tqKZ, 'EdgeColor','none')
% title('Knee Torque, bifemlh, Z axis')
% ylabel('Knee Flexion, Degrees')
% xlabel('Hip Flexion, Degrees')
% zlabel('Torque, N*m')
% colorbar
% colormap(GetBichromaticColormap())