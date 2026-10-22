clear; clc; close all;
% Robotic Legs Kinematics
%Units in m, N, rad/s unless otherwise stated.

load AttachPoints.mat  %Attachment points for muscles with stationary attachment points
load MovePoints.mat %Moving attachment points and knee data
load MuscProp.mat      %Muscle specific properties
IndexPrep %Matrices for length, force, and moment arm  

%Range of motion, radians
range_hz = [ -20*c, 140*c];    %Hip (Wikipedia)
thetamin_hz = range_hz(1);      %20 deg extention
thetamax_hz = range_hz(2);      %140 deg flexion
range_k = [-2.0943951, 0.17453293];    %Knee
thetamin_k = range_k(1);       %120 deg flexion
thetamax_k = range_k(2);       %10 deg extension
thetamin_a = -60*c;       %ankle (OpenSim values)
thetamax_a = 60*c;
thetamin_b = -90*c;       %back  (OpenSim values)
thetamax_b = 90*c;
k = 180/pi;     %convert radians to degrees
c = pi/180;     %convert degrees to radians

%Indexing prep, see IndexPrep.m
s = 100;
t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)

for j=1:t
 for i = 1:s
%Degree of freedom, position or range
theta_h_z = thetamin_hz + (thetamax_hz-thetamin_hz)*(j-1)/s; %Hip, right, extension-/flexion+
thetamin_hy = c*(-32.5 - (theta_h_z*k/160)*20); %Max External Rotation (linear interpolation of Wikipedia RoM)
thetamax_hy = c*(40);  %Max Internal Rotation (Wikipedia)
theta_h_y = 0; %Hip, right, external-/internal+ rotation
thetamin_hx = c*(-53.75 - (theta_h_z*k/160)*30); %Max Abduction (linear interpolation of Wikipedia RoM)
thetamax_hx = c*(28.75 - (theta_h_z*k/160)*10); %Max Adduction (linear interpolation of Wikipedia RoM)
theta_h_x = thetamin_hx + (thetamax_hx-thetamin_hx)*(i-1)/s; %Hip, right, adduction-/abduction+
theta_k = 0; %Knee, right, flexion-/extension+ 
theta_b_x = 0; %Lumbar bending left-/right+
theta_b_y = 0; %Lumbar rotation cw-/ccw+
theta_b_z = 0; %Back flexion-/extension+
theta_a = 0; %Ankle, right, extension-/flexion+
theta_s = 0; %Subtalar, right, external-/internal+ rotation
theta_m = 0; %Metatarsophalangeal (mtp), right, flexion-/extension+

%Index for plotting
theta_hX(i,j) = theta_h_x*k; %Convert from radians to degrees
%theta_hY(i,1) = (thetamin_hy + (thetamax_hy-thetamin_hy)*(i-1)/s)*k;
theta_hZ(j,1) = theta_h_z*k;
%theta_K(i,1) = (thetamin_k + (thetamax_k-thetamin_k)*(i-1)/s)*k;
%theta_bX(i,1) = (thetamin_b + (thetamax_b-thetamin_b)*(i-1)/s)*k;
%theta_bY(i,1) = (thetamin_b + (thetamax_b-thetamin_b)*(i-1)/s)*k;
%theta_bZ(i,1) = (thetamin_b + (thetamax_b-thetamin_b)*(i-1)/s)*k;
%theta_A(i,1) = (thetamin_a + (thetamax_a-thetamin_a)*(i-1)/s)*k;
%theta_S(i,1) = (thetamin_b + (thetamax_b-thetamin_b)*(i-1)/s)*k;
%theta_M(i,1) = (thetamin_a + (thetamax_a-thetamin_a)*(i-1)/s)*k;

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
omg_a = [-.10501355; -0.17402245; 097912632]; %Axis of rotation
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
p_k = [f1(theta_k); f2(theta_k); 0]; %Knee, function of knee angle (see MovePoints)
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
glut_med1_length(i,j) = norm(VecTrans(T_h,glut_med1_i) - glut_med1_o); %Total length
glut_med1_ma(i,j) = norm(cross(glut_med1_o - p_h, (VecTrans(T_h,glut_med1_i)-glut_med1_o)/norm((VecTrans(T_h,glut_med1_i)-glut_med1_o)))); %Moment arm
% glut_med1_f(i,j) = forz(glut_med1_length(i,j), glut_med1_mif, glut_med1_ofl, glut_med1_tsl, glut_med1_pa); %Force
% glut_med1_torque(i,j) = glut_med1_f(i,j)*glut_med1_ma(i,j); %Torque
%Gluteus Medius, 2
glut_med2_length(i,j) = norm(VecTrans(T_h,glut_med2_i) - glut_med2_o); %Total length
glut_med2_ma(i,j) = norm(cross(glut_med2_o - p_h, (VecTrans(T_h,glut_med2_i)-glut_med2_o)/norm((VecTrans(T_h,glut_med2_i)-glut_med2_o)))); %Moment arm
% glut_med2_f(i,j) = forz(glut_med2_length(i,j), glut_med2_mif, glut_med2_ofl, glut_med2_tsl, glut_med2_pa); %Force
% glut_med2_torque(i,j) = glut_med2_f(i,j)*glut_med2_ma(i,j); %Torque
%Gluteus Medius, 3
glut_med3_length(i,j) = norm(VecTrans(T_h,glut_med3_i) - glut_med3_o); %Total length
glut_med3_ma(i,j) = norm(cross(glut_med3_o - p_h, (VecTrans(T_h,glut_med3_i)-glut_med3_o)/norm((VecTrans(T_h,glut_med3_i)-glut_med3_o)))); %Moment arm
% glut_med3_f(i,j) = forz(glut_med3_length(i,j), glut_med3_mif, glut_med3_ofl, glut_med3_tsl, glut_med3_pa); %Force
% glut_med3_torque(i,j) = glut_med3_f(i,j)*glut_med3_ma(i,j); %Torque
%Gluteus Minimus, 1
glut_min1_length(i,j) = norm(VecTrans(T_h,glut_min1_i) - glut_min1_o); %Total length
glut_min1_ma(i,j) = norm(cross(glut_min1_o - p_h, (VecTrans(T_h,glut_min1_i)-glut_min1_o)/norm((VecTrans(T_h,glut_min1_i)-glut_min1_o)))); %Moment arm
% glut_min1_f(i,j) = forz(glut_min1_length(i,j), glut_min1_mif, glut_min1_ofl, glut_min1_tsl, glut_min1_pa); %Force
% glut_min1_torque(i,j) = glut_min1_f(i,j)*glut_min1_ma(i,j); %Torque
%Gluteus Minimus, 2
glut_min2_length(i,j) = norm(VecTrans(T_h,glut_min2_i) - glut_min2_o); %Total length
glut_min2_ma(i,j) = norm(cross(glut_min2_o - p_h, (VecTrans(T_h,glut_min2_i)-glut_min2_o)/norm((VecTrans(T_h,glut_min2_i)-glut_min2_o)))); %Moment arm
% glut_min2_f(i,j) = forz(glut_min2_length(i,j), glut_min2_mif, glut_min2_ofl, glut_min2_tsl, glut_min2_pa); %Force
% glut_min2_torque(i,j) = glut_min2_f(i,j)*glut_min2_ma(i,j); %Torque
%Gluteus Minimus, 3
glut_min3_length(i,j) = norm(VecTrans(T_h,glut_min3_i) - glut_min3_o); %Total length
glut_min3_ma(i,j) = norm(cross(glut_min3_o - p_h, (VecTrans(T_h,glut_min3_i)-glut_min3_o)/norm((VecTrans(T_h,glut_min3_i)-glut_min3_o)))); %Moment arm
% glut_min3_f(i,j) = forz(glut_min3_length(i,j), glut_min3_mif, glut_min3_ofl, glut_min3_tsl, glut_min3_pa); %Force
% glut_min3_torque(i,j) = glut_min3_f(i,j)*glut_min3_ma(i,j); %Torque
%Gluteus Maximus, 1
glut_max1_length(i,j) = norm(glut_max1_wr1-glut_max1_o)+norm(VecTrans(T_h,glut_max1_wr2) - glut_max1_wr1)+norm(glut_max1_i-glut_max1_wr2); %Total length
glut_max1_ma(i,j) = norm(cross(glut_max1_o - p_h, (VecTrans(T_h,glut_max1_i)-glut_max1_o)/norm((VecTrans(T_h,glut_max1_i)-glut_max1_o)))); %Moment arm
glut_max1_f(i,j) = forz(glut_max1_length(i,j), glut_max1_mif, glut_max1_ofl, glut_max1_tsl, glut_max1_pa); %Force
glut_max1_torque(i,j) = glut_max1_f(i,j)*glut_max1_ma(i,j); %Torque
%Gluteus Maximus, 2
glut_max2_length(i,j) = norm(glut_max2_wr1-glut_max2_o)+norm(VecTrans(T_h,glut_max2_wr2) - glut_max2_wr1)+norm(glut_max2_i-glut_max2_wr2); %Total length
glut_max2_ma(i,j) = norm(cross(glut_max2_o - p_h, (VecTrans(T_h,glut_max2_i)-glut_max2_o)/norm((VecTrans(T_h,glut_max2_i)-glut_max2_o)))); %Moment arm
% glut_max2_f(i,j) = forz(glut_max2_length(i,j), glut_max2_mif, glut_max2_ofl, glut_max2_tsl, glut_max2_pa); %Force
% glut_max2_torque(i,j) = glut_max2_f(i,j)*glut_max2_ma(i,j); %Torque
%Gluteus Maximus, 3
glut_max3_length(i,j) = norm(glut_max3_wr1-glut_max3_o)+norm(VecTrans(T_h,glut_max3_wr2) - glut_max3_wr1)+norm(glut_max3_i-glut_max3_wr2); %Total length
glut_max3_ma(i,j) = norm(cross(glut_max3_o - p_h, (VecTrans(T_h,glut_max3_i)-glut_max3_o)/norm((VecTrans(T_h,glut_max3_i)-glut_max3_o)))); %Moment arm
% glut_max3_f(i,j) = forz(glut_max3_length(i,j), glut_max3_mif, glut_max3_ofl, glut_max3_tsl, glut_max3_pa); %Force
% glut_max3_torque(i,j) = glut_max3_f(i,j)*glut_max3_ma(i,j); %Torque
%Semimembranosus
semimem_length(i,j) = norm(VecTrans(T_k,semimem_wr1)-semimem_o)+norm(semimem_i-semimem_wr1); %Total length
semimem_ma(i,j) = norm(cross(semimem_o - p_k, (VecTrans(T_k,semimem_i)-semimem_o)/norm((VecTrans(T_k,semimem_i)-semimem_o)))); %Moment arm
% semimem_f(i,j) = forz(semimem_length(i,j), semimem_mif, semimem_ofl, semimem_tsl, semimem_pa); %Force
% semimem_torque(i,j) = semimem_f(i,j)*semimem_ma(i,j); %Torque
%Semitendinosus
semiten_length(i,j) = norm(VecTrans(T_k,semiten_wr1)-semiten_o)+norm(semiten_wr2-semiten_wr1)+norm(semiten_wr3-semiten_wr2)+norm(semiten_wr3-semiten_i); %Total length
semiten_ma(i,j) = norm(cross(semiten_o - p_k, (VecTrans(T_k,semiten_i)-semiten_o)/norm((VecTrans(T_k,semiten_i)-semiten_o)))); %Moment arm
% semiten_f(i,j) = forz(semiten_length(i,j), semiten_mif, semiten_ofl, semiten_tsl, semiten_pa); %Force
% semiten_torque(i,j) = semiten_f(i,j)*semiten_ma(i,j); %Torque
%Bicep Femoris, Long Head
bifemlh_length(i,j) = norm(VecTrans(T_k,bifemlh_wr1)-bifemlh_o)+norm(bifemlh_i-bifemlh_wr1); %Total length
bifemlh_ma(i,j) = norm(cross(bifemlh_o - p_k, (VecTrans(T_k,bifemlh_i)-bifemlh_o)/norm(VecTrans(T_k,bifemlh_i)-bifemlh_o))); %Moment arm
% bifemlh_f(i,j) = forz(bifemlh_length(i,j), bifemlh_mif, bifemlh_ofl, bifemlh_tsl, bifemlh_pa); %Force
% bifemlh_torque(i,j) = bifemlh_f(i,j)*bifemlh_ma(i,j); %Torque
%Bicep Femoris, Short Head
bifemsh_length(i,j) = norm(VecTrans(T_k,bifemsh_wr1)-bifemsh_o)+norm(bifemsh_i-bifemsh_wr1); %Total length, right
bifemsh_ma(i,j) = norm(cross(bifemsh_o - p_k, (VecTrans(T_k,bifemsh_i)-bifemsh_o)/norm(VecTrans(T_k,bifemsh_i)-bifemsh_o))); %Moment arm
% bifemsh_f(i,j) = forz(bifemsh_length(i,j), bifemsh_mif, bifemsh_ofl, bifemsh_tsl, bifemsh_pa); %Force
% bifemsh_torque(i,j) = bifemsh_f(i,j)*bifemsh_ma(i,j); %Torque
%Sartorius
sar_length(i,j) = norm(VecTrans(T_h,sar_wr1)-sar_o)+norm(VecTrans(T_k,sar_wr2)-sar_wr1)+norm(sar_wr3-sar_wr2)+norm(sar_wr3-sar_i); %Total length, right
sar_ma(i,j) = norm(cross(sar_o - p_h, (VecTrans(T_h*T_k,sar_i)-sar_o)/norm(VecTrans(T_h*T_k,sar_i)-sar_o))); %Moment arm
% sar_f(i,j) = forz(sar_length(i,j), sar_mif, sar_ofl, sar_tsl, sar_pa); %Force
% sar_torque(i,j) = sar_f(i,j)*sar_ma(i,j); %Torque
%Adductor Longus
add_long_length(i,j) = norm(VecTrans(T_h,add_long_i) - add_long_o); %Total length, right
add_long_ma(i,j) = norm(cross(add_long_o - p_h, (VecTrans(T_h,add_long_i)-add_long_o)/norm(VecTrans(T_h,add_long_i)-add_long_o))); %Moment arm
% add_long_f(i,j) = forz(add_long_length(i,j), add_long_mif, add_long_ofl, add_long_tsl, add_long_pa); %Force
% add_long_torque(i,j) = add_long_f(i,j)*add_long_ma(i,j); %Torque
%Adductor Brevis
add_brev_length(i,j) = norm(VecTrans(T_h,add_brev_i) - add_brev_o); %Total length
add_brev_ma(i,j) = norm(cross(add_brev_o - p_h, (VecTrans(T_h,add_brev_i)-add_brev_o)/norm(VecTrans(T_h,add_brev_i)-add_brev_o))); %Moment arm
% add_brev_f(i,j) = forz(add_brev_length(i,j), add_brev_mif, add_brev_ofl, add_brev_tsl, add_brev_pa); %Force
% add_brev_torque(i,j) = add_brev_f(i,j)*add_brev_ma(i,j); %Torque
%Adductor Magnus, 1
add_mag1_length(i,j) = norm(VecTrans(T_h,add_mag1_i) - add_mag1_o); %Total length
add_mag1_ma(i,j) = norm(cross(add_mag1_o - p_h, (VecTrans(T_h,add_mag1_i)-add_mag1_o)/norm(VecTrans(T_h,add_mag1_i)-add_mag1_o))); %Moment arm
% add_mag1_f(i,j) = forz(add_mag1_length(i,j), add_mag1_mif, add_mag1_ofl, add_mag1_tsl, add_mag1_pa); %Force
% add_mag1_torque(i,j) = add_mag1_f(i,j)*add_mag1_ma(i,j); %Torque
%Adductor Magnus, 2
add_mag2_length(i,j) = norm(VecTrans(T_h,add_mag2_i) - add_mag2_o); %Total length
add_mag2_ma(i,j) = norm(cross(add_mag2_o - p_h, (VecTrans(T_h,add_mag2_i)-add_mag2_o)/norm(VecTrans(T_h,add_mag2_i)-add_mag2_o))); %Moment arm
% add_mag2_f(i,j) = forz(add_mag2_length(i,j), add_mag2_mif, add_mag2_ofl, add_mag2_tsl, add_mag2_pa); %Force
% add_mag2_torque(i,j) = add_mag2_f(i,j)*add_mag2_ma(i,j); %Torque
%Adductor Magnus, 3
add_mag3_length(i,j) = norm(VecTrans(T_h,add_mag3_i) - add_mag3_o); %Total length
add_mag3_ma(i,j) = norm(cross(add_mag3_o - p_h, (VecTrans(T_h,add_mag3_i)-add_mag3_o)/norm(VecTrans(T_h,add_mag3_i)-add_mag3_o))); %Moment arm
% add_mag3_f(i,j) = forz(add_mag3_length(i,j), add_mag3_mif, add_mag3_ofl, add_mag3_tsl, add_mag3_pa); %Force
% add_mag3_torque(i,j) = add_mag3_f(i,j)*add_mag3_ma(i,j); %Torque
%Tensor Fasciae Latae
tfl_length(i,j) = norm(VecTrans(T_h,tfl_wr1)-tfl_o)+norm(tfl_wr2 - tfl_wr1)+norm(tfl_wr2-VecTrans(T_k,tfl_i)); %Total length, right
tfl_ma(i,j) = norm(cross(tfl_o - p_h, (VecTrans(T_h*T_k,tfl_i)-tfl_o)/norm(VecTrans(T_h*T_k,tfl_i)-tfl_o))); %Moment arm
% tfl_f(i,j) = forz(tfl_length(i,j), tfl_mif, tfl_ofl, tfl_tsl, tfl_pa); %Force
% tfl_torque(i,j) = tfl_f(i,j)*tfl_ma(i,j); %Torque
%Pectineus
pect_length(i,j) = norm(VecTrans(T_h,pect_i) - pect_o); %Total length
pect_ma(i,j) = norm(cross(pect_o - p_h, (VecTrans(T_h,pect_i)-pect_o)/norm(VecTrans(T_h,pect_i)-pect_o))); %Moment arm
% pect_f(i,j) = forz(pect_length(i,j), pect_mif, pect_ofl, pect_tsl, pect_pa); %Force
% pect_torque(i,j) = pect_f(i,j)*pect_ma(i,j); %Torque
%Gracilis
grac_length(i,j) = norm(VecTrans(T_h*T_k,grac_wr1)-grac_o)+norm(grac_wr2 - grac_wr1)+norm(grac_i-grac_wr2); %Total length, right
grac_ma(i,j) = norm(cross(grac_o - p_h, (VecTrans(T_h*T_k,grac_i)-grac_o)/norm(VecTrans(T_h*T_k,grac_i)-grac_o))); %Moment arm
% grac_f(i,j) = forz(grac_length(i,j), grac_mif, grac_ofl, grac_tsl, grac_pa); %Force
% grac_torque(i,j) = grac_f(i,j)*grac_ma(i,j); %Torque
%Iliacus
iliacus_length(i,j) = norm(iliacus_wr1-iliacus_o)+norm(iliacus_wr2-iliacus_wr1)+norm(VecTrans(T_h,iliacus_wr3)-iliacus_wr2)+norm(iliacus_wr3-iliacus_i); %Total length
iliacus_ma(i,j) = norm(cross(iliacus_o - p_h, (VecTrans(T_h,iliacus_i)-iliacus_o)/norm(VecTrans(T_h,iliacus_i)-iliacus_o))); %Moment arm
% iliacus_f(i,j) = forz(iliacus_length(i,j), iliacus_mif, iliacus_ofl, iliacus_tsl, iliacus_pa); %Force
% iliacus_torque(i,j) = iliacus_f(i,j)*iliacus_ma(i,j); %Torque
%Psoas
psoas_length(i,j) = norm(psoas_wr1-psoas_o)+norm(psoas_wr2-psoas_wr1)+norm(VecTrans(T_h,psoas_wr3)-psoas_wr2)+norm(psoas_wr3-psoas_i); %Total length, right
psoas_ma(i,j) = norm(cross(psoas_o - p_h, (VecTrans(T_h,psoas_i)-psoas_o)/norm(VecTrans(T_h,psoas_i)-psoas_o))); %Moment arm
% psoas_f(i,j) = forz(psoas_length(i,j), psoas_mif, psoas_ofl, psoas_tsl, psoas_pa); %Force
% psoas_torque(i,j) = psoas_f(i,j)*psoas_ma(i,j); %Torque
%Quadriceps Femoris
quad_fem_length(i,j) = norm(VecTrans(T_h,quad_fem_i) - quad_fem_o); %Total length
quad_fem_ma(i,j) = norm(cross(quad_fem_o - p_h, (VecTrans(T_h,quad_fem_i)-quad_fem_o)/norm(VecTrans(T_h,quad_fem_i)-quad_fem_o))); %Moment arm
% quad_fem_f(i,j) = forz(quad_fem_length(i,j), quad_fem_mif, quad_fem_ofl, quad_fem_tsl, quad_fem_pa); %Force
% quad_fem_torque(i,j) = quad_fem_f(i,j)*quad_fem_ma(i,j); %Torque
%Gemelli
gem_length(i,j) = norm(VecTrans(T_h,gem_i) - gem_o); %Total length
gem_ma(i,j) = norm(cross(gem_o - p_h, (VecTrans(T_h,gem_i)-gem_o)/norm(VecTrans(T_h,gem_i)-gem_o))); %Moment arm
% gem_f(i,j) = forz(gem_length(i,j), gem_mif, gem_ofl, gem_tsl, gem_pa); %Force
% gem_torque(i,j) = gem_f(i,j)*gem_ma(i,j); %Torque
%Piriformis
peri_length(i,j) = norm(VecTrans(T_h,peri_wr1)-peri_o)+norm(peri_i-peri_wr1); %Total length
peri_ma(i,j) = norm(cross(peri_o - p_h, (VecTrans(T_h,peri_i)-peri_o)/norm(VecTrans(T_h,peri_i)-peri_o))); %Moment arm
% peri_f(i,j) = forz(peri_length(i,j), peri_mif, peri_ofl, peri_tsl, peri_pa); %Force
% peri_torque(i,j) = peri_f(i,j)*peri_ma(i,j); %Torque
%Rectus Femoris
rect_fem_o = [-0.03; -0.031; 0.097]; %Pelvis, right, origin
rect_fem_wr1 = [0.033; -0.403; 0.002]; %Femur, right, wrapping point 1
rect_fem_i = [0.062; 0.021; 0.001]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
rect_fem_length(i,j) = norm(VecTrans(T_h,rect_fem_wr1)-rect_fem_o)+norm(VecTrans(T_k,rect_fem_i)-rect_fem_wr1); %Total length
rect_fem_ma(i,j) = norm(cross(rect_fem_o - p_h, (VecTrans(T_h*T_k,rect_fem_i)-rect_fem_o)/norm(VecTrans(T_h*T_k,rect_fem_i)-rect_fem_o))); %Moment arm
% rect_fem_f(i,j) = forz(rect_fem_length(i,j), rect_fem_mif, rect_fem_ofl, rect_fem_tsl, rect_fem_pa); %Force
% rect_fem_torque(i,j) = rect_fem_f(i,j)*rect_fem_ma(i,j); %Torque
%Vastus Medialis
vas_med_o = [0.014; -0.21; 0.019]; %Femur, right, origin
vas_med_wr1 = [0.036; -0.277; 0.001]; %Femur, right, wrapping point 1
vas_med_wr2 = [0.037; -0.405; -0.012]; %Femur, right, wrapping point 2
vas_med_wr3 = [0.027; -0.426; -0.013]; %Femur, right, wrapping point 3
vas_med_i = [0.056; 0.022; -0.015]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
vas_med_length(i,j) = norm(vas_med_wr1-vas_med_o)+norm(vas_med_wr2-vas_med_wr1)+norm(vas_med_wr3-vas_med_wr2)+norm(vas_med_wr3-VecTrans(T_k,vas_med_i)); %Total length
vas_med_ma(i,j) = norm(cross(vas_med_o - p_k, (VecTrans(T_k,vas_med_i)-vas_med_o)/norm(VecTrans(T_k,vas_med_i)-vas_med_o))); %Moment arm
% vas_med_f(i,j) = forz(vas_med_length(i,j), vas_med_mif, vas_med_ofl, vas_med_tsl, vas_med_pa); %Force
% vas_med_torque(i,j) = vas_med_f(i,j)*vas_med_ma(i,j); %Torque
%Vastus Intermedius
vas_int_o = [0.029; -0.192; 0.031]; %Femur, right, origin
vas_int_wr1 = [0.034; -0.208; 0.028]; %Femur, right, wrapping point 1
vas_int_wr2 = [0.034; -0.403; 0.006]; %Femur, right, wrapping point 2
vas_int_i = [0.055; 0.025; 0.002]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
vas_int_length(i,j) = norm(vas_int_wr1-vas_int_o)+norm(vas_int_wr2-vas_int_wr1)+norm(vas_int_wr2-VecTrans(T_k,vas_int_i)); %Total length, right
vas_int_ma(i,j) = norm(cross(vas_int_o - p_k, (VecTrans(T_k,vas_int_i)-vas_int_o)/norm(VecTrans(T_k,vas_int_i)-vas_int_o))); %Moment arm
% vas_int_f(i,j) = forz(vas_int_length(i,j), vas_int_mif, vas_int_ofl, vas_int_tsl, vas_int_pa); %Force
% vas_int_torque(i,j) = vas_int_f(i,j)*vas_int_ma(i,j); %Torque
%Vastus Lateralis
vas_lat_o = [0.005; -0.185; 0.035]; %Femur, right, origin
vas_lat_wr1 = [0.027; -0.259; 0.041]; %Femur, right, wrapping point 1
vas_lat_wr2 = [0.036; -0.403; 0.002]; %Femur, right, wrapping point 2
vas_lat_wr3 = [0.025; -0.424; 0.018]; %Femur, right, wrapping point 3
vas_lat_i = [0.06; 0.02; 0.016]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
vas_lat_length(i,j) = norm(vas_lat_wr1-vas_lat_o)+norm(vas_lat_wr2-vas_lat_wr1)+norm(vas_lat_wr3-vas_lat_wr2)+norm(vas_lat_wr3-(VecTrans(T_k,vas_lat_i))); %Total length
vas_lat_ma(i,j) = norm(cross(vas_lat_o - p_k, (VecTrans(T_k,vas_lat_i)-vas_lat_o)/norm(VecTrans(T_k,vas_lat_i)-vas_lat_o))); %Moment arm
% vas_lat_f(i,j) = forz(vas_lat_length(i,j), vas_lat_mif, vas_lat_ofl, vas_lat_tsl, vas_lat_pa); %Force
% vas_lat_torque(i,j) = vas_lat_f(i,j)*vas_lat_ma(i,j); %Torque
%Medial Gastroocnemius
%Lateral Gastrocenemius
%Soleus
%Tibialis Posterior
%Flexor Digitorum Longus
%Flexor Hallucis Longus
%Peroneus Brevis
%Peroneus Longus
%Peroneus Tertius
%Extensor Digitorum Longus
%Extensor Hallucis Longus
%Erector Spinae
ercspn_length(i,j) = norm(VecTrans(T_b,ercspn_i) - ercspn_o); %Total length, right
ercspn_ma(i,j) = norm(cross(ercspn_o - p_b, (VecTrans(T_b,ercspn_i)-ercspn_o)/norm(VecTrans(T_b,ercspn_i)-ercspn_o))); %Moment arm
% ercspn_f(i,j) = forz(ercspn_length(i,j), ercspn_mif, ercspn_ofl, ercspn_tsl, ercspn_pa); %Force
% ercspn_torque(i,j) = ercspn_f(i,j)*ercspn_ma(i,j); %Torque
%Internal Oblique
intobl_length(i,j) = norm(VecTrans(T_b,intobl_i) - intobl_o); %Total length, right
intobl_ma(i,j) = norm(cross(intobl_o - p_b, (VecTrans(T_b,intobl_i)-intobl_o)/norm(VecTrans(T_b,intobl_i)-intobl_o))); %Moment arm
% intobl_f(i,j) = forz(intobl_length(i,j), intobl_mif, intobl_ofl, intobl_tsl, intobl_pa); %Force
% intobl_torque(i,j) = intobl_f(i,j)*intobl_ma(i,j); %Torque
%External Oblique extobl
extobl_length(i,j) = norm(VecTrans(T_b,extobl_i) - extobl_o); %Total length, right
extobl_ma(i,j) = norm(cross(extobl_o - p_b, (VecTrans(T_b,extobl_i)-extobl_o)/norm(VecTrans(T_b,extobl_i)-extobl_o))); %Moment arm
% extobl_f(i,j) = forz(extobl_length(i,j), extobl_mif, extobl_ofl, extobl_tsl, extobl_pa); %Force
% extobl_torque(i,j) = extobl_f(i,j)*extobl_ma(i,j); %Torque
end
end

figure
surf(theta_hZ, theta_hX, glut_max1_length)
title('Length vs Hip flex and adduction')
ylabel('Hip Adduction/Abduction, Degrees')
xlabel('Hip Extension/Flexion, Degrees')
zlabel('Length, meters')
colorbar

figure
surf(theta_hZ, theta_hX, glut_max1_ma) 
title('Moment vs Hip flex and adduction')
ylabel('Hip Adduction/Abduction, Degrees')
xlabel('Hip Extension/Flexion, Degrees')
zlabel('Moment arm, meters')
colorbar

figure
surf(theta_hZ1, theta_hX, glut_max1_f)
title('Force vs Hip flex and adduction')
ylabel('Hip Adduction/Abduction, Degrees')
xlabel('Hip Extension/Flexion, Degrees')
zlabel('Force, Newtons')
colorbar

figure
surf(theta_hZ, theta_hX, glut_max1_torque)
title('Torque vs Hip flex and adduction')
ylabel('Hip Adduction/Abduction, Degrees')
xlabel('Hip Extension/Flexion, Degrees')
zlabel('Torque, N-m')
colorbar

[glut_med1_dia, glut_med1_size]= Size(glut_med1_length, glut_med1_mif);
[glut_med2_dia, glut_med2_size]= Size(glut_med2_length, glut_med2_mif);
[glut_med3_dia, glut_med3_size]= Size(glut_med3_length, glut_med3_mif);
[glut_min1_dia, glut_min1_size]= Size(glut_min1_length, glut_min1_mif);
[glut_min2_dia, glut_min2_size]= Size(glut_min2_length, glut_min2_mif);
[glut_min3_dia, glut_min3_size]= Size(glut_min3_length, glut_min3_mif);
[glut_max1_dia, glut_max1_size]= Size(glut_max1_length, glut_max1_mif);
[glut_max2_dia, glut_max2_size]= Size(glut_max2_length, glut_max2_mif);
[glut_max3_dia, glut_max3_size]= Size(glut_max3_length, glut_max3_mif);
[semimem_dia, semimem_size]= Size(semimem_length, semimem_mif);
[semiten_dia, semiten_size]= Size(semiten_length, semiten_mif);
[bifemlh_dia, bifemlh_size]= Size(bifemlh_length, bifemlh_mif);
[bifemsh_dia, bifemsh_size]= Size(bifemsh_length, bifemsh_mif);
[sar_dia, sar_size]= Size(sar_length, sar_mif);
[add_long_dia, add_long_size]= Size(add_long_length, add_long_mif);
[add_brev_dia, add_brev_size]= Size(add_brev_length, add_brev_mif);
[add_mag1_dia, add_mag1_size]= Size(add_mag1_length, add_mag1_mif);
[add_mag2_dia, add_mag2_size]= Size(add_mag2_length, add_mag2_mif);
[add_mag3_dia, add_mag3_size]= Size(add_mag3_length, add_mag3_mif);
[tfl_dia, tfl_size]= Size(tfl_length, tfl_mif);
[pect_dia, pect_size]= Size(pect_length, pect_mif);
[grac_dia, grac_size]= Size(grac_length, grac_mif);
[iliacus_dia, iliacus_size]= Size(iliacus_length, iliacus_mif);
[psoas_dia, psoas_size]= Size(psoas_length, psoas_mif);
[quad_fem_dia, quad_fem_size]= Size(quad_fem_length, quad_fem_mif);
[gem_dia, gem_size]= Size(gem_length, gem_mif);
[peri_dia, peri_size]= Size(peri_length, peri_mif);
[rect_fem_dia, rect_fem_size]= Size(rect_fem_length, rect_fem_mif);
[vas_med_dia, vas_med_size]= Size(vas_med_length, vas_med_mif);
[vas_int_dia, vas_int_size]= Size(vas_int_length, vas_int_mif);
[vas_lat_dia, vas_lat_size]= Size(vas_lat_length, vas_lat_mif);
[ercspn_dia, ercspn_size]= Size(ercspn_length, ercspn_mif);
[intobl_dia, intobl_size]= Size(intobl_length, intobl_mif);
[extobl_dia, extobl_size]= Size(extobl_length, extobl_mif);

