clear;clc;

theta_h_r_x = 0; %Hip, right, abduction/adduction
theta_h_r_y = 0; %Hip, right, external/internal rotation
theta_h_r_z = pi/2; %Hip, right, extension/flexion
theta_h_l_x = 0; %Hip, left, adduction/abduction
theta_h_l_y = 0; %Hip, left, internal/external rotation
theta_h_l_z = pi/2; %Hip, left, extension/flexion
theta_k_r = 0; %Knee, right, flexion/extension
theta_k_l = 0; %Knee, left, flexion/extension
theta_b_x = 0; %lumbar bending left/right
theta_b_y = 0; %lumbar rotation cw/ccw
theta_b_z = 0; %back flexion/extension

%Rotation matrices
%Hip, right
R_h_r_x = [1, 0, 0; 0, cos(theta_h_r_x), -sin(theta_h_r_x); 0, sin(theta_h_r_x), cos(theta_h_r_x)];
R_h_r_y = [cos(theta_h_r_y), 0, sin(theta_h_r_y); 0, 1, 0; -sin(theta_h_r_y), 0, cos(theta_h_r_y)];
R_h_r_z = [cos(theta_h_r_z), -sin(theta_h_r_z), 0; sin(theta_h_r_z), cos(theta_h_r_z), 0; 0, 0, 1];
R_hR = R_h_r_x*R_h_r_y*R_h_r_z;
%Hip, left
R_h_l_x = [1, 0, 0; 0, cos(theta_h_l_x), -sin(theta_h_l_x); 0, sin(theta_h_l_x), cos(theta_h_l_x)];
R_h_l_y = [cos(theta_h_l_y), 0, sin(theta_h_l_y); 0, 1, 0; -sin(theta_h_l_y), 0, cos(theta_h_l_y)];
R_h_l_z = [cos(theta_h_l_z), -sin(theta_h_l_z), 0; sin(theta_h_l_z), cos(theta_h_l_z), 0; 0, 0, 1];
R_hL = R_h_l_x*R_h_l_y*R_h_l_z;
%Knee, right
R_kR = [cos(theta_k_r), -sin(theta_k_r), 0; sin(theta_k_r), cos(theta_k_r), 0; 0, 0, 1];
%Knee, left
R_kL = [cos(theta_k_l), -sin(theta_k_l), 0; sin(theta_k_l), cos(theta_k_l), 0; 0, 0, 1];
%Back
R_b_x = [1, 0, 0; 0, cos(theta_b_x), -sin(theta_b_x); 0, sin(theta_b_x), cos(theta_b_x)];
R_b_y = [cos(theta_b_y), 0, sin(theta_b_y); 0, 1, 0; -sin(theta_b_y), 0, cos(theta_b_y)];
R_b_z = [cos(theta_b_z), -sin(theta_b_z), 0; sin(theta_b_z), cos(theta_b_z), 0; 0, 0, 1];
R_b = R_b_x*R_b_y*R_b_z;

%Joint location in home position
p_hR = [-0.0707; -0.0661; 0.0835]; %Hip, right
p_hL = [-0.0707; -0.0661; -0.0835]; %Hip, left
p_kR = [-0.00449986; -0.395821; 0]; %Knee, function of knee angle
p_kL = [-0.00449986; -0.395821; 0]; %Knee, function of knee angle
p_b = [-0.1007; 0.0815; 0]; %Back

%Transformation matrices
T_p = [1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1]; %Pelvis
T_hR = RpToTrans(R_hR,p_hR); %Hip to Pelvis, right
T_hL = RpToTrans(R_hL,p_hL); %Hip to Pelvis, left
T_kR = RpToTrans(R_kR,p_kR); %Knee to Femur, right
T_kL = RpToTrans(R_kL,p_kL); %Knee to Femur, left
T_b = RpToTrans(R_b,p_b); %Back to Pelvis

glut_max3_r_o = [-0.156; -0.031; 0.006]; %Pelvis, right, origin
glut_max3_r_wr1 = [-0.153; -0.105; 0.04]; %Pelvis, right, wrapping point 1
glut_max3_r_wr2 = [-0.03; -0.104; 0.014]; %Femur, right, wrapping point 2
glut_max3_r_i = [-0.006; -0.142; 0.041]; %Femur, right, insertion
glut_max3_r_wr2T = VecTrans(T_hR,glut_max3_r_wr2) %Transform attachment point
glut_max3_r_wr1-glut_max3_r_o
norm(glut_max3_r_wr1-glut_max3_r_o)
glut_max3_r_wr2T - glut_max3_r_wr1
norm(glut_max3_r_wr2T - glut_max3_r_wr1)
glut_max3_r_i-glut_max3_r_wr2T
norm(glut_max3_r_i-glut_max3_r_wr2T)
glut_max3_r_length = norm(-glut_max3_r_wr1+glut_max3_r_o)+norm(-glut_max3_r_wr2T + glut_max3_r_wr1)+norm(-glut_max3_r_i+glut_max3_r_wr2)