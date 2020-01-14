clear
clc
close all

%Units in m, N, rad/s unless otherwise stated.

load AttachPoints.mat  %Attachment points for muscles with stationary attachment points
load MovePoints.mat    %Moving attachment points and knee data
load MuscProp.mat      %Muscle specific properties

k = 180/pi;  %convert radians to degrees
c = pi/180;  %convert degrees to radians
    
if exist('divisions', 'var') == 0
    divisions = 100;
end

if exist('ChooseJoint', 'var') == 0
    ChooseJoint = 'Toe';
end

if isequal(ChooseJoint, 'Back')
    %Range of motion, radians
    thetamin_bx = -20*c;       %Lumbar bending  (I.A. Kapandji, Phsiology of the Joints. Churchill Livingstone, 1986)
    thetamax_bx = 20*c;
    thetamin_by = -5*c;       %Lumbar rotation  (Kapandji 1986)
    thetamax_by = 5*c;
    thetamin_bz = -60*c;       %Back flexion-/extension+  (Kapandji 1986)
    thetamax_bz = 20*c;

    %Indexing prep, see IndexPrep.m
    s = 100;
    t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)
   
    %T matrix, length, moment arm, force, and torque calculations
    for j=1:t
     for i = 1:s
    %Degree of freedom, position or range
    theta_b_x = thetamin_bx + (thetamax_bx-thetamin_bx)*(i-1)/s; %Lumbar bending left-/right+
    theta_b_y = 0; %thetamin_by + (thetamax_by-thetamin_by)*(j-1)/s; %Lumbar rotation cw-/ccw+
    theta_b_z = thetamin_bz + (thetamax_bz-thetamin_bz)*(j-1)/s; %Back flexion-/extension+

    %Index for plotting (for 2 dof, j index to x axis, i index to y axis)
    theta_bX(i,1) = theta_b_x*k;
    theta_bZ(j,1) = theta_b_z*k;

    %Rotation matrices
    %Back
    R_b_x = [1, 0, 0; 0, cos(theta_b_x), -sin(theta_b_x); 0, sin(theta_b_x), cos(theta_b_x)];
    R_b_y = [cos(theta_b_y), 0, sin(theta_b_y); 0, 1, 0; -sin(theta_b_y), 0, cos(theta_b_y)];
    R_b_z = [cos(theta_b_z), -sin(theta_b_z), 0; sin(theta_b_z), cos(theta_b_z), 0; 0, 0, 1];
    R_b = R_b_x*R_b_y*R_b_z;

    %Joint location in home position
    p_b = [-0.1007; 0.0815; 0]; %Back

    %Transformation matrices
    T_b = RpToTrans(R_b,p_b); %Back to Pelvis

    %Muscle length, moment arm, force, and torque
    %Erector Spinae
    ercspn_length(i,j) = norm(VecTrans(T_b,ercspn_i) - ercspn_o); %Total length, right
    ercspn_maX(i,j) = CrossProd(ercspn_i, (VecTrans(T_b\eye(4),ercspn_o)-ercspn_i)/norm(VecTrans(T_b\eye(4),ercspn_o)-ercspn_i),1); %Moment arm, x axis
    ercspn_maY(i,j) = CrossProd(ercspn_i, (VecTrans(T_b\eye(4),ercspn_o)-ercspn_i)/norm(VecTrans(T_b\eye(4),ercspn_o)-ercspn_i),2); %Moment arm, y axis
    ercspn_maZ(i,j) = CrossProd(ercspn_i, (VecTrans(T_b\eye(4),ercspn_o)-ercspn_i)/norm(VecTrans(T_b\eye(4),ercspn_o)-ercspn_i),3); %Moment arm, z axis
    ercspn_f(i,j) = forz(ercspn_length(i,j), ercspn_mif, ercspn_ofl, ercspn_tsl, ercspn_pa); %Force
    ercspn_tqX(i,j) = ercspn_f(i,j)*ercspn_maX(i,j); %Torque, x axis
    ercspn_tqY(i,j) = ercspn_f(i,j)*ercspn_maY(i,j); %Torque, y axis
    ercspn_tqZ(i,j) = ercspn_f(i,j)*ercspn_maZ(i,j); %Torque, z axis
    end
    end
    
    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    HumanAxis1 = theta_bZ;
    HumanAxis1Label = 'Back Flexion, Degrees';
    HumanAxis2 = theta_bX;
    HumanAxis2Label = 'Lumbar Bending, Degrees';
    HumanTorque1 = ercspn_tqZ;
    HumanTorque2 = ercspn_tqY;
    HumanTorque3 = ercspn_tqX;
    HumanTitle1 = 'Back Torque, Erector Spinae, Z Axis';
    HumanTitle2 = 'Back Torque, Erector Spinae, Y Axis';
    HumanTitle3 = 'Back Torque, Erector Spinae, X Axis';
    
    save(strcat('Human_', ChooseJoint, '_Data.mat'), 'HumanAxis1', 'HumanAxis1Label', 'HumanAxis2', 'HumanAxis2Label', 'HumanTorque1', 'HumanTorque2', 'HumanTorque3', 'HumanTitle1', 'HumanTitle2', 'HumanTitle3')

elseif isequal(ChooseJoint, 'Bi_Hip')
    %Range of motion, radians
    range_hz = c*[ -20, 140];    %Hip (Wikipedia)
    thetamin_hz = range_hz(1);      %20 deg extention
    thetamax_hz = range_hz(2);      %140 deg flexion
    range_k = [-2.0943951, 0.17453293];    %Knee
    thetamin_k = range_k(1);       %120 deg flexion
    thetamax_k = range_k(2);       %10 deg extension

    %Indexing prep, see IndexPrep.m
    s = 100;
    t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)

    for j=1:t
     for i = 1:s
    %Degree of freedom, position or range
    theta_h_z = thetamin_hz + (thetamax_hz-thetamin_hz)*(j-1)/s; %Hip, right, extension-/flexion+
    thetamin_hy = c*(-32.5 - (theta_h_z/(range_hz(2)-range_hz(1)))*20); %Max External Rotation (linear interpolation of Wikipedia RoM)
    thetamax_hy = c*(40);  %Max Internal Rotation (Wikipedia)
    theta_h_y = 0; %Hip, right, external-/internal+ rotation
    thetamin_hx = c*(-53.75 - (theta_h_z/(range_hz(2)-range_hz(1)))*30); %Max Abduction (linear interpolation of Wikipedia RoM)
    thetamax_hx = c*(28.75 - (theta_h_z/(range_hz(2)-range_hz(1)))*10); %Max Adduction (linear interpolation of Wikipedia RoM)
    theta_h_x = 0; %thetamin_hx + (thetamax_hx-thetamin_hx)*(i-1)/s; %Hip, right, adduction-/abduction+
    theta_k = thetamin_k + (thetamax_k-thetamin_k)*(i-1)/s; %Knee, right, flexion-/extension+ 

    %Index for plotting (for 2 dof, j index to x axis, i index to y axis)
    theta_hZ(j,1) = theta_h_z*k; %Convert from radians to degrees
    theta_K(i,1) = theta_k*k;

    %Rotation matrices
    %Hip
    R_h_x = [1, 0, 0; 0, cos(theta_h_x), -sin(theta_h_x); 0, sin(theta_h_x), cos(theta_h_x)];
    R_h_y = [cos(theta_h_y), 0, sin(theta_h_y); 0, 1, 0; -sin(theta_h_y), 0, cos(theta_h_y)];
    R_h_z = [cos(theta_h_z), -sin(theta_h_z), 0; sin(theta_h_z), cos(theta_h_z), 0; 0, 0, 1]; 
    R_h = R_h_x*R_h_y*R_h_z;
    %Knee
    R_k = [cos(theta_k), -sin(theta_k), 0; sin(theta_k), cos(theta_k), 0; 0, 0, 1];

    %Joint location in home position
    p_h = [-0.0707; -0.0661; 0.0835]; %Hip, right
    p_k = [fcn1(theta_k); fcn2(theta_k); 0]; %Knee, function of knee angle (see MovePoints)

    %Transformation matrices
    T_h = RpToTrans(R_h,p_h); %Hip to Pelvis, right
    T_k = RpToTrans(R_k,p_k); %Knee to Hip, right

    %Muscle length, moment arm, force, and torque
    %Bicep Femoris, Long Head
    bifemlh_length(i,j) = norm(VecTrans(T_h*T_k,bifemlh_wr1)-bifemlh_o)+norm(bifemlh_i-bifemlh_wr1); %Total length
    bifemlh_maHX(i,j) = CrossProd(VecTrans(T_k,bifemlh_wr1), (VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1))/norm(VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1)),1); %Moment arm, Hip, x axis
    bifemlh_maHY(i,j) = CrossProd(VecTrans(T_k,bifemlh_wr1), (VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1))/norm(VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1)),2); %Moment arm, Hip, y axis
    bifemlh_maHZ(i,j) = CrossProd(VecTrans(T_k,bifemlh_wr1), (VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1))/norm(VecTrans(T_h\eye(4),bifemlh_o)-VecTrans(T_k,bifemlh_wr1)),3); %Moment arm, Hip, z axis
    bifemlh_maKZ(i,j) = CrossProd(bifemlh_wr1, (VecTrans((T_h*T_k)\eye(4),bifemlh_o)-bifemlh_wr1)/norm(VecTrans((T_h*T_k)\eye(4),bifemlh_o)-bifemlh_wr1),3); %Moment arm, Knee, z axis
    bifemlh_f(i,j) = forz(bifemlh_length(i,j), bifemlh_mif, bifemlh_ofl, bifemlh_tsl, bifemlh_pa); %Force
    bifemlh_tqHX(i,j) = bifemlh_f(i,j)*bifemlh_maHX(i,j); %Torque, Hip, x axis
    bifemlh_tqHY(i,j) = bifemlh_f(i,j)*bifemlh_maHY(i,j); %Torque, Hip, y axis
    bifemlh_tqHZ(i,j) = bifemlh_f(i,j)*bifemlh_maHZ(i,j); %Torque, Hip, z axis
    bifemlh_tqKZ(i,j) = bifemlh_f(i,j)*bifemlh_maKZ(i,j); %Torque, Knee, z axis

    end
    end
    
    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    HumanAxis1 = theta_hZ;
    HumanAxis1Label = 'Hip Flexion, Degrees';
    HumanAxis2 = theta_K;
    HumanAxis2Label = 'Knee Flexion, Degrees';
    HumanTorque1 = bifemlh_tqHX;
    HumanTorque2 = bifemlh_tqHY;
    HumanTorque3 = bifemlh_tqHZ;
    HumanTorque4 = bifemlh_tqKZ;
    HumanTitle1 = 'Hip Torque, bifemlh, X axis';
    HumanTitle2 = 'Hip Torque, bifemlh, Y axis';
    HumanTitle3 = 'Hip Torque, bifemlh, Z axis';
    HumanTitle4 = 'Knee Torque, bifemlh, Z axis';
    
    save(strcat('Human_', ChooseJoint, '_Data.mat'), 'HumanAxis1', 'HumanAxis1Label', 'HumanAxis2', 'HumanAxis2Label', 'HumanTorque1', 'HumanTorque2', 'HumanTorque3', 'HumanTorque4', 'HumanTitle1', 'HumanTitle2', 'HumanTitle3', 'HumanTitle4')

elseif isequal(ChooseJoint, 'Calves')
    %Range of motion, radians
    range_k = c*[-120, 10];    %Knee
    thetamin_k = range_k(1);       %120 deg flexion
    thetamax_k = range_k(2);       %10 deg extension
    thetamin_a = -50*c;       %ankle
    thetamax_a = 20*c;
    thetamin_s = -25*c;       %subtalar
    thetamax_s = 35*c;
    theta_s = 0;

    %Indexing prep, see IndexPrep.m
    s = 100;
    t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)
    theta_h_x = 0; %thetamin_hx + (thetamax_hx-thetamin_hx)*(i-1)/s; %Hip, right, adduction-/abduction+

    for j=1:t
     for i = 1:s
    %Degree of freedom, position or range
    theta_k = thetamin_k + (thetamax_k-thetamin_k)*(i-1)/s; %Knee, right, flexion-/extension+ 
    theta_a = thetamin_a + (thetamax_a-thetamin_a)*(j-1)/s; %Ankle, right, plantar flexion-/dorsiflexion+

    %Index for plotting (for 2 dof, j index to x axis, i index to y axis)
    theta_K(i,1) = theta_k*k;
    theta_A(j,1) = theta_a*k;

    %Rotation matrices
    %Knee
    R_k = [cos(theta_k), -sin(theta_k), 0; sin(theta_k), cos(theta_k), 0; 0, 0, 1];

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


    %Joint location in home position
    p_k = [fcn1(theta_k); fcn2(theta_k); 0]; %Knee, function of knee angle (see MovePoints)
    p_a = [0; -0.43; 0]; %Ankle, right
    p_s = [-0.04877; -0.04195; 0.00792]; %Subtalar, right

    %Transformation matrices
    T_k = RpToTrans(R_k,p_k); %Knee to Hip, right
    T_a = RpToTrans(R_a,p_a); %Ankle to Knee, right
    T_s = RpToTrans(R_s,p_s); %Subtalar to Ankle, right

    %Muscle length, moment arm, force, and torque
    % %Medial Gastroocnemius
    med_gas_length(i,j) = norm(VecTrans(T_k*T_a*T_s,med_gas_i)-med_gas_o); %Total length
    med_gas_maKZ(i,j) = CrossProd(VecTrans(T_a*T_s,med_gas_i), (VecTrans(T_k\eye(4),med_gas_o)-VecTrans(T_a*T_s,med_gas_i))/norm(VecTrans(T_k\eye(4),med_gas_o)-VecTrans(T_a*T_s,med_gas_i)),3); %Moment arm, Knee, z axis
    med_gas_maAZ(i,j) = CrossProd(VecTrans(T_s,med_gas_i), (VecTrans((T_k*T_a)\eye(4),med_gas_o)-VecTrans(T_s,med_gas_i))/norm(VecTrans((T_k*T_a)\eye(4),med_gas_o)-VecTrans(T_s,med_gas_i)),3,R_a); %Moment arm, Ankle, z' axis
    med_gas_maSX(i,j) = CrossProd(med_gas_i, (VecTrans((T_k*T_a*T_s)\eye(4),med_gas_o)-med_gas_i)/norm(VecTrans((T_k*T_a*T_s)\eye(4),med_gas_o)-med_gas_i),1,R_s); %Moment arm, Subtalar, x' axis
    med_gas_f(i,j) = forz(med_gas_length(i,j), med_gas_mif, med_gas_ofl, med_gas_tsl, med_gas_pa); %Force
    med_gas_tqKZ(i,j) = med_gas_f(i,j)*med_gas_maKZ(i,j); %Torque, Knee, z axis
    med_gas_tqAZ(i,j) = med_gas_f(i,j)*med_gas_maAZ(i,j); %Torque, Ankle, z' axis
    med_gas_tqSX(i,j) = med_gas_f(i,j)*med_gas_maSX(i,j); %Torque, Subtalar, x' axis
    %Lateral Gastrocenemius
    lat_gas_length(i,j) = norm(VecTrans(T_k*T_a*T_s,lat_gas_i)-lat_gas_o); %Total length
    lat_gas_maKZ(i,j) = CrossProd(VecTrans(T_a*T_s,lat_gas_i), (VecTrans(T_k\eye(4),lat_gas_o)-VecTrans(T_a*T_s,lat_gas_i))/norm(VecTrans(T_k\eye(4),lat_gas_o)-VecTrans(T_a*T_s,lat_gas_i)),3); %Moment arm, Knee, z axis
    lat_gas_maAZ(i,j) = CrossProd(VecTrans(T_s,lat_gas_i), (VecTrans((T_k*T_a)\eye(4),lat_gas_o)-VecTrans(T_s,lat_gas_i))/norm(VecTrans((T_k*T_a)\eye(4),lat_gas_o)-VecTrans(T_s,lat_gas_i)),3,R_a); %Moment arm, Ankle, z' axis
    lat_gas_maSX(i,j) = CrossProd(lat_gas_i, (VecTrans((T_k*T_a*T_s)\eye(4),lat_gas_o)-lat_gas_i)/norm(VecTrans((T_k*T_a*T_s)\eye(4),lat_gas_o)-lat_gas_i),1,R_s); %Moment arm, Subtalar, x' axis
    lat_gas_f(i,j) = forz(lat_gas_length(i,j), lat_gas_mif, lat_gas_ofl, lat_gas_tsl, lat_gas_pa); %Force
    lat_gas_tqKZ(i,j) = lat_gas_f(i,j)*lat_gas_maKZ(i,j); %Torque, Knee, z axis
    lat_gas_tqAZ(i,j) = lat_gas_f(i,j)*lat_gas_maAZ(i,j); %Torque, Ankle, z' axis
    lat_gas_tqSX(i,j) = lat_gas_f(i,j)*lat_gas_maSX(i,j); %Torque, Subtalar, x' axis
    %Soleus
    soleus_length(i,j) = norm(VecTrans(T_a*T_s,soleus_i)-soleus_o); %Total length
    soleus_maAZ(i,j) = CrossProd(VecTrans(T_s,soleus_i), (VecTrans(T_a\eye(4),soleus_o)-VecTrans(T_s,soleus_i))/norm(VecTrans(T_a\eye(4),soleus_o)-VecTrans(T_s,soleus_i)),3,R_a); %Moment arm, Ankle, z' axis
    soleus_maSX(i,j) = CrossProd(soleus_i, (VecTrans((T_a*T_s)\eye(4),soleus_o)-soleus_i)/norm(VecTrans((T_a*T_s)\eye(4),soleus_o)-soleus_i),1,R_s); %Moment arm, Subtalar, x' axis
    soleus_f(i,j) = forz(soleus_length(i,j), soleus_mif, soleus_ofl, soleus_tsl, soleus_pa); %Force
    soleus_tqAZ(i,j) = soleus_f(i,j)*soleus_maAZ(i,j); %Torque, Ankle, z' axis
    soleus_tqSX(i,j) = soleus_f(i,j)*soleus_maSX(i,j); %Torque, Subtalar, x' axis

    end
    end

    calves_tqKZ = med_gas_tqKZ + lat_gas_tqKZ;
    calves_tqAZ = med_gas_tqAZ + lat_gas_tqAZ + soleus_tqAZ;
    calves_tqSX = med_gas_tqSX + lat_gas_tqSX + soleus_tqSX;
    
    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    HumanAxis1 = theta_A;
    HumanAxis1Label = 'Ankle Flexion, Degrees';
    HumanAxis2 = theta_K;
    HumanAxis2Label = 'Knee Flexion, Degrees';
    HumanTorque1 = real(calves_tqKZ);
    HumanTorque2 = real(med_gas_tqAZ + lat_gas_tqAZ);
    HumanTorque3 = real(soleus_tqAZ);
    HumanTitle1 = 'Knee Torque, Gastrocnemius, Z axis';
    HumanTitle2 = 'Ankle Torque, Gastrocnemius, Z` axis';
    HumanTitle3 = 'Ankle Torque, Soleus, Z` axis';
    
    save(strcat('Human_', ChooseJoint, '_Data.mat'), 'HumanAxis1', 'HumanAxis1Label', 'HumanAxis2', 'HumanAxis2Label', 'HumanTorque1', 'HumanTorque2', 'HumanTorque3', 'HumanTitle1', 'HumanTitle2', 'HumanTitle3')

elseif isequal(ChooseJoint, 'Foot')
    %Range of motion, radians
    thetamin_a = -50*c;       %ankle
    thetamax_a = 20*c;
    thetamin_s = -25*c;       %subtalar
    thetamax_s = 35*c;

    %Indexing prep, see IndexPrep.m
    s = 100;
    t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)

    for j=1:t
     for i = 1:s
    %Degree of freedom, position or range
    theta_a = thetamin_a + (thetamax_a-thetamin_a)*(j-1)/s; %Ankle, right, plantar flexion-/dorsiflexion+
    theta_s = thetamin_s + (thetamax_s-thetamin_s)*(i-1)/s; %Subtalar, right, inversion-/eversion+ rotation

    %Index for plotting (for 2 dof, j index to x axis, i index to y axis)
    theta_A(j,1) = theta_a*k;
    theta_S(i,1) = theta_s*k;

    %Rotation matrices
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

    %Joint location in home position
    p_a = [0; -0.43; 0]; %Ankle, right
    p_s = [-0.04877; -0.04195; 0.00792]; %Subtalar, right

    %Transformation matrices
    T_a = RpToTrans(R_a,p_a); %Ankle to Knee, right
    T_s = RpToTrans(R_s,p_s); %Subtalar to Ankle, right

    %Muscle length, moment arm, force, and torque
    %Tibialis Posterior
    tib_post_length(i,j) = norm(tib_post_wr1-tib_post_o)+norm(VecTrans(T_a*T_s,tib_post_wr2)-tib_post_wr1)+norm(tib_post_i-tib_post_wr2); %Total length
    tib_post_maAZ(i,j) = CrossProd(VecTrans(T_s,tib_post_wr2), (VecTrans(T_a\eye(4),tib_post_wr1)-VecTrans(T_s,tib_post_wr2))/norm(VecTrans(T_a\eye(4),tib_post_wr1)-VecTrans(T_s,tib_post_wr2)),3,R_a); %Moment arm, Ankle, z' axis
    tib_post_maSX(i,j) = CrossProd(tib_post_wr2, (VecTrans((T_a*T_s)\eye(4),tib_post_wr1)-tib_post_wr2)/norm(VecTrans((T_a*T_s)\eye(4),tib_post_wr1)-tib_post_wr2),1,R_s); %Moment arm, subtalar, x' axis
    tib_post_f(i,j) = forz(tib_post_length(i,j), tib_post_mif, tib_post_ofl, tib_post_tsl, tib_post_pa); %Force
    tib_post_tqAZ(i,j) = tib_post_f(i,j)*tib_post_maAZ(i,j); %Torque, Ankle, z' axis
    tib_post_tqSX(i,j) = tib_post_f(i,j)*tib_post_maSX(i,j); %Torque, Subtalar, x' axis
    %Tibialis Anterior
    tib_ant_length(i,j) = norm(tib_ant_wr1-tib_ant_o)+norm(VecTrans(T_a*T_s,tib_ant_i)-tib_ant_wr1); %Total length
    tib_ant_maAZ(i,j) = CrossProd(VecTrans(T_s,tib_ant_i), (VecTrans(T_a\eye(4),tib_ant_wr1)-VecTrans(T_s,tib_ant_i))/norm(VecTrans(T_a\eye(4),tib_ant_wr1)-VecTrans(T_s,tib_ant_i)),3,R_a); %Moment arm, Ankle, z' axis
    tib_ant_maSX(i,j) = CrossProd(tib_ant_i, (VecTrans((T_a*T_s)\eye(4),tib_ant_wr1)-tib_ant_i)/norm(VecTrans((T_a*T_s)\eye(4),tib_ant_wr1)-tib_ant_i),1,R_s); %Moment arm, subtalar, x' axis
    tib_ant_f(i,j) = forz(tib_ant_length(i,j), tib_ant_mif, tib_ant_ofl, tib_ant_tsl, tib_ant_pa); %Force
    tib_ant_tqAZ(i,j) = tib_ant_f(i,j)*tib_ant_maAZ(i,j); %Torque, Ankle, z' axis
    tib_ant_tqSX(i,j) = tib_ant_f(i,j)*tib_ant_maSX(i,j); %Torque, Subtalar, x' axis
    %Peroneus Brevis
    per_brev_wr2 = [-0.014; -0.429; 0.029]; %Tibia, right, wrapping point 2
    per_brev_wr3 = [0.047; 0.027; 0.023]; %Calcn, right, wrapping point 3
    per_brev_length(i,j) = norm(per_brev_o-per_brev_wr1)+norm(per_brev_wr1-per_brev_wr2)+norm(VecTrans(T_a*T_s,per_brev_wr3)-per_brev_wr2)+norm(per_brev_wr3-per_brev_i); %Total length
    per_brev_maAZ(i,j) = CrossProd(VecTrans(T_s,per_brev_wr3), (VecTrans(T_a\eye(4),per_brev_wr2)-VecTrans(T_s,per_brev_wr3))/norm(VecTrans(T_a\eye(4),per_brev_wr2)-VecTrans(T_s,per_brev_wr3)),3,R_a); %Moment arm, Ankle, z' axis
    per_brev_maSX(i,j) = CrossProd(per_brev_wr3, (VecTrans((T_a*T_s)\eye(4),per_brev_wr2)-per_brev_wr3)/norm(VecTrans((T_a*T_s)\eye(4),per_brev_wr2)-per_brev_wr3),1,R_s); %Moment arm, subtalar, x' axis
    per_brev_f(i,j) = forz(per_brev_length(i,j), per_brev_mif, per_brev_ofl, per_brev_tsl, per_brev_pa); %Force
    per_brev_tqAZ(i,j) = per_brev_f(i,j)*per_brev_maAZ(i,j); %Torque, Ankle, z' axis
    per_brev_tqSX(i,j) = per_brev_f(i,j)*per_brev_maSX(i,j); %Torque, Subtalar, x' axis
    %Peroneus Longus
    per_long_length(i,j) = norm(per_long_o-per_long_wr1)+norm(per_long_wr1-per_long_wr2)+norm(VecTrans(T_a*T_s,per_long_wr3)-per_long_wr2)+norm(per_long_wr3-per_long_wr4)+norm(per_long_wr5-per_long_wr4)+norm(per_long_wr5-per_long_i); %Total length
    per_long_maAZ(i,j) = CrossProd(VecTrans(T_s,per_long_wr3), (VecTrans(T_a\eye(4),per_long_wr2)-VecTrans(T_s,per_long_wr3))/norm(VecTrans(T_a\eye(4),per_long_wr2)-VecTrans(T_s,per_long_wr3)),3,R_a); %Moment arm, Ankle, z' axis
    per_long_maSX(i,j) = CrossProd(per_long_wr3, (VecTrans((T_a*T_s)\eye(4),per_long_wr2)-per_long_wr3)/norm(VecTrans((T_a*T_s)\eye(4),per_long_wr2)-per_long_wr3),1,R_s); %Moment arm, subtalar, x' axis
    per_long_f(i,j) = forz(per_long_length(i,j), per_long_mif, per_long_ofl, per_long_tsl, per_long_pa); %Force
    per_long_tqAZ(i,j) = per_long_f(i,j)*per_long_maAZ(i,j); %Torque, Ankle, z' axis
    per_long_tqSX(i,j) = per_long_f(i,j)*per_long_maSX(i,j); %Torque, Subtalar, x' axis
    %Peroneus Tertius
    per_tert_length(i,j) = norm(per_tert_wr1-per_tert_o)+norm(VecTrans(T_a*T_s,per_tert_i)-per_tert_wr1) ; %Total length
    per_tert_maAZ(i,j) = CrossProd(VecTrans(T_s,per_tert_i), (VecTrans(T_a\eye(4),per_tert_wr1)-VecTrans(T_s,per_tert_i))/norm(VecTrans(T_a\eye(4),per_tert_wr1)-VecTrans(T_s,per_tert_i)),3,R_a); %Moment arm, Ankle, z' axis
    per_tert_maSX(i,j) = CrossProd(per_tert_i, (VecTrans((T_a*T_s)\eye(4),per_tert_wr1)-per_tert_i)/norm(VecTrans((T_a*T_s)\eye(4),per_tert_wr1)-per_tert_i),1,R_s); %Moment arm, subtalar, x' axis
    per_tert_f(i,j) = forz(per_tert_length(i,j), per_tert_mif, per_tert_ofl, per_tert_tsl, per_tert_pa); %Force
    per_tert_tqAZ(i,j) = per_tert_f(i,j)*per_tert_maAZ(i,j); %Torque, Ankle, z' axis
    per_tert_tqSX(i,j) = per_tert_f(i,j)*per_tert_maSX(i,j); %Torque, Subtalar, x' axis

    end
    end

    foot_tqAZ = tib_post_tqAZ + tib_ant_tqAZ + per_brev_tqAZ + per_long_tqAZ + per_tert_tqAZ;
    foot_tqSX = tib_post_tqSX + tib_ant_tqSX + per_brev_tqSX + per_long_tqSX + per_tert_tqSX;
    
    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    HumanAxis1 = theta_A;
    HumanAxis1Label = 'Ankle Flexion, Degrees';
    HumanAxis2 = theta_S;
    HumanAxis2Label = 'Subtalar Eversion, Degrees';
    HumanTorque1 = real(foot_tqAZ);
    HumanTorque2 = real(foot_tqSX);
    HumanTitle1 = 'Knee Torque, Gastrocnemius, Z axis';
    HumanTitle2 = 'Ankle Torque, Gastrocnemius, Z` axis';
    
    save(strcat('Human_', ChooseJoint, '_Data.mat'), 'HumanAxis1', 'HumanAxis1Label', 'HumanAxis2', 'HumanAxis2Label', 'HumanTorque1', 'HumanTorque2', 'HumanTitle1', 'HumanTitle2')

elseif isequal(ChooseJoint, 'Toe')
    %Only part of the script that follows new code
    
    %Create a Joint object that calculates things like transformation
    %matrix
    %Note: Currently only goes for the matrices at minimum theta. Will
    %implement some iteration process to go from minimum to maximum. 
        
    Raxis = [-0.10501355; -0.17402245; 0.97912632];
    MaxTheta = 20*pi/180;
    MinTheta = -50*pi/180;
    Home = [0; -0.43; 0];                   %Home position from the Tibia
    Ankle = JointData('Ankle', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    Raxis = [0.7871796; 0.60474746; -0.12094949];
    MaxTheta = 0;
    MinTheta = 0;               %Ben's script lists the min and max, however only 0 is used. Why?
    Home = [-0.04877; -0.04195; 0.00792];   %Home position from the Talus
    Subtalar = JointData('Subtalar', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    Raxis = [-0.5809544; 0; 0.81393611];
    MaxTheta = 80*pi/180;
    MinTheta = -30*pi/180;
    Home = [0.1788; -0.002; 0.00108];       %Home Position from the Calcn
    MTP = JointData('MTP', Raxis, MaxTheta, MinTheta, Home, divisions);
    
    T_a = Ankle.TransformationMat(:, :, :);
    T_s = Subtalar.TransformationMat(:, :, :);
    T_m = MTP.TransformationMat(:, :, :);
    
    %Create a Mesh of the Ankle and MTP. First two dimensions are the
    %values of the transformation matrics. The third dimension determines
    %the joint being observed, and the fourth dimension is its location
    %in the mesh, which will need to be converted from a one dimensional
    %array to two dimensions
    T = zeros(4, 4, 3, divisions^2);
    ii = 1;
    for i = 1:divisions
        for j = 1:divisions
            T(:, :, 1, ii) = T_a(:, :, i);            %Keep the ankle at one point as the MTP changes
            T(:, :, 2, ii) = T_s(:, :, i);            %Subtalar doesn't change
            T(:, :, 3, ii) = T_m(:, :, j);            %Change the MTP while the subtalar stays the same
            ii = ii + 1;
        end
    end
    
    %The Location matrix is created from AttachPoints.m. in this case, the
    %column vectors are concatenated into a 3xn matrix, where n is the
    %number of wrapping points + origin point + insertion point. With this
    %method, we do lose some data as to the general location of where the
    %wrapping point is located, i.e. 'Calcn, right'
    
    %Flexus Digitorum Longus
    Location1 = [-0.008, -0.015, 0.044, 0.071, 0.166, -0.002, 0.028, 0.044;
                 -0.205, -0.405, 0.032, 0.018, -0.008, -0.008, -0.007, -0.006;
                 0.002, -0.02, -0.028, -0.026, 0.012, 0.015, 0.022, 0.024];
    ViaPoints1 = [3, 3, 6];                    %Via points are the points where a transformation matrix is needed. Typically wrap point + 1
    Axis1 = [3 1 3];                           %The axis of interest when calculating the moment arm about each joint
    Muscle1 = MuscleData('Flexor Digitorum Longus', Location1, ViaPoints1, 310, 0.034, 0.4, 0.12217305, T, Axis1);

    
    %Flexor Hallucis Longus
    Location2 = [-0.008, -0.019, 0.037, 0.104, 0.173, 0.016, 0.056;
        -0.233, -0.408, 0.028, 0.007, -0.005, -0.006, -0.01;
        0.024, -0.017, -0.024, -0.026, -0.027, -0.026, -0.018];
    ViaPoints2 = [3, 3, 6];                    %Via points are the points where a transformation matrix is needed. Typically wrap point + 1
    Axis2 = [3 1 3];
    Muscle2 = MuscleData('Flexor Hallucis Longus', Location2, ViaPoints2, 322, 0.043, 0.38, 0.17453293, T, Axis2);
    
    %Extensor Digitorum Longus
    Location3 = [0.003, 0.029, 0.092, 0.162, 0, 0.044;
                -0.138, -0.401, 0.039, 0.006, 0.005, 0.0;
                0.028, 0.007, 0.0, 0.013, 0.015, 0.025];
    ViaPoints3 = [3, 3, 5];
    Axis3 = [3 1 3];
    Muscle3 = MuscleData('Extensor Digitorum Longus', Location3, ViaPoints3, 512, 0.102, 0.345, 0.13962634, T, Axis3);
    
    %Extensor Hallucis Longus
    Location4 = [0.001, 0.033, 0.097, 0.129, 0.173, 0.03, 0.056;
                -0.177, -0.398, 0.039, 0.031, 0.014, 0.004, 0.003;
                0.023, -0.008, -0.021, -0.026, -0.028, -0.024, -0.019];
    ViaPoints4 = [3, 3, 6];
    Axis4 = [3 1 3];
    Muscle4 = MuscleData('Extensor Hallucis Longus', Location4, ViaPoints4, 162, 0.111, 0.305, 0.10471976, T, Axis4);


    AnkleTorque = Muscle1.Torque(1, :, :) + Muscle2.Torque(1, :, :) + Muscle3.Torque(1, :, :) + Muscle4.Torque(1, :, :);
    SubtalarTorque = Muscle1.Torque(2, :, :) + Muscle2.Torque(2, :, :) + Muscle3.Torque(2, :, :) + Muscle4.Torque(2, :, :);
    MTPTorque = Muscle1.Torque(3, :, :) + Muscle2.Torque(3, :, :) + Muscle3.Torque(3, :, :) + Muscle4.Torque(3, :, :);
    
    %Create the Mesh of Torques to corespond with the joint angles
    for i = 1:divisions
        AnkleTorqueM(:, i) = AnkleTorque(((i-1)*divisions)+1:i*divisions);
        SubtalarTorqueM(:, i) = SubtalarTorque(((i-1)*divisions)+1:i*divisions);
        MTPTorqueM(:, i) = MTPTorque(((i-1)*divisions)+1:i*divisions);
    end
    
    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    HumanAxis1 = Ankle.Theta;
    HumanAxis1Label = 'Ankle Flexion, Degrees';
    HumanAxis2 = MTP.Theta;
    HumanAxis2Label = 'MTP Flexion, Degrees';
    HumanTorque1 = AnkleTorqueM;
    HumanTorque2 = SubtalarTorqueM;
    HumanTorque3 = MTPTorqueM;
    HumanTitle1 = 'Ankle Torque, Z Axis';
    HumanTitle2 = 'Subtalar Torque, X Axis';
    HumanTitle3 = 'MTP Torque, Z Axis';
    
    save(strcat('Human_', ChooseJoint, '_Data.mat'), 'HumanAxis1', 'HumanAxis1Label', 'HumanAxis2', 'HumanAxis2Label', 'HumanTorque1', 'HumanTorque2', 'HumanTorque3', 'HumanTitle1', 'HumanTitle2', 'HumanTitle3')

    
elseif isequal(ChooseJoint, 'Uni_Hip')
    %Range of motion, radians
    range_hz = c*[ -20, 140];    %Hip (Wikipedia)
    thetamin_hz = range_hz(1);      %20 deg extention
    thetamax_hz = range_hz(2);      %140 deg flexion
    thetamin_bx = -20*c;       %Lumbar bending  (I.A. Kapandji, Phsiology of the Joints. Churchill Livingstone, 1986)
    thetamax_bx = 20*c;
    thetamin_by = -5*c;       %Lumbar rotation  (Kapandji 1986)
    thetamax_by = 5*c;
    thetamin_bz = -60*c;       %Back flexion-/extension+  (Kapandji 1986)
    thetamax_bz = 20*c;

    %Indexing prep, see IndexPrep.m
    s = 100;
    t = s; %Either the same value as s (3-D plot) or set equal to 1 (2-D plot)

    for j=1:t
     for i = 1:s
    %Degree of freedom, position or range
    theta_h_z = thetamin_hz + (thetamax_hz-thetamin_hz)*(j-1)/s; %Hip, right, extension-/flexion+
    thetamin_hy = c*(-32.5 - (theta_h_z/(range_hz(2)-range_hz(1)))*20); %Max External Rotation (linear interpolation of Wikipedia RoM)
    thetamax_hy = c*(40);  %Max Internal Rotation (Wikipedia)
    theta_h_y = 0; %Hip, right, external-/internal+ rotation
    thetamin_hx = c*(-53.75 - (theta_h_z/(range_hz(2)-range_hz(1)))*30); %Max Abduction (linear interpolation of Wikipedia RoM)
    thetamax_hx = c*(28.75 - (theta_h_z/(range_hz(2)-range_hz(1)))*10); %Max Adduction (linear interpolation of Wikipedia RoM)
    theta_h_x = thetamin_hx + (thetamax_hx-thetamin_hx)*(i-1)/s; %Hip, right, adduction-/abduction+
    theta_b_x = 0; %thetamin_bx + (thetamax_bx-thetamin_bx)*(i-1)/s; %Lumbar bending left-/right+
    theta_b_y = 0; %thetamin_by + (thetamax_by-thetamin_by)*(j-1)/s; %Lumbar rotation cw-/ccw+
    theta_b_z = 0; %thetamin_bz + (thetamax_bz-thetamin_bz)*(j-1)/s; %Back flexion-/extension+

    %Index for plotting (for 2 dof, j index to x axis, i index to y axis)
    theta_hZ(j,1) = theta_h_z*k; %Convert from radians to degrees
    theta_hX(i,1) = theta_h_x*k;

    %Rotation matrices
    %Hip
    R_h_x = [1, 0, 0; 0, cos(theta_h_x), -sin(theta_h_x); 0, sin(theta_h_x), cos(theta_h_x)];
    R_h_y = [cos(theta_h_y), 0, sin(theta_h_y); 0, 1, 0; -sin(theta_h_y), 0, cos(theta_h_y)];
    R_h_z = [cos(theta_h_z), -sin(theta_h_z), 0; sin(theta_h_z), cos(theta_h_z), 0; 0, 0, 1]; 
    R_h = R_h_x*R_h_y*R_h_z;

    %Back
    R_b_x = [1, 0, 0; 0, cos(theta_b_x), -sin(theta_b_x); 0, sin(theta_b_x), cos(theta_b_x)];
    R_b_y = [cos(theta_b_y), 0, sin(theta_b_y); 0, 1, 0; -sin(theta_b_y), 0, cos(theta_b_y)];
    R_b_z = [cos(theta_b_z), -sin(theta_b_z), 0; sin(theta_b_z), cos(theta_b_z), 0; 0, 0, 1];
    R_b = R_b_x*R_b_y*R_b_z;

    %Joint location in home position
    p_h = [-0.0707; -0.0661; 0.0835]; %Hip, right
    p_b = [-0.1007; 0.0815; 0]; %Back

    %Transformation matrices
    T_p = [1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1]; %Pelvis
    T_h = RpToTrans(R_h,p_h); %Hip to Pelvis, right

    %Muscle length, moment arm, force, and torque
    % %Gluteus Medius, 1
    glut_med1_length(i,j) = norm(VecTrans(T_h,glut_med1_i) - glut_med1_o); %Total length
    glut_med1_maX(i,j) = CrossProd(glut_med1_i, (VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i)/norm(VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i),1); %Moment arm, x axis
    glut_med1_maY(i,j) = CrossProd(glut_med1_i, (VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i)/norm(VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i),2); %Moment arm, y axis
    glut_med1_maZ(i,j) = CrossProd(glut_med1_i, (VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i)/norm(VecTrans(T_h\eye(4),glut_med1_o)-glut_med1_i),3); %Moment arm, z axis
    glut_med1_f(i,j) = forz(glut_med1_length(i,j), glut_med1_mif, glut_med1_ofl, glut_med1_tsl, glut_med1_pa); %Force
    glut_med1_tqX(i,j) = glut_med1_f(i,j)*glut_med1_maX(i,j); %Torque, x axis
    glut_med1_tqY(i,j) = glut_med1_f(i,j)*glut_med1_maY(i,j); %Torque, y axis
    glut_med1_tqZ(i,j) = glut_med1_f(i,j)*glut_med1_maZ(i,j); %Torque, z axis
    % %Gluteus Medius, 2
    glut_med2_length(i,j) = norm(VecTrans(T_h,glut_med2_i) - glut_med2_o); %Total length
    glut_med2_maX(i,j) = CrossProd(glut_med2_i, (VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i)/norm(VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i),1); %Moment arm, x axis
    glut_med2_maY(i,j) = CrossProd(glut_med2_i, (VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i)/norm(VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i),2); %Moment arm, y axis
    glut_med2_maZ(i,j) = CrossProd(glut_med2_i, (VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i)/norm(VecTrans(T_h\eye(4),glut_med2_o)-glut_med2_i),3); %Moment arm, z axis
    glut_med2_f(i,j) = forz(glut_med2_length(i,j), glut_med2_mif, glut_med2_ofl, glut_med2_tsl, glut_med2_pa); %Force
    glut_med2_tqX(i,j) = glut_med2_f(i,j)*glut_med2_maX(i,j); %Torque, x axis
    glut_med2_tqY(i,j) = glut_med2_f(i,j)*glut_med2_maY(i,j); %Torque, y axis
    glut_med2_tqZ(i,j) = glut_med2_f(i,j)*glut_med2_maZ(i,j); %Torque, z axis
    % %Gluteus Medius, 3
    glut_med3_length(i,j) = norm(VecTrans(T_h,glut_med3_i) - glut_med3_o); %Total length
    glut_med3_maX(i,j) = CrossProd(glut_med3_i, (VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i)/norm(VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i),1); %Moment arm, x axis
    glut_med3_maY(i,j) = CrossProd(glut_med3_i, (VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i)/norm(VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i),2); %Moment arm, y axis
    glut_med3_maZ(i,j) = CrossProd(glut_med3_i, (VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i)/norm(VecTrans(T_h\eye(4),glut_med3_o)-glut_med3_i),3); %Moment arm, z axis
    glut_med3_f(i,j) = forz(glut_med3_length(i,j), glut_med3_mif, glut_med3_ofl, glut_med3_tsl, glut_med3_pa); %Force
    glut_med3_tqX(i,j) = glut_med3_f(i,j)*glut_med3_maX(i,j); %Torque, x axis
    glut_med3_tqY(i,j) = glut_med3_f(i,j)*glut_med3_maY(i,j); %Torque, y axis
    glut_med3_tqZ(i,j) = glut_med3_f(i,j)*glut_med3_maZ(i,j); %Torque, z axis
    % %Gluteus Minimus, 1
    glut_min1_length(i,j) = norm(VecTrans(T_h,glut_min1_i) - glut_min1_o); %Total length
    glut_min1_maX(i,j) = CrossProd(glut_min1_i, (VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i)/norm(VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i),1); %Moment arm, x axis
    glut_min1_maY(i,j) = CrossProd(glut_min1_i, (VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i)/norm(VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i),2); %Moment arm, y axis
    glut_min1_maZ(i,j) = CrossProd(glut_min1_i, (VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i)/norm(VecTrans(T_h\eye(4),glut_min1_o)-glut_min1_i),3); %Moment arm, z axis
    glut_min1_f(i,j) = forz(glut_min1_length(i,j), glut_min1_mif, glut_min1_ofl, glut_min1_tsl, glut_min1_pa); %Force
    glut_min1_tqX(i,j) = glut_min1_f(i,j)*glut_min1_maX(i,j); %Torque, x axis
    glut_min1_tqY(i,j) = glut_min1_f(i,j)*glut_min1_maY(i,j); %Torque, y axis
    glut_min1_tqZ(i,j) = glut_min1_f(i,j)*glut_min1_maZ(i,j); %Torque, z axis
    % %Gluteus Minimus, 2
    glut_min2_length(i,j) = norm(VecTrans(T_h,glut_min2_i) - glut_min2_o); %Total length
    glut_min2_maX(i,j) = CrossProd(glut_min2_i, (VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i)/norm(VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i),1); %Moment arm, x axis
    glut_min2_maY(i,j) = CrossProd(glut_min2_i, (VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i)/norm(VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i),2); %Moment arm, y axis
    glut_min2_maZ(i,j) = CrossProd(glut_min2_i, (VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i)/norm(VecTrans(T_h\eye(4),glut_min2_o)-glut_min2_i),3); %Moment arm, z axis
    glut_min2_f(i,j) = forz(glut_min2_length(i,j), glut_min2_mif, glut_min2_ofl, glut_min2_tsl, glut_min2_pa); %Force
    glut_min2_tqX(i,j) = glut_min2_f(i,j)*glut_min2_maX(i,j); %Torque, x axis
    glut_min2_tqY(i,j) = glut_min2_f(i,j)*glut_min2_maY(i,j); %Torque, y axis
    glut_min2_tqZ(i,j) = glut_min2_f(i,j)*glut_min2_maZ(i,j); %Torque, z axis
    % %Gluteus Minimus, 3
    glut_min3_length(i,j) = norm(VecTrans(T_h,glut_min3_i) - glut_min3_o); %Total length
    glut_min3_maX(i,j) = CrossProd(glut_min3_i, (VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i)/norm(VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i),1); %Moment arm, x axis
    glut_min3_maY(i,j) = CrossProd(glut_min3_i, (VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i)/norm(VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i),2); %Moment arm, y axis
    glut_min3_maZ(i,j) = CrossProd(glut_min3_i, (VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i)/norm(VecTrans(T_h\eye(4),glut_min3_o)-glut_min3_i),3); %Moment arm, z axis
    glut_min3_f(i,j) = forz(glut_min3_length(i,j), glut_min3_mif, glut_min3_ofl, glut_min3_tsl, glut_min3_pa); %Force
    glut_min3_tqX(i,j) = glut_min3_f(i,j)*glut_min3_maX(i,j); %Torque, x axis
    glut_min3_tqY(i,j) = glut_min3_f(i,j)*glut_min3_maY(i,j); %Torque, y axis
    glut_min3_tqZ(i,j) = glut_min3_f(i,j)*glut_min3_maZ(i,j); %Torque, z axis
    % %Gluteus Maximus, 1
    glut_max1_length(i,j) = norm(glut_max1_wr1-glut_max1_o)+norm(VecTrans(T_h,glut_max1_wr2) - glut_max1_wr1)+norm(glut_max1_i-glut_max1_wr2); %Total length
    glut_max1_maX(i,j) = CrossProd(glut_max1_wr2, (VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2)/norm(VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2),1); %Moment arm, x axis
    glut_max1_maY(i,j) = CrossProd(glut_max1_wr2, (VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2)/norm(VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2),2); %Moment arm, y axis
    glut_max1_maZ(i,j) = CrossProd(glut_max1_wr2, (VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2)/norm(VecTrans(T_h\eye(4),glut_max1_wr1)-glut_max1_wr2),3); %Moment arm, z axis
    glut_max1_f(i,j) = forz(glut_max1_length(i,j), glut_max1_mif, glut_max1_ofl, glut_max1_tsl, glut_max1_pa); %Force
    glut_max1_tqX(i,j) = glut_max1_f(i,j)*glut_max1_maX(i,j); %Torque, x axis
    glut_max1_tqY(i,j) = glut_max1_f(i,j)*glut_max1_maY(i,j); %Torque, y axis
    glut_max1_tqZ(i,j) = glut_max1_f(i,j)*glut_max1_maZ(i,j); %Torque, z axis
    % %Gluteus Maximus, 2
    glut_max2_length(i,j) = norm(glut_max2_wr1-glut_max2_o)+norm(VecTrans(T_h,glut_max2_wr2) - glut_max2_wr1)+norm(glut_max2_i-glut_max2_wr2); %Total length
    glut_max2_maX(i,j) = CrossProd(glut_max2_wr2, (VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2)/norm(VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2),1); %Moment arm, x axis
    glut_max2_maY(i,j) = CrossProd(glut_max2_wr2, (VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2)/norm(VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2),2); %Moment arm, y axis
    glut_max2_maZ(i,j) = CrossProd(glut_max2_wr2, (VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2)/norm(VecTrans(T_h\eye(4),glut_max2_wr1)-glut_max2_wr2),3); %Moment arm, z axis
    glut_max2_f(i,j) = forz(glut_max2_length(i,j), glut_max2_mif, glut_max2_ofl, glut_max2_tsl, glut_max2_pa); %Force
    glut_max2_tqX(i,j) = glut_max2_f(i,j)*glut_max2_maX(i,j); %Torque, x axis
    glut_max2_tqY(i,j) = glut_max2_f(i,j)*glut_max2_maY(i,j); %Torque, y axis
    glut_max2_tqZ(i,j) = glut_max2_f(i,j)*glut_max2_maZ(i,j); %Torque, z axis
    % %Gluteus Maximus, 3
    glut_max3_length(i,j) = norm(glut_max3_wr1-glut_max3_o)+norm(VecTrans(T_h,glut_max3_wr2) - glut_max3_wr1)+norm(glut_max3_i-glut_max3_wr2); %Total length
    glut_max3_maX(i,j) = CrossProd(glut_max3_wr2, (VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2)/norm(VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2),1); %Moment arm, x axis
    glut_max3_maY(i,j) = CrossProd(glut_max3_wr2, (VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2)/norm(VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2),2); %Moment arm, y axis
    glut_max3_maZ(i,j) = CrossProd(glut_max3_wr2, (VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2)/norm(VecTrans(T_h\eye(4),glut_max3_wr1)-glut_max3_wr2),3); %Moment arm, z axis
    glut_max3_f(i,j) = forz(glut_max3_length(i,j), glut_max3_mif, glut_max3_ofl, glut_max3_tsl, glut_max3_pa); %Force
    glut_max3_tqX(i,j) = glut_max3_f(i,j)*glut_max3_maX(i,j); %Torque, x axis
    glut_max3_tqY(i,j) = glut_max3_f(i,j)*glut_max3_maY(i,j); %Torque, y axis
    glut_max3_tqZ(i,j) = glut_max3_f(i,j)*glut_max3_maZ(i,j); %Torque, z axis
    % %Adductor Longus
    add_long_length(i,j) = norm(VecTrans(T_h,add_long_i) - add_long_o); %Total length, right
    add_long_maX(i,j) = CrossProd(add_long_i, (VecTrans(T_h\eye(4),add_long_o)-add_long_i)/norm(VecTrans(T_h\eye(4),add_long_o)-add_long_i),1); %Moment arm, x axis
    add_long_maY(i,j) = CrossProd(add_long_i, (VecTrans(T_h\eye(4),add_long_o)-add_long_i)/norm(VecTrans(T_h\eye(4),add_long_o)-add_long_i),2); %Moment arm, y axis
    add_long_maZ(i,j) = CrossProd(add_long_i, (VecTrans(T_h\eye(4),add_long_o)-add_long_i)/norm(VecTrans(T_h\eye(4),add_long_o)-add_long_i),3); %Moment arm, z axis
    add_long_f(i,j) = forz(add_long_length(i,j), add_long_mif, add_long_ofl, add_long_tsl, add_long_pa); %Force
    add_long_tqX(i,j) = add_long_f(i,j)*add_long_maX(i,j); %Torque, x axis
    add_long_tqY(i,j) = add_long_f(i,j)*add_long_maY(i,j); %Torque, y axis
    add_long_tqZ(i,j) = add_long_f(i,j)*add_long_maZ(i,j); %Torque, z axis
    % %Adductor Brevis
    add_brev_length(i,j) = norm(VecTrans(T_h,add_brev_i) - add_brev_o); %Total length
    add_brev_maX(i,j) = CrossProd(add_brev_i, (VecTrans(T_h\eye(4),add_brev_o)-add_brev_i)/norm(VecTrans(T_h\eye(4),add_brev_o)-add_brev_i),1); %Moment arm, x axis
    add_brev_maY(i,j) = CrossProd(add_brev_i, (VecTrans(T_h\eye(4),add_brev_o)-add_brev_i)/norm(VecTrans(T_h\eye(4),add_brev_o)-add_brev_i),2); %Moment arm, y axis
    add_brev_maZ(i,j) = CrossProd(add_brev_i, (VecTrans(T_h\eye(4),add_brev_o)-add_brev_i)/norm(VecTrans(T_h\eye(4),add_brev_o)-add_brev_i),3); %Moment arm, z axis
    add_brev_f(i,j) = forz(add_brev_length(i,j), add_brev_mif, add_brev_ofl, add_brev_tsl, add_brev_pa); %Force
    add_brev_tqX(i,j) = add_brev_f(i,j)*add_brev_maX(i,j); %Torque, x axis
    add_brev_tqY(i,j) = add_brev_f(i,j)*add_brev_maY(i,j); %Torque, y axis
    add_brev_tqZ(i,j) = add_brev_f(i,j)*add_brev_maZ(i,j); %Torque, z axis
    % %Adductor Magnus, 1
    add_mag1_length(i,j) = norm(VecTrans(T_h,add_mag1_i) - add_mag1_o); %Total length
    add_mag1_maX(i,j) = CrossProd(add_mag1_i, (VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i)/norm(VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i),1); %Moment arm, x axis
    add_mag1_maY(i,j) = CrossProd(add_mag1_i, (VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i)/norm(VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i),2); %Moment arm, y axis
    add_mag1_maZ(i,j) = CrossProd(add_mag1_i, (VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i)/norm(VecTrans(T_h\eye(4),add_mag1_o)-add_mag1_i),3); %Moment arm, z axis
    add_mag1_f(i,j) = forz(add_mag1_length(i,j), add_mag1_mif, add_mag1_ofl, add_mag1_tsl, add_mag1_pa); %Force
    add_mag1_tqX(i,j) = add_mag1_f(i,j)*add_mag1_maX(i,j); %Torque, x axis
    add_mag1_tqY(i,j) = add_mag1_f(i,j)*add_mag1_maY(i,j); %Torque, y axis
    add_mag1_tqZ(i,j) = add_mag1_f(i,j)*add_mag1_maZ(i,j); %Torque, z axis
    % %Adductor Magnus, 2
    add_mag2_length(i,j) = norm(VecTrans(T_h,add_mag2_i) - add_mag2_o); %Total length
    add_mag2_maX(i,j) = CrossProd(add_mag2_i, (VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i)/norm(VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i),1); %Moment arm, x axis
    add_mag2_maY(i,j) = CrossProd(add_mag2_i, (VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i)/norm(VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i),2); %Moment arm, y axis
    add_mag2_maZ(i,j) = CrossProd(add_mag2_i, (VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i)/norm(VecTrans(T_h\eye(4),add_mag2_o)-add_mag2_i),3); %Moment arm, z axis
    add_mag2_f(i,j) = forz(add_mag2_length(i,j), add_mag2_mif, add_mag2_ofl, add_mag2_tsl, add_mag2_pa); %Force
    add_mag2_tqX(i,j) = add_mag2_f(i,j)*add_mag2_maX(i,j); %Torque, x axis
    add_mag2_tqY(i,j) = add_mag2_f(i,j)*add_mag2_maY(i,j); %Torque, y axis
    add_mag2_tqZ(i,j) = add_mag2_f(i,j)*add_mag2_maZ(i,j); %Torque, z axis
    % %Adductor Magnus, 3
    add_mag3_length(i,j) = norm(VecTrans(T_h,add_mag3_i) - add_mag3_o); %Total length
    add_mag3_maX(i,j) = CrossProd(add_mag3_i, (VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i)/norm(VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i),1); %Moment arm, x axis
    add_mag3_maY(i,j) = CrossProd(add_mag3_i, (VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i)/norm(VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i),2); %Moment arm, y axis
    add_mag3_maZ(i,j) = CrossProd(add_mag3_i, (VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i)/norm(VecTrans(T_h\eye(4),add_mag3_o)-add_mag3_i),3); %Moment arm, z axis
    add_mag3_f(i,j) = forz(add_mag3_length(i,j), add_mag3_mif, add_mag3_ofl, add_mag3_tsl, add_mag3_pa); %Force
    add_mag3_tqX(i,j) = add_mag3_f(i,j)*add_mag3_maX(i,j); %Torque, x axis
    add_mag3_tqY(i,j) = add_mag3_f(i,j)*add_mag3_maY(i,j); %Torque, y axis
    add_mag3_tqZ(i,j) = add_mag3_f(i,j)*add_mag3_maZ(i,j); %Torque, z axis
    %Pectineus
    pect_length(i,j) = norm(VecTrans(T_h,pect_i) - pect_o); %Total length
    pect_maX(i,j) = CrossProd(pect_i, (VecTrans(T_h\eye(4),pect_o)-pect_i)/norm(VecTrans(T_h\eye(4),pect_o)-pect_i),1); %Moment arm, x axis
    pect_maY(i,j) = CrossProd(pect_i, (VecTrans(T_h\eye(4),pect_o)-pect_i)/norm(VecTrans(T_h\eye(4),pect_o)-pect_i),2); %Moment arm, y axis
    pect_maZ(i,j) = CrossProd(pect_i, (VecTrans(T_h\eye(4),pect_o)-pect_i)/norm(VecTrans(T_h\eye(4),pect_o)-pect_i),3); %Moment arm, z axis
    pect_f(i,j) = forz(pect_length(i,j), pect_mif, pect_ofl, pect_tsl, pect_pa); %Force
    pect_tqX(i,j) = pect_f(i,j)*pect_maX(i,j); %Torque, x axis
    pect_tqY(i,j) = pect_f(i,j)*pect_maX(i,j); %Torque, y axis
    pect_tqZ(i,j) = pect_f(i,j)*pect_maX(i,j); %Torque, z axis
    %Iliacus
    iliacus_length(i,j) = norm(iliacus_wr1-iliacus_o)+norm(VecTrans(T_h,iliacus_wr3)-iliacus_wr1)+norm(iliacus_wr3-iliacus_i); %Total length
    iliacus_maX(i,j) = CrossProd(iliacus_wr3, (VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3),1); %Moment arm, x axis
    iliacus_maY(i,j) = CrossProd(iliacus_wr3, (VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3),2); %Moment arm, y axis
    iliacus_maZ(i,j) = CrossProd(iliacus_wr3, (VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3)/norm(VecTrans(T_h\eye(4),iliacus_wr1)-iliacus_wr3),3); %Moment arm, z axis
    iliacus_f(i,j) = forz(iliacus_length(i,j), iliacus_mif, iliacus_ofl, iliacus_tsl, iliacus_pa); %Force
    iliacus_tqX(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, x axis
    iliacus_tqY(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, y axis
    iliacus_tqZ(i,j) = iliacus_f(i,j)*iliacus_maX(i,j); %Torque, z axis
    %Psoas
    psoas_length(i,j) = norm(psoas_wr1-psoas_o)+norm(VecTrans(T_h,psoas_wr3)-psoas_wr1)+norm(psoas_wr3-psoas_i); %Total length, right
    psoas_maHX(i,j) = CrossProd(psoas_wr3, (VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3)/norm(VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3),1); %Moment arm, x axis
    psoas_maHY(i,j) = CrossProd(psoas_wr3, (VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3)/norm(VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3),2); %Moment arm, y axis
    psoas_maHZ(i,j) = CrossProd(psoas_wr3, (VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3)/norm(VecTrans(T_h\eye(4),psoas_wr1)-psoas_wr3),3); %Moment arm, z axis
    psoas_f(i,j) = forz(psoas_length(i,j), psoas_mif, psoas_ofl, psoas_tsl, psoas_pa); %Force
    psoas_tqHX(i,j) = psoas_f(i,j)*psoas_maHX(i,j); %Torque, x axis
    psoas_tqHY(i,j) = psoas_f(i,j)*psoas_maHY(i,j); %Torque, y axis
    psoas_tqHZ(i,j) = psoas_f(i,j)*psoas_maHZ(i,j); %Torque, z axis
    %Quadriceps Femoris
    quad_fem_length(i,j) = norm(VecTrans(T_h,quad_fem_i) - quad_fem_o); %Total length
    quad_fem_maX(i,j) = CrossProd(quad_fem_i, (VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i)/norm(VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i),1); %Moment arm, x axis
    quad_fem_maY(i,j) = CrossProd(quad_fem_i, (VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i)/norm(VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i),2); %Moment arm, y axis
    quad_fem_maZ(i,j) = CrossProd(quad_fem_i, (VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i)/norm(VecTrans(T_h\eye(4),quad_fem_o)-quad_fem_i),3); %Moment arm, z axis
    quad_fem_f(i,j) = forz(quad_fem_length(i,j), quad_fem_mif, quad_fem_ofl, quad_fem_tsl, quad_fem_pa); %Force
    quad_fem_tqX(i,j) = quad_fem_f(i,j)*quad_fem_maX(i,j); %Torque, x axis
    quad_fem_tqY(i,j) = quad_fem_f(i,j)*quad_fem_maY(i,j); %Torque, x axis
    quad_fem_tqZ(i,j) = quad_fem_f(i,j)*quad_fem_maZ(i,j); %Torque, x axis
    %Gemelli
    gem_length(i,j) = norm(VecTrans(T_h,gem_i) - gem_o); %Total length
    gem_maX(i,j) = CrossProd(gem_i, (VecTrans(T_h\eye(4),gem_o)-gem_i)/norm(VecTrans(T_h\eye(4),gem_o)-gem_i),1); %Moment arm, x axis
    gem_maY(i,j) = CrossProd(gem_i, (VecTrans(T_h\eye(4),gem_o)-gem_i)/norm(VecTrans(T_h\eye(4),gem_o)-gem_i),2); %Moment arm, y axis
    gem_maZ(i,j) = CrossProd(gem_i, (VecTrans(T_h\eye(4),gem_o)-gem_i)/norm(VecTrans(T_h\eye(4),gem_o)-gem_i),3); %Moment arm, z axis
    gem_f(i,j) = forz(gem_length(i,j), gem_mif, gem_ofl, gem_tsl, gem_pa); %Force
    gem_tqX(i,j) = gem_f(i,j)*gem_maX(i,j); %Torque, x axis
    gem_tqY(i,j) = gem_f(i,j)*gem_maX(i,j); %Torque, y axis
    gem_tqZ(i,j) = gem_f(i,j)*gem_maX(i,j); %Torque, z axis
    %Piriformis
    peri_length(i,j) = norm(VecTrans(T_h,peri_i)-peri_wr1)+norm(peri_o-peri_wr1); %Total length
    peri_maX(i,j) = CrossProd(peri_i, (VecTrans(T_h\eye(4),peri_wr1)-peri_i)/norm(VecTrans(T_h\eye(4),peri_wr1)-peri_i),1); %Moment arm, x axis
    peri_maY(i,j) = CrossProd(peri_i, (VecTrans(T_h\eye(4),peri_wr1)-peri_i)/norm(VecTrans(T_h\eye(4),peri_wr1)-peri_i),2); %Moment arm, y axis
    peri_maZ(i,j) = CrossProd(peri_i, (VecTrans(T_h\eye(4),peri_wr1)-peri_i)/norm(VecTrans(T_h\eye(4),peri_wr1)-peri_i),3); %Moment arm, z axis
    peri_f(i,j) = forz(peri_length(i,j), peri_mif, peri_ofl, peri_tsl, peri_pa); %Force
    peri_tqX(i,j) = peri_f(i,j)*peri_maX(i,j); %Torque, x axis
    peri_tqY(i,j) = peri_f(i,j)*peri_maY(i,j); %Torque, y axis
    peri_tqZ(i,j) = peri_f(i,j)*peri_maZ(i,j); %Torque, z axis
%     %Rectus Femoris
%     rect_fem_o = [-0.03; -0.031; 0.097]; %Pelvis, right, origin
%     rect_fem_i = [fcn3(theta_k); fcn4(theta_k); 0.001]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% 
%     %Vastus Medialis
%     vas_med_o = [0.014; -0.21; 0.019]; %Femur, right, origin
%     vas_med_wr1 = [0.036; -0.277; 0.001]; %Femur, right, wrapping point 1
%     vas_med_i = [fcn5(theta_k); fcn6(theta_k); -0.015]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% 
%     %Vastus Intermedius
%     vas_int_o = [0.029; -0.192; 0.031]; %Femur, right, origin
%     vas_int_wr1 = [0.034; -0.208; 0.028]; %Femur, right, wrapping point 1
%     vas_int_i = [fcn7(theta_k); fcn8(theta_k); 0.002]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)
% 
%     %Vastus Lateralis
%     vas_lat_o = [0.005; -0.185; 0.035]; %Femur, right, origin
%     vas_lat_wr1 = [0.027; -0.259; 0.041]; %Femur, right, wrapping point 1
% 
%     vas_lat_i = [fcn9(theta_k); fcn10(theta_k); 0.016]; %Tibia, right, insertion (left in, not a stationary point in OpenSim due to patella)

    end
    end


    uni_hip_tqZ = glut_min1_tqZ+glut_min2_tqZ+glut_min3_tqZ+glut_med1_tqZ+glut_med2_tqZ+glut_med3_tqZ+glut_max1_tqZ+glut_max2_tqZ+glut_max3_tqZ+add_long_tqZ+add_brev_tqZ+add_mag1_tqZ+add_mag2_tqZ+add_mag3_tqZ+pect_tqZ+iliacus_tqZ+psoas_tqHZ+quad_fem_tqZ+gem_tqZ+peri_tqZ;
    uni_hip_tqY = glut_min1_tqY+glut_min2_tqY+glut_min3_tqY+glut_med1_tqY+glut_med2_tqY+glut_med3_tqY+glut_max1_tqY+glut_max2_tqY+glut_max3_tqY+add_long_tqY+add_brev_tqY+add_mag1_tqY+add_mag2_tqY+add_mag3_tqY+pect_tqY+iliacus_tqY+psoas_tqHY+quad_fem_tqY+gem_tqY+peri_tqY;
    uni_hip_tqX = glut_min1_tqX+glut_min2_tqX+glut_min3_tqX+glut_med1_tqX+glut_med2_tqX+glut_med3_tqX+glut_max1_tqX+glut_max2_tqX+glut_max3_tqX+add_long_tqX+add_brev_tqX+add_mag1_tqX+add_mag2_tqX+add_mag3_tqX+pect_tqX+iliacus_tqX+psoas_tqHX+quad_fem_tqX+gem_tqX+peri_tqX;

    %Including a way to make the data generic, so that things can be
    %plotted in one location instead of spread out through all of the
    %different sections
    HumanAxis1 = theta_hZ;
    HumanAxis1Label = 'Hip Flexion, Degrees';
    HumanAxis2 = theta_hX;
    HumanAxis2Label = 'Adduction/Abduction, Degrees';
    HumanTorque1 = real(uni_hip_tqZ);
    HumanTorque2 = real(uni_hip_tqY);
    HumanTorque3 = real(uni_hip_tqX);
    HumanTitle1 = 'Hip Torque, uniarticular, Z axis';
    HumanTitle2 = 'Hip Torque, uniarticular, Y axis';
    HumanTitle3 = 'Hip Torque, uniarticular, X axis';
    
    save(strcat('Human_', ChooseJoint, '_Data.mat'), 'HumanAxis1', 'HumanAxis1Label', 'HumanAxis2', 'HumanAxis2Label', 'HumanTorque1', 'HumanTorque2', 'HumanTorque3', 'HumanTitle1', 'HumanTitle2', 'HumanTitle3')
else
    fprintf("Incorrect Joint Selected. Check spelling\n")
end

caxisRange = [-40 150];

figure
hold on
surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque1, 'EdgeColor', 'none')
title(HumanTitle1); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label); zlabel('Torque, N*m')
colorbar; caxis(caxisRange)
hold off

figure
surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque2, 'EdgeColor', 'none')
title(HumanTitle2); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label); zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('HumanTorque3', 'var') == 1
    figure
    surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque3, 'EdgeColor', 'none')
    title(HumanTitle3); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label); zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('HumanTorque4', 'var') == 1
    figure
    surf(HumanAxis1*180/pi, HumanAxis2*180/pi, HumanTorque4, 'EdgeColor', 'none')
    title(HumanTitle4); xlabel(HumanAxis1Label); ylabel(HumanAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end
