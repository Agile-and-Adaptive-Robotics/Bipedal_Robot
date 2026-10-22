%Fit for moving attachment points and knee coordinate system
%Data from OpenSim 4.3
%% Setup
clear, clc

c = pi/180; %Convert from degrees to radians
d = 180/pi; %Convert from radians to degrees

kneeMin = -2.0943951;                           %Knee max flexion
kneeMax = 0.17453293;                           %Knee max extension
theta = linspace(kneeMin, kneeMax, 100);        %Knee Range of Motion

%% Insertion points and functions
%Knee (OpenSim data)
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');

%Rectus Femoris 
%(Insertion point as a function of knee angle due to patella)
rect_fem_x = [0.0156367; 0.0179948; 0.0274274; 0.029683; 0.0306; 0.0366; 0.0422; 0.0451; 0.0484; 0.0533; 0.0617; 0.0634; 0.067; 0.0733];
rect_fem_xD = c*[-120.118; -114.871; -90.068; -83.532; -80; -60; -40; -30; -20; -10; 0; 1.6; 5; 10];
fcn3 = fit(rect_fem_xD,rect_fem_x,'smoothingspline');
rect_fem_y =[0.0234; 0.0238; 0.0251; 0.0253; 0.025284; 0.0249; 0.0243; 0.0239; 0.0234; 0.0228; 0.0210; 0.0206; 0.0192; 0.0160];
rect_fem_yD = c*[-120; -114.6; -90; -83.5; -80.01; -60; -40; -30; -20; -10; 0; 1.6; 5; 10 ];
fcn4 = fit(rect_fem_yD,rect_fem_y,'smoothingspline');

%Vastus Medialis
%(Insertion point as a function of knee angle due to patella)
fcn5 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.009811),'smoothingspline');
fcn6 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.02242),'smoothingspline');

%Vastus Intermedius
%(Insertion point as a function of knee angle due to patella)
fcn7 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.008285),'smoothingspline');
fcn8 = fit(rect_fem_yD,rect_fem_y+(0.0256239-0.02346),'smoothingspline');

%Vastus Lateralis
%(Insertion point as a function of knee angle due to patella)
fcn9 = fit(rect_fem_xD,rect_fem_x-(0.0156367-0.0142881),'smoothingspline');
fcn10 = fit(rect_fem_yD,rect_fem_y-(0.02346-0.0215281),'smoothingspline');

%% Plotting

for i = 1:length(theta)
    xo(i) = fcn1(theta(i));  %Knee
    yo(i) = fcn2(theta(i));
    x1(i) = fcn3(theta(i));  %Rectus Femoris
    y1(i) = fcn4(theta(i));
    x2(i) = fcn5(theta(i));  %Vastus Medialis
    y2(i) = fcn6(theta(i));
    x3(i) = fcn7(theta(i));  %Vastus Intermedius
    y3(i) = fcn8(theta(i));
    x4(i) = fcn9(theta(i));  %Vastus Lateralis
    y4(i) = fcn10(theta(i));
end

figure    %Plot human knee
t0 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(t0, 'Knee Position vs. Knee Angle')
xlabel(t0, 'Knee angle, Rad')
ylabel(t0, 'Distance (mm)')
nexttile
plot(theta, xo, knee_angle_x, knee_x, 'o')
ylabel('X position')
nexttile
plot(theta, yo, knee_angle_y, knee_y, 'o')
ylabel('Y position')

figure    %Plot Rectus Femoris
t1 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(t1, 'Rectus Femoris Insertion vs. Knee Angle')
xlabel(t1, 'Knee angle, Rad')
ylabel(t1, 'Distance (mm)')
nexttile
plot(theta*d, x1, rect_fem_xD*d,rect_fem_x, 'o')
ylabel('X position')
nexttile
plot(theta*d, y1, rect_fem_yD*d,rect_fem_y, 'o')
ylabel('Y position')

figure    %Plot Vastus Medialis
t2 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(t2, 'Vastus Medialis Insertion vs. Knee Angle')
xlabel(t2, 'Knee angle, Rad')
ylabel(t2, 'Distance (mm)')
nexttile
plot(theta*d, x2, rect_fem_xD*d,rect_fem_x-(0.0156367-0.009811), 'o')
ylabel('X position')
nexttile
plot(theta*d, y2, rect_fem_yD*d,rect_fem_y-(0.02346-0.02242), 'o')
ylabel('Y position')

figure    %Plot Vastus Intermedius
t3 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(t3, 'Vastus Intermedius Insertion vs. Knee Angle')
xlabel(t3, 'Knee angle, Rad')
ylabel(t3, 'Distance (mm)')
nexttile
plot(theta*d, x3, rect_fem_xD*d,rect_fem_x-(0.0156367-0.008285), 'o')
ylabel('X position')
nexttile
plot(theta*d, y3, rect_fem_yD*d,rect_fem_y+(0.0256239-0.02346), 'o')
ylabel('Y position')

figure    %Plot Vastus Lateralis
t4 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(t4, 'Vastus Lateralis Insertion vs. Knee Angle')
xlabel(t4, 'Knee angle, Rad')
ylabel(t4, 'Distance (mm)')
nexttile
plot(theta*d, x4, rect_fem_xD*d,rect_fem_x-(0.0156367-0.0142881), 'o')
ylabel('X position')
nexttile
plot(theta*d, y4, rect_fem_yD*d,rect_fem_y-(0.02346-0.0215281), 'o')
ylabel('Y position')

save MovePoints.mat