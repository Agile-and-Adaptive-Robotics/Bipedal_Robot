clear
%close all
clc

%% Open Sim model

kneeMin = -2.0943951;
kneeMax = 0.17453293;
theta = linspace(kneeMin, kneeMax, 100);

%Human
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');

%Robot
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     [0.0010	0.0027	0.0038	0.0045	0.0064	0.0084	0.0105	0.0164	0.0213	0.0246	0.0255	0.0239	0.0197	0.0132	0.0052	-0.0036]';
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     [-0.3982	-0.3969	-0.3962	-0.3958	-0.3950	-0.3944	-0.3942	-0.3951	-0.3984	-0.4035	-0.4099	-0.4167	-0.4228	-0.4274	-0.4298	-0.4292]';
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

for i = 1:length(theta)
    xo(i) = fcn1(theta(i));
    yo(i) = fcn2(theta(i));
    x1(i) = fcn3(theta(i));
    y1(i) = fcn4(theta(i));
end

figure    %Plot human knee
hold on
subplot(2, 1, 1)
plot(theta, xo, knee_angle_x, knee_x, 'o')

subplot(2, 1, 2)
plot(theta, yo, knee_angle_y, knee_y, 'o')
hold off

figure    %Plot robot knee
hold on
subplot(2, 1, 1)
plot(theta, x1, knee_angle, knee_x_Pam, 'o')

subplot(2, 1, 2)
plot(theta, y1, knee_angle, knee_y_Pam, 'o')
hold off