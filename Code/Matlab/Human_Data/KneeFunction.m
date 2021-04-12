clear
close all
clc

%% Open Sim model

knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');

for i = 1:length(theta)
    xo(i) = fcn1(theta(i));
    yo(i) = fcn2(theta(i));
end

figure
hold on
subplot(2, 1, 1)
plot(theta, xo, knee_angle_x, knee_x, 'o')

subplot(2, 1, 2)
plot(theta, yo, knee_angle_y, knee_y, 'o')
hold off