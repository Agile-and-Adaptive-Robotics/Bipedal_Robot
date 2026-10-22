
clear, close('all'), clc

load('Test.mat'), load('Test2.mat')

figure, hold on, rotate3d on
h1 = surf(theta_hZ, theta_hX, uni_hip_tqZ, 'Edgecolor', 'none', 'FaceColor', 'r', 'Facealpha', 0.6);
h2 = surf(xs, ys, zs-200, 'Edgecolor', 'none', 'FaceColor', 'b', 'Facealpha', 0.6);
title('Hip Joint RoM')
ylabel('Adduction/Abduction, Degrees')
xlabel('Flexion/Extension, Degrees')
view(-90, 90)
legend('Robot', 'OpenSim')