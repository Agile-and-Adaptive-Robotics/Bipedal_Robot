clc                 %Command Line Clear
clear               %Clear the workspace of stored variables
close all           %Close all open figures

% x = [0, 0.1, 0.17, 0.25]';
% y = [630, 300, 150, 0]';
% 
% figure
% plot(x, y)
% 
% 
% tenFit = fit(x, y, 'exp2');
% 
% xFit = linspace(0, 0.25);
% yFit = tenFit(xFit);
% 
% figure
% hold on
% plot(xc, yc, '.')
% plot(x, y, 'r.')
% plot(xFit, yFit)
% hold off

% x = [0, 0.07, 0.11, 0.15, 0.25]';
% y = [1400, 800, 600, 400, 0]';
% 
% figure
% plot(x, y)
% 
% 
% twentyFit = fit(x, y, 'poly2');
% 
% xFit = linspace(0, 0.25);
% yFit = twentyFit(xFit);
% 
% figure
% hold on
% plot(x, y, 'r.')
% plot(xFit, yFit)
% hold off

x = [0, 0.06, 0.12, 0.15, 0.25]';
y = [6000, 3500, 2000, 1500, 0]';

figure
plot(x, y)


twentyFit = fit(x, y, 'poly2');

xFit = linspace(0, 0.25);
yFit = twentyFit(xFit);

figure
hold on
plot(x, y, 'r.')
plot(xFit, yFit)
hold off