%% Example
clear;
clc;
close all;


k1 = 100;
F = linspace(0,500);
x1 = F./k1;

k2 = 1000;
xT = F./k1+F./k2;

figure
hold on
plot(F,x1,'DisplayName','spring 2 infinite')
plot(F,xT,'DisplayName','spring 2 >> spring 1')
hold off
xlabel('Force')
ylabel('Distance')
lgd = legend;

figure
hold on
plot(x1,F,'DisplayName','spring 2 infinite')
plot(xT,F,'DisplayName','spring 2 >> spring 1')
hold off
xlabel('Distance')
ylabel('Force')
lgd = legend;