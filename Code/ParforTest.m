clear
close all
clc

addpath('Functions');

run('RobotPAMCalculationOptimization.m');
 

beginOptimization = 1;

ep = 0.1;

ID = 0;
size = 0;

parfor iiii = 1:10
    [ID, size] = Size(Muscle1.LongestL, Muscle1.MIF);
end


    