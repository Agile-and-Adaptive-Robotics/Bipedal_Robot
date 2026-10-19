%dump_ExtAll.m — ground-truth dump of all nine pinned-extensor evaluator tests:
% actual rest/tendon lengths + GoF (zero-correction and adopted row-32 terms).
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(fullfile(root, 'Testing_Data/2022_02_Festo'));
warning('off', 'all');
g = [-6.378022e-03, 3.997873e+04, 1.473375e+04, 1.582887e-01];
tests = 1:9;
[f0, b0] = minimizeExtX3(0, Inf, Inf, 0, tests);
[f1, b1] = minimizeExtX3(g(1), g(2), g(3), g(4), tests);
fprintf('\nidx  rest(m)  tendon(m)  Kmax   unopt(RMSE/FVU/Max)      opt(RMSE/FVU/Max)\n');
for i = 1:9
    fprintf('%2d   %.4f   %.4f   %.4f   %8.4f %8.4f %8.4f   %8.4f %8.4f %8.4f\n', ...
        i, b0(i).rest, b0(i).ten, b0(i).Kmax, f0(i,1), f0(i,2), f0(i,3), f1(i,1), f1(i,2), f1(i,3));
end
disp('DONE_dumpExtAll');
