%rerun_FlxBiomimetic.m — refresh the biomimetic FLEXOR data chain per Ben
%2026-10-19: run Knee_Flexor_Data_20mm + Knee_Flexor_Data (old mats
%archived), save KneeFlx_20mm_42cm.mat / KneeFlx_10mm_42cm.mat, then run
%the Plot functions in Testing_Data. Prints the Xi each script actually
%uses so the terms of record can be checked.
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Human_Data'));
addpath(genpath(fullfile(root, 'Code/Matlab/Mesh_Optimization')));
warning('off', 'all');
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
resdir = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results');

%% 20 mm flexor data -> KneeFlx_20mm_42cm.mat
cd(resdir);
m20 = fullfile(festodir, 'KneeFlx_20mm_42cm.mat');
a20 = fullfile(festodir, 'KneeFlx_20mm_42cm_stale_20261019.mat');
if exist(m20, 'file') && ~exist(a20, 'file')
    copyfile(m20, a20); fprintf('ARCHIVED %s\n', a20);
end
run('Knee_Flexor_Data_20mm.m');
m20 = fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo', 'KneeFlx_20mm_42cm.mat');
save(m20);
fprintf('WROTE %s\n', m20);
close all;

%% 10 mm flexor data -> KneeFlx_10mm_42cm.mat
cd(fullfile('D:/GitHub/Bipedal_Robot', 'Code/Matlab/Mesh_Optimization/Results'));
m10 = fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo', 'KneeFlx_10mm_42cm.mat');
a10 = fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo', 'KneeFlx_10mm_42cm_stale_20261019.mat');
if exist(m10, 'file') && ~exist(a10, 'file')
    copyfile(m10, a10); fprintf('ARCHIVED %s\n', a10);
end
run('Knee_Flexor_Data.m');
m10 = fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo', 'KneeFlx_10mm_42cm.mat');
save(m10);
fprintf('WROTE %s\n', m10);
close all;

%% Plot functions (Testing_Data section)
cd(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo'));
addpath(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo/Flx_10mm'));
addpath(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo/Flx_20mm'));
try
    run('Plot_KneeFlx_20mm_42cm.m'); fprintf('PLOT 20mm OK\n');
catch e
    fprintf('PLOT 20mm FAILED: %s\n', e.message);
end
try
    run('Plot_KneeFlx_10mm_42cm.m'); fprintf('PLOT 10mm OK\n');
catch e
    fprintf('PLOT 10mm FAILED: %s\n', e.message);
end
rmpath(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo/Flx_10mm'));
rmpath(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo/Flx_20mm'));
close all;
disp('DONE_FlxBiomimetic');
