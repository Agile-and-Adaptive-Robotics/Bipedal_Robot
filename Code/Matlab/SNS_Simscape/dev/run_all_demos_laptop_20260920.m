%% run_all_demos_laptop_20260920.m — verify all SNS demo models run on R2025b laptop
% Runs the four demo runners in sequence and prints their summary lines.
% Log: logs\run_all_demos_laptop_20260920.log

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));  % dev -> SNS_Simscape -> Matlab
sns  = fullfile(root, 'SNS_Simscape');
addpath(sns); addpath(fullfile(sns, 'demos'));

fprintf('=== MATLAB %s on %s ===\n', version, computer);

try
    % 1. KneeReflexDemo
    run(fullfile(sns, 'demos', 'sns_run_demo.m'));
    fprintf('DEMO 1 KneeReflexDemo: OK\n');
catch ME
    fprintf('DEMO 1 KneeReflexDemo: FAILED — %s\n', ME.message);
end

try
    % 2. BPACPGLegDemo
    run(fullfile(sns, 'demos', 'sns_run_cpg_demo.m'));
    fprintf('DEMO 2 BPACPGLegDemo: OK\n');
catch ME
    fprintf('DEMO 2 BPACPGLegDemo: FAILED — %s\n', ME.message);
end

try
    % 3. BeerCupReflexDemo
    run(fullfile(sns, 'demos', 'sns_run_beer_demo.m'));
    fprintf('DEMO 3 BeerCupReflexDemo: OK\n');
catch ME
    fprintf('DEMO 3 BeerCupReflexDemo: FAILED — %s\n', ME.message);
end

try
    % 4. Deng CPG demo
    sns_run_deng_demo();
    fprintf('DEMO 4 SNS_Deng_CPGDemo: OK\n');
catch ME
    fprintf('DEMO 4 SNS_Deng_CPGDemo: FAILED — %s\n', ME.message);
end

fprintf('=== run_all_demos DONE ===\n');
