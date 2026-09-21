%% export_final_20260920.m — build KneeReflexCircuit view + final exports
%% (EPS chain now uses an in-folder temp crop file).

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
cd(sns); addpath(sns); addpath(fullfile(sns, 'demos'));
fprintf('=== MATLAB %s ===\n', version);

%% 1. build the dissertation circuit view
run(fullfile(sns, 'demos', 'sns_build_circuit_view.m'));

%% 2. export it + KneeReflexDemo (full model, for reference) with EPS
sns_export_diagram(fullfile(sns, 'figures'), {'KneeReflexCircuit', 'KneeReflexDemo'}, false);

fprintf('=== EXPORT FINAL DONE ===\n');
