%% figs_and_exports2_20260920.m — exports round 2 (pdfcrop+mgs EPS chain),
%% circuit redraw, Deng re-verify, R2025a library copy.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
cd(sns); addpath(sns); addpath(fullfile(sns, 'demos'));
fprintf('=== MATLAB %s ===\n', version);

%% 1. dissertation exports
sns_export_diagram(fullfile(sns, 'figures'), ...
    {'KneeReflexDemo', 'SNS_Deng_RG', 'BeerCupReflexDemo'}, false);

%% 2. redraw the vector circuit figure at 16.5 cm / 10-11 pt
run(fullfile(sns, 'sns_draw_circuit.m'));

%% 3. re-verify Deng demo (base-MATLAB Pearson)
sns_run_deng_demo();
fprintf('DENG RE-VERIFY DONE\n');

%% 4. R2025a copy of the rebuilt library
Simulink.exportToVersion('SNS_Library', 'SNS_Library_R2025a', 'R2025a');
fprintf('EXPORTED SNS_Library_R2025a.slx\n');
close_system('SNS_Library', 0);
close_system('SNS_Library_R2025a', 0);
fprintf('=== FIGS2 DONE ===\n');
