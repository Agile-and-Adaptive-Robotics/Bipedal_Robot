%% figs_and_exports_20260920.m — icon preview + dissertation exports + Deng
%% re-verify, after the Rybak-style library rebuild + font-param fix.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
cd(sns); addpath(sns); addpath(fullfile(sns, 'demos'));
fprintf('=== MATLAB %s ===\n', version);

%% 1. icon preview: neuron -> E-syn / I-syn -> neuron chain + Ia + BPA_20mm
mdl = 'sns_icon_preview';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); open_system(mdl);
add_block('SNS_Library/NonSpikingNeuron', [mdl '/pre_E'],  'Position', [60 60 130 130]);
add_block('SNS_Library/NonSpikingSynapse', [mdl '/synE'], 'Position', [200 60 250 110], ...
    'Esyn', '0', 'gmax', '1');
add_block('SNS_Library/NonSpikingNeuron', [mdl '/post_E'], 'Position', [330 60 400 130]);
add_block('SNS_Library/NonSpikingNeuron', [mdl '/pre_I'],  'Position', [60 220 130 290]);
add_block('SNS_Library/NonSpikingSynapse', [mdl '/synI'], 'Position', [200 220 250 270], ...
    'Esyn', '-72', 'gmax', '1');
add_block('SNS_Library/NonSpikingNeuron', [mdl '/post_I'], 'Position', [330 220 400 290]);
add_block('SNS_Library/IaMuscleSpindle', [mdl '/Ia'], 'Position', [60 380 130 450]);
add_block('SNS_Library/BPA_20mm', [mdl '/BPA'], 'Position', [330 380 400 450]);
add_line(mdl, 'pre_E/1',  'synE/1', 'autorouting', 'on');
add_line(mdl, 'post_E/1', 'synE/2', 'autorouting', 'on');
add_line(mdl, 'synE/1',   'post_E/1', 'autorouting', 'on');
add_line(mdl, 'pre_I/1',  'synI/1', 'autorouting', 'on');
add_line(mdl, 'post_I/1', 'synI/2', 'autorouting', 'on');
add_line(mdl, 'synI/1',   'post_I/1', 'autorouting', 'on');
set_param(mdl, 'ScreenColor', 'white');
blks = find_system(mdl, 'Type', 'Block');
for b = 1:numel(blks)
    try, set_param(blks{b}, 'FontSize', '12'); catch, end
end
print(['-s' mdl], '-dpng', '-r200', fullfile(sns, 'figures', 'sns_icon_preview.png'));
fprintf('ICON PREVIEW WRITTEN\n');
close_system(mdl, 0);

%% 2. dissertation exports (eps/pdf/png, 12 pt, synapse names hidden)
sns_export_diagram(fullfile(sns, 'figures'), {'KneeReflexDemo', 'SNS_Deng_RG'}, false);

%% 3. redraw the vector circuit figure at 16.5 cm / 10-11 pt
run(fullfile(sns, 'sns_draw_circuit.m'));

%% 4. re-verify Deng demo (corr -> base-MATLAB Pearson)
sns_run_deng_demo();
fprintf('DENG RE-VERIFY DONE\n');

%% 5. R2025a copy of the rebuilt library
Simulink.exportToVersion('SNS_Library', 'SNS_Library_R2025a', 'R2025a');
fprintf('EXPORTED SNS_Library_R2025a.slx\n');
close_system('SNS_Library', 0);
close_system('SNS_Library_R2025a', 0);
fprintf('=== FIGS DONE ===\n');
