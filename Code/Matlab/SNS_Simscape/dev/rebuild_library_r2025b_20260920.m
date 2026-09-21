%% rebuild_library_r2025b_20260920.m — rebuild SNS_Library with Rybak-style
% synapse icons (pass-through axon + postsynaptic terminal), re-validate, and
% render an icon preview. Canonical SNS_Library.slx becomes R2025b-native
% (regenerate on R2025a machines via sns_build_library + sns_build_actuators);
% an _R2025a copy is exported alongside for inspection.

sns = fileparts(fileparts(mfilename('fullpath')));   % dev -> SNS_Simscape
cd(sns); addpath(sns); addpath(fullfile(sns, 'demos'));
fprintf('=== MATLAB %s ===\n', version);

%% 0. lint the edited scripts
for f = {'sns_build_library.m', 'sns_export_diagram.m', fullfile('demos','sns_run_deng_demo.m')}
    msgs = checkcode(fullfile(sns, f{1}), '-severity');
    bad = msgs(~cellfun(@isempty, regexpi({msgs.message}, 'error|unterminated|might be missing')));
    fprintf('checkcode %-35s : %d messages (%d severe)\n', f{1}, numel(msgs), numel(bad));
    for k = 1:numel(bad), fprintf('  SEVERE L%d: %s\n', bad(k).line, bad(k).message); end
end

%% 1. rebuild library (7 blocks) + actuators (11 total)
run(fullfile(sns, 'sns_build_library.m'));
run(fullfile(sns, 'sns_build_actuators.m'));
n = numel(find_system('SNS_Library', 'SearchDepth', 1, 'Type', 'Block'));
fprintf('LIBRARY REBUILT: %d top-level blocks\n', n);

%% 2. re-validate numerics (BPA/BioMuscle vs Ben's equations)
run(fullfile(sns, 'sns_test_actuators.m'));

%% 3. icon preview: neuron -> E-syn -> neuron / neuron -> I-syn -> neuron chain
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
add_line(mdl, 'post_E/1', 'synE/2', 'autorouting', 'on');   % Vpost feedback wire
add_line(mdl, 'synE/1',   'post_E/1', 'autorouting', 'on');
add_line(mdl, 'pre_I/1',  'synI/1', 'autorouting', 'on');
add_line(mdl, 'post_I/1', 'synI/2', 'autorouting', 'on');
add_line(mdl, 'synI/1',   'post_I/1', 'autorouting', 'on');
set_param(mdl, 'FontSize', '12', 'ScreenColor', 'white');
print(['-s' mdl], '-dpng', '-r200', fullfile(sns, 'figures', 'sns_icon_preview.png'));
fprintf('icon preview written\n');
close_system(mdl, 0);

%% 4. R2025a copy of the rebuilt library (for R2025a machines; canonical
%%    regeneration route there = re-run the two build scripts)
Simulink.exportToVersion('SNS_Library', 'SNS_Library_R2025a', 'R2025a');
fprintf('EXPORTED SNS_Library_R2025a.slx\n');
close_system('SNS_Library', 0);
close_system('SNS_Library_R2025a', 0);
fprintf('=== REBUILD DONE ===\n');
