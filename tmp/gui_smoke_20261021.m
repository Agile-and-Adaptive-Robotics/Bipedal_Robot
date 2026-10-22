% GUI smoke test: build the demo circuit, capture screenshots, measure
% load/build timings, and verify the end-to-end build path (GUI JSON ->
% sns_build_from_json -> .slx). Runs headless in -batch (uifigure renders
% offscreen; exportapp captures it).
here = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape';
addpath(here); addpath(fullfile(here, 'gui')); addpath(fullfile(here, 'demos'));
compDir = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\60-sns\components';
spinal = 'C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal';

app = SNSGui();
app.demoKneeReflex();
drawnow; pause(0.5);
exportapp(app.f, fullfile(compDir, 'sns_gui_demo_circuit.png'));
fprintf('shot 1 ok (demo circuit)\n');

% ---- clicks-to-circuit (derived from the real interaction model) ----
% place N blocks: 1 palette click + 1 canvas click each
% C connections: 1 connect-mode toggle + 2 block clicks each
% 1 param edit (Esyn flip): 1 select + 1 field commit
clicksDemo = 2*numel(app.blocks) + (1 + 2*numel(app.edges)) + 2;
fprintf('CLICKS-TO-CIRCUIT (knee reflex, %d blocks %d connections): %d\n', ...
    numel(app.blocks), numel(app.edges), clicksDemo);

% ---- big-network load: the tuned production export ----
t0 = tic;
app.loadJson(fullfile(spinal, 'spinal_net_export.json'));
tLoad410 = toc(t0);
drawnow; pause(1.0);
exportapp(app.f, fullfile(compDir, 'sns_gui_production_network.png'));
fprintf('shot 2 ok (production network, %.2f s load)\n', tLoad410);

% ---- winner slice load ----
t0 = tic;
app.loadJson(fullfile(spinal, 'reports_20261016', 'sns_simulink', 'slice_export.json'));
tSlice = toc(t0);
fprintf('winner-slice load: %.2f s\n', tSlice);

% ---- end-to-end build from the GUI model (demo circuit again) ----
app.demoKneeReflex();
app.cbBuild();
drawnow; pause(0.3);
slxPath = fullfile(here, 'gui', 'SNS_GuiCircuit.slx');
fprintf('build output exists: %d\n', isfile(slxPath));
t0 = tic; load_system(slxPath); tOpen = toc(t0);
nB = numel(find_system('SNS_GuiCircuit', 'SearchDepth', 1, 'Type', 'Block'));
fprintf('GUI-built model opens in %.2f s, %d top-level blocks\n', tOpen, nB);
close_system('SNS_GuiCircuit', 0);

fprintf('GUI SMOKE DONE. nClicks instrumented=%d, load410=%.2fs, build=%.2fs\n', ...
    app.nClicks, app.timings.loadJson, app.timings.build);
delete(app.f);
