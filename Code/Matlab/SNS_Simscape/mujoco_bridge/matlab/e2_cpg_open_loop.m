function e2_cpg_open_loop()
% E2: the FULL tuned spinal network (SNS_SpinalNetwork, 410 neurons) driving
% ALL 92 muscles of the MuJoCo gait2392 model through the bridge,
% open-loop (no afferents), DRIVE = 2.5 nA constant.
% Solver: fixed-step ode1 @ 2 ms = the production Euler semantics the v4b
% tuning was optimized under (network steps at its native dt; Rate
% Transition ZOHs S into the 5 ms plant).
% Outputs: exp_out\e2_cpg.mat + e2_cpg.png + harness model .slx.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));               % SNS_Simscape (SNS_Library)
addpath(fullfile(here, '..', '..', 'results'));    % SNS_SpinalNetwork.slx
xml = strrep(fullfile(here, '..', '..', '..', '..', '..', 'Solid_Models', ...
    'OpenSim', 'Gait2392_Robotbody', 'mjc', 'gait2392_simbody', ...
    'gait2392_simbody_cvt3_simbridge2.xml'), '\', '/');  % 2 ms, as the runner

mdl = 'sns_bridge_cpg';
if bdIsLoaded(mdl), close_system(mdl, 0); end
outDir = fullfile(here, 'exp_out');
if ~isfolder(outDir), mkdir(outDir); end
if exist(fullfile(outDir, [mdl '.slx']), 'file')
    delete(fullfile(outDir, [mdl '.slx']));
end
new_system(mdl);

%% ---------------- network (model reference) + drive -----------------------
load_system('SNS_SpinalNetwork');
% make the saved default the production semantics (see README: chaos note)
set_param('SNS_SpinalNetwork', 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002');
save_system('SNS_SpinalNetwork');
add_block('simulink/Ports & Subsystems/Model', [mdl '/cpg'], ...
          'ModelFile', 'SNS_SpinalNetwork.slx', ...
          'Position', [200 100 320 220]);
u = zeros(376, 1); u(1) = 2.5;                     % DRIVE port is u(1)
add_block('simulink/Sources/Constant', [mdl '/uDrive'], 'Value', 'uvec', ...
          'Position', [60 130 130 170]);
add_line(mdl, 'uDrive/1', 'cpg/1', 'autorouting', 'on');

%% ---------------- rate transition + plant --------------------------------
add_block('simulink/Signal Attributes/Rate Transition', [mdl '/RT'], ...
          'Position', [420 130 460 170]);
add_line(mdl, 'cpg/1', 'RT/1', 'autorouting', 'on');
add_block('mjLib/MuJoCo Plant', [mdl '/Plant'], ...
    'xmlFileRel', xml, 'renderingType', 'None', ...
    'rgbOutOption', 'off', 'depthOutOption', 'off', ...
    'Position', [540 100 700 240]);
add_line(mdl, 'RT/1', 'Plant/1', 'autorouting', 'on');

%% ---------------- logging --------------------------------------------------
% S vector (92 MN drives) + sensor bus
ph = get_param([mdl '/cpg'], 'PortHandles');
set_param(ph.Outport(1), 'DataLogging', 'on', ...
          'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'S');
add_block('simulink/Signal Routing/Bus Selector', [mdl '/sens'], ...
    'OutputSignals', 'knee_r_pos,vas_med_r_len,vas_med_r_frc', ...
    'Position', [760 100 780 170]);
add_line(mdl, 'Plant/1', 'sens/1', 'autorouting', 'on');
for k = 1:3
    ph2 = get_param([mdl '/sens'], 'PortHandles');
    set_param(ph2.Outport(k), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', char("sens" + k));
end

%% ---------------- run ------------------------------------------------------
assignin('base', 'uvec', u);
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002', 'StopTime', '5', ...
          'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');

S  = out.sigs.get('S').Values;      %#ok<NASGU>
Sd = squeeze(S.Data);
if iscolumn(Sd)                      % logger may flatten the 92-wide signal
    Sd = reshape(Sd, [], 92);
end
knee = out.sigs.get('sens1').Values;
len  = out.sigs.get('sens2').Values;
frc  = out.sigs.get('sens3').Values;
save_system(mdl, fullfile(outDir, [mdl '.slx']));
close_system(mdl, 0);

%% ---------------- report + figure -----------------------------------------
fprintf('E2 (full CPG -> 92 muscles, DRIVE=2.5, 5 s):\n');
fprintf(['  S drive: mean %.3f, max %.3f; fraction of muscles with ' ...
         'S>0.1 at some point: %.0f%%\n'], mean(mean(Sd)), max(Sd(:)), ...
        100 * mean(any(Sd > 0.1, 1)));
fprintf('  knee_r: %.1f -> %.1f deg (min %.1f, max %.1f)\n', ...
        rad2deg(knee.Data(1)), rad2deg(knee.Data(end)), ...
        min(rad2deg(knee.Data(:))), max(rad2deg(knee.Data(:))));
fprintf('  vas_med_r |frc| max %.0f N; all S finite: %d\n', ...
        max(abs(frc.Data(:))), all(isfinite(Sd(:))));

% vas_med_r drive index (actuator 28 -> S column 29)
fig = figure('Visible', 'off', 'Position', [100 100 900 620]);
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact');
nexttile(tl); plot(S.Time, Sd(:, 29), 'LineWidth', 1.1);
grid on; ylabel('S vas\_med\_r'); title('E2: tuned CPG drives MuJoCo gait2392 (open loop)');
nexttile(tl); plot(knee.Time, rad2deg(knee.Data(:)), 'LineWidth', 1.1);
grid on; ylabel('knee\_r (deg)');
nexttile(tl); plot(frc.Time, frc.Data(:), 'LineWidth', 1.1);
grid on; ylabel('vas\_med\_r force (N)'); xlabel('t (s)');
exportgraphics(fig, fullfile(outDir, 'e2_cpg.png'), 'Resolution', 130);
save(fullfile(outDir, 'e2_cpg.mat'), 'S', 'knee', 'len', 'frc');
fprintf('saved %s\n', fullfile(outDir, 'e2_cpg.mat'));
end
