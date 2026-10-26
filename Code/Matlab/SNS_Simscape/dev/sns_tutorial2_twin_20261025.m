function sns_tutorial2_twin()
% Simulink twin of sns-toolbox Tutorial 2: 5-neuron network with
% excitatory (Erev 50), disynaptic inhibitory (Dest2In -40 onto Dest2),
% and modulatory (Dest1 -> Dest3, Erev 0) synapses, all gmax 1.0 uS,
% toolbox-default neurons (tau 50 ms), 20 nA constant input, 50 ms.
% Mapping (units-test convention): Vrest 0, Gm 0.1 uS, Cm 5 nF,
% ThrPre 0 / SlopePre 20 (toolbox synapse defaults e_lo/e_hi 0/20).
% Reference: results/tutorial2_ref.mat (numpy backend, dt 0.01 ms).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
cd(root); addpath(root);
load_system('SNS_Library');

mdl = 'sns_tutorial2_twin';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(fullfile(root, 'results', [mdl '.slx']), 'file')
    delete(fullfile(root, 'results', [mdl '.slx']));
end
new_system(mdl);

N = @(nm, y) add_block('SNS_Library/NonSpikingNeuron', [mdl '/' nm], ...
    'Vrest', '0', 'Gm', '1', 'Cm', '5', 'Thr', '0', 'Slope', '20', ...
    'Position', [220 y 320 y + 90]);
N('SourceNrn', 40); N('Dest1', 160); N('Dest2', 280); N('Dest2In', 400); N('Dest3', 520);
S = @(nm, esyn, y) add_block('SNS_Library/NonSpikingSynapse', [mdl '/' nm], ...
    'gmax', '1', 'Esyn', num2str(esyn), 'ThrPre', '0', 'SlopePre', '20', ...
    'ShowName', 'off', 'Position', [150 y 195 y + 34]);
S('s_src_d1', 50, 100); S('s_src_d2', 50, 220); S('s_d2_in', 50, 340);
S('s_inh_d2', -40, 430); S('s_src_d3', 50, 480); S('s_mod_d3', 0, 580);
add_block('simulink/Sources/Constant', [mdl '/I20'], 'Value', '20', ...
    'Position', [60 40 120 80]);
add_line(mdl, 'I20/1', 'SourceNrn/1', 'autorouting', 'on');
add_line(mdl, 'SourceNrn/1', 's_src_d1/1', 'autorouting', 'on');
add_line(mdl, 's_src_d1/1', 'Dest1/2', 'autorouting', 'on');
add_line(mdl, 'SourceNrn/1', 's_src_d2/1', 'autorouting', 'on');
add_line(mdl, 's_src_d2/1', 'Dest2/2', 'autorouting', 'on');
add_line(mdl, 'Dest2/1', 's_d2_in/1', 'autorouting', 'on');
add_line(mdl, 's_d2_in/1', 'Dest2In/2', 'autorouting', 'on');
add_line(mdl, 'Dest2In/1', 's_inh_d2/1', 'autorouting', 'on');
add_line(mdl, 's_inh_d2/1', 'Dest2/3', 'autorouting', 'on');
add_line(mdl, 'SourceNrn/1', 's_src_d3/1', 'autorouting', 'on');
add_line(mdl, 's_src_d3/1', 'Dest3/2', 'autorouting', 'on');
add_line(mdl, 'Dest1/1', 's_mod_d3/1', 'autorouting', 'on');
add_line(mdl, 's_mod_d3/1', 'Dest3/3', 'autorouting', 'on');

logn = {'SourceNrn', 'Dest1', 'Dest2', 'Dest2In', 'Dest3'};
for k = 1:numel(logn)
    ph = get_param([mdl '/' logn{k}], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', 'DataLoggingName', logn{k});
end

set_param(mdl, 'Solver', 'ode45', 'RelTol', '1e-8', 'AbsTol', '1e-10', ...
          'StopTime', '0.05', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs', 'UnconnectedInputMsg', 'none');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
save_system(mdl, fullfile(root, 'results', [mdl '.slx']));
close_system(mdl, 0);

ref = load(fullfile(root, 'results', 'tutorial2_ref.mat'));
refmap = {'SourceNrn', 'Src'; 'Dest1', 'Dest1'; 'Dest2', 'Dest2'; ...
          'Dest2In', 'Dest2In'; 'Dest3', 'Dest3'};
dmax = 0;
fprintf('TUTORIAL-2 TWIN final values (mV):\n');
for k = 1:size(refmap, 1)
    nm = refmap{k, 1};
    v = squeeze(out.sigs.get(nm).Values.Data);
    r = double(ref.(refmap{k, 2})(:));
    d = abs(v(end) - r(end));
    dmax = max(dmax, d);
    fprintf('  %-10s twin %8.4f | ref %8.4f | |d| %.2e\n', ...
            nm, v(end), r(end), d);
end
if dmax < 0.05
    fprintf('TUTORIAL-2 TWIN: PASS (all finals within 0.05 mV of toolbox)\n');
else
    fprintf('TUTORIAL-2 TWIN: CHECK FAILED (max |d| %.3e mV)\n', dmax);
end
end
