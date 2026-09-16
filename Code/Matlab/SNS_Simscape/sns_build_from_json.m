function sns_build_from_json()
% Generate the editable Simulink spinal-network model from SNS_Library
% blocks, driven by Code\MuJoCo_SNS\spinal\spinal_net_export.json
% (the tuned runner --fitted --best configuration, exported 2026-09-12).
%
% Units mapping (verified by sns_units_test_2n, max dev 6e-4 mV):
%   toolbox NonSpikingNeuron C=tau uF, Gm=1 uS, Vrest=0  ->  Cm=1000*tau nF
%   synapse e_lo/e_hi = 0/5 mV  ->  ThrPre=0, SlopePre=5
%
% Model layout: one flat diagram, columns by cell class, one wide input
% inport u [376 x 1] (Selector per input port, ordering = JSON inputs[]),
% per-neuron Sum blocks combining incoming synapse currents + external
% input current, and a 92-wide Mux of the MN S(V) outputs ordered by
% actuator id (= MuJoCo ctrl order) -> outport S.

here = fileparts(mfilename('fullpath'));
code = fullfile(here, '..', '..');            % ...\Code
jtxt = fileread(fullfile(code, 'MuJoCo_SNS', 'spinal', ...
                         'spinal_net_export.json'));
J = jsondecode(jtxt);

neurons = J.neurons;
synapses = J.synapses;
inputs = J.inputs;
outputs = J.outputs;
nN = numel(neurons);
fprintf('JSON: %d neurons, %d synapses, %d inputs, %d MN outputs\n', ...
        nN, numel(synapses), numel(inputs), numel(outputs));

% ---- index helpers --------------------------------------------------------
ipos = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:nN
    ipos(neurons(k).name) = k;
end
% category column for layout
cat = cell(nN, 1);
for k = 1:nN
    nm = neurons(k).name;
    if contains(nm, "MN_"), c = 6;
    elseif nm(1) == "I" && (startsWith(nm, "Ia_") || startsWith(nm, "II_") || startsWith(nm, "Ib_")), c = 5;
    elseif startsWith(nm, "IBEXC"), c = 4;
    elseif startsWith(nm, "PF"), c = 3;
    elseif startsWith(nm, "RG") || startsWith(nm, "ADAP") || startsWith(nm, "PRESET"), c = 2;
    else, c = 1;   % DRIVE / POSTURE / BAL_* shared cells
    end
    cat{k} = c;
end
ncol = 6;
colX = [560 720 880 1040 1200 1360];
yCtr = zeros(ncol, 1);
Y0 = 60; dY = 96;

mdl = 'SNS_SpinalNetwork';
if bdIsLoaded(mdl), close_system(mdl, 0); end
slxDst = fullfile(here, 'results', [mdl '.slx']);
if exist(slxDst, 'file'), delete(slxDst); end
new_system(mdl);

% ---- inport + one 376-way demux (port k of u) -----------------------------
add_block('simulink/Sources/In1', [mdl '/u'], 'PortDimensions', ...
          num2str(numel(inputs)), 'Position', [40 40 70 120]);
add_block('simulink/Signal Routing/Demux', [mdl '/demux_u'], ...
          'Outputs', num2str(numel(inputs)), 'Position', ...
          [150 40 155 40 + 12 * numel(inputs)]);
add_line(mdl, 'u/1', 'demux_u/1', 'autorouting', 'on');
for k = 1:numel(inputs)
    d = ipos(inputs(k).dst);
    y = Y0 + (yCtr(cat{d})) * dY; yCtr(cat{d}) = yCtr(cat{d}) + 1;
end

% ---- neurons --------------------------------------------------------------
ny = zeros(nN, 1);
for k = 1:nN
    c = cat{k};
    y = Y0 + yCtr(c) * dY; yCtr(c) = yCtr(c) + 1; ny(k) = y;
    add_block('SNS_Library/NonSpikingNeuron', [mdl '/' neurons(k).name], ...
        'Vrest', num2str(neurons(k).Vrest_mV), ...
        'Gm', num2str(neurons(k).Gm_uS), ...
        'Cm', num2str(neurons(k).Cm_nF), ...
        'Thr', '0', 'Slope', '5', ...
        'Position', [colX(c) y colX(c) + 80 y + 80]);
end

% ---- synapses (band left of the neuron columns, y at destination) --------
for k = 1:numel(synapses)
    s = synapses(k);
    d = ipos(s.dst);
    y = ny(d) + 8;
    nm = sprintf('syn_%d_%s_to_%s', k, s.src, s.dst);
    add_block('SNS_Library/NonSpikingSynapse', [mdl '/' nm], ...
        'gmax', num2str(s.g_uS, 12), ...
        'Esyn', num2str(s.Esyn_mV, 12), ...
        'ThrPre', num2str(s.ThrPre_mV, 12), ...
        'SlopePre', num2str(s.SlopePre_mV, 12), ...
        'Position', [340 y 440 y + 80]);
end

% ---- per-neuron sums + wiring --------------------------------------------
incoming = cell(nN, 1);
for k = 1:nN, incoming{k} = {}; end
for k = 1:numel(inputs)
    d = ipos(inputs(k).dst);
    incoming{d}{end + 1} = sprintf('demux_u/%d', k); %#ok<SAGROW>
end
for k = 1:numel(synapses)
    s = synapses(k);
    d = ipos(s.dst);
    incoming{d}{end + 1} = sprintf('syn_%d_%s_to_%s/1', k, s.src, s.dst); %#ok<SAGROW>
end
for k = 1:nN
    nm = neurons(k).name;
    inc = incoming{k};
    if isempty(inc)
        continue                       % dangling Isyn (none expected)
    elseif numel(inc) == 1
        add_line(mdl, inc{1}, [nm '/1'], 'autorouting', 'on');
    else
        add_block('simulink/Math Operations/Sum', [mdl '/sum_' nm], ...
                  'Inputs', repmat('+', 1, numel(inc)), ...
                  'Position', [colX(cat{k}) - 90 ny(k) + 24 ...
                               colX(cat{k}) - 60 ny(k) + 56]);
        for j = 1:numel(inc)
            add_line(mdl, inc{j}, sprintf('sum_%s/%d', nm, j), ...
                     'autorouting', 'on');
        end
        add_line(mdl, ['sum_' nm '/1'], [nm '/1'], 'autorouting', 'on');
    end
end

% ---- synapse Vpre/Vpost wiring -------------------------------------------
for k = 1:numel(synapses)
    s = synapses(k);
    synbl = sprintf('syn_%d_%s_to_%s', k, s.src, s.dst);
    add_line(mdl, [s.src '/1'], [synbl '/1'], 'autorouting', 'on');  % V -> Vpre
    add_line(mdl, [s.dst '/1'], [synbl '/2'], 'autorouting', 'on');  % V -> Vpost
end

% ---- MN S outputs -> Mux -> outport (actuator-id order) ------------------
[~, ord] = sort([outputs.actuator_id]);
add_block('simulink/Signal Routing/Mux', [mdl '/S_mux'], ...
          'Inputs', num2str(numel(outputs)), 'Position', ...
          [colX(ncol) + 160 60 colX(ncol) + 163 60 + 20 * numel(outputs)]);
for j = 1:numel(ord)
    o = outputs(ord(j));
    add_line(mdl, [o.mn '/2'], sprintf('S_mux/%d', j), 'autorouting', 'on');
end
add_block('simulink/Sinks/Out1', [mdl '/S'], 'Position', ...
          [colX(ncol) + 260 60 colX(ncol) + 290 140]);
add_line(mdl, 'S_mux/1', 'S/1', 'autorouting', 'on');

% ---- model annotation + config -------------------------------------------
try
    a = Simulink.Annotation(mdl, sprintf( ...
        ['Tuned gait2392 spinal network (SNS_Library blocks).\n' ...
         'Source: %s\nu [376 x 1] input currents (nA), ordering = ' ...
         'spinal_net_export.json inputs[] (1:8 DRIVE..BAL_LAT_L, then ' ...
         'per-muscle POST_/Ia_/II_/Ib_).\nS [92 x 1] MN drives ordered by ' ...
         'MuJoCo actuator id (ctrl order).\nUnits mapping verified by ' ...
         'sns_units_test_2n (2026-09-12).'], J.meta.source));
    a.Position = [40 200; 40 200];
catch
    % annotation API varies; model is complete without it
end
set_param(mdl, 'Description', sprintf( ...
    ['Tuned gait2392 spinal network generated from ' ...
     'spinal_net_export.json (%s). u = input currents nA; S = MN drives ' ...
     'in actuator order.'], J.meta.source));
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002');   % production Euler semantics (chaos note)

save_system(mdl, slxDst);
fprintf('saved %s\n', slxDst);

% ---- smoke run: compile + 0.1 s with zero input ---------------------------
set_param(mdl, 'StopTime', '0.1', 'SignalLogging', 'off');
sim(mdl);
fprintf('smoke run 0.1 s OK (zero input)\n');
close_system(mdl, 0);
end
