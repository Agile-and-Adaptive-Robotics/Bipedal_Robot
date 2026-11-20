function sns_build_from_json(jsonPath, mdlName, slxDst)
% Generate the editable Simulink spinal-network model from SNS_Library
% blocks, driven by Code\MuJoCo_SNS\spinal\spinal_net_export.json
% (the tuned runner --fitted --best configuration, exported 2026-09-12).
%
% sns_build_from_json()                    -- the default tuned export
% sns_build_from_json(json, mdl, slxDst)   -- GUI/graph caller form:
%                                           any network JSON of the same
%                                           schema; inputs/outputs may be
%                                           absent (the inport/S-mux
%                                           sections are skipped).
%
% Units mapping (verified by sns_units_test_2n, max dev 6e-4 mV):
%   toolbox NonSpikingNeuron C=tau uF, Gm=1 uS, Vrest=0  ->  Cm=1000*tau nF
%   synapse e_lo/e_hi = 0/5 mV  ->  ThrPre=0, SlopePre=5
%
% Model layout: one flat diagram, columns by cell class, one wide input
% inport u [376 x 1] (Demux per input port, ordering = JSON inputs[]).
% 2026-09-22 architecture: neurons sum their synaptic inputs INTERNALLY —
% each NonSpikingSynapse (one input: Vpre) lands on a syn1..syn6 port of
% its postsynaptic neuron; external input currents are Muxed onto the
% neuron's Iapp port (element-wise summed inside); neurons with MORE than
% 6 incoming synapses (up to 17 here) route the excess through chained
% SynSum junction blocks into the last syn port. MN S(V) outputs are
% Muxed (92 wide) by actuator id (= MuJoCo ctrl order) -> outport S.

here = fileparts(mfilename('fullpath'));
code = fullfile(here, '..', '..');            % ...\Code
if nargin < 1 || isempty(jsonPath)
    jsonPath = fullfile(code, 'MuJoCo_SNS', 'spinal', ...
                        'spinal_net_export.json');
end
if nargin < 2 || isempty(mdlName), mdlName = 'SNS_SpinalNetwork'; end
if nargin < 3 || isempty(slxDst)
    slxDst = fullfile(here, 'results', [mdlName '.slx']);
end
jtxt = fileread(jsonPath);
J = jsondecode(jtxt);

neurons = J.neurons;
synapses = J.synapses;
if isfield(J, 'inputs'), inputs = J.inputs; else, inputs = struct([]); end
if isfield(J, 'outputs'), outputs = J.outputs; else, outputs = struct([]); end
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

mdl = mdlName;
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(slxDst, 'file'), delete(slxDst); end
new_system(mdl);
set_param(mdl, 'UnconnectedInputMsg', 'none');   % unused syn ports ground

% ---- inport + one N-way demux (port k of u) -------------------------------
if ~isempty(inputs)
add_block('simulink/Sources/In1', [mdl '/u'], 'PortDimensions', ...
          num2str(numel(inputs)), 'Position', [40 40 70 120]);
add_block('simulink/Signal Routing/Demux', [mdl '/demux_u'], ...
          'Outputs', num2str(numel(inputs)), 'Position', ...
          [150 40 155 40 + 12 * numel(inputs)]);
add_line(mdl, 'u/1', 'demux_u/1', 'autorouting', 'on');
end
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

% ---- synapses (band left of the destination column, y at destination) ----
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
        'Position', [colX(cat{d}) - 170 y colX(cat{d}) - 126 y + 32]);
    set_param([mdl '/' nm], 'ShowName', 'off');
end

% ---- wiring: inputs -> Iapp (mux), synapses -> syn ports (+SynSum chain) --
incomingIn = cell(nN, 1);   % external input lines per neuron
incomingSyn = cell(nN, 1);  % synapse block names per neuron
for k = 1:nN
    incomingIn{k} = {};
    incomingSyn{k} = {};
end
for k = 1:numel(inputs)
    d = ipos(inputs(k).dst);
    incomingIn{d}{end + 1} = sprintf('demux_u/%d', k); %#ok<SAGROW>
end
for k = 1:numel(synapses)
    s = synapses(k);
    d = ipos(s.dst);
    incomingSyn{d}{end + 1} = sprintf('syn_%d_%s_to_%s', k, s.src, s.dst); %#ok<SAGROW>
end

NSYNPORT = 6;   % neuron syn1..syn6 ports
nSS = 0;
for k = 1:nN
    nm = neurons(k).name;
    x = colX(cat{k});

    % external currents -> Iapp port 1 (vector inputs sum element-wise)
    if ~isempty(incomingIn{k})
        if numel(incomingIn{k}) == 1
            add_line(mdl, incomingIn{k}{1}, [nm '/1'], 'autorouting', 'on');
        else
            muxn = ['muxIn_' nm];
            add_block('simulink/Signal Routing/Mux', [mdl '/' muxn], ...
                'Inputs', num2str(numel(incomingIn{k})), ...
                'Position', [x - 120 ny(k) - 6 x - 117 ny(k) + 30]);
            set_param([mdl '/' muxn], 'ShowName', 'off');
            for j = 1:numel(incomingIn{k})
                add_line(mdl, incomingIn{k}{j}, sprintf('%s/%d', muxn, j), ...
                    'autorouting', 'on');
            end
            add_line(mdl, [muxn '/1'], [nm '/1'], 'autorouting', 'on');
        end
    end

    % synapses -> syn ports; >6 -> chained SynSum into the last port
    ns = numel(incomingSyn{k});
    if ns == 0
        continue
    end
    ndirect = min(ns, 5);                 % keep >=1 port for a SynSum if needed
    if ns <= NSYNPORT
        ndirect = ns;
    end
    for j = 1:ndirect
        add_line(mdl, [incomingSyn{k}{j} '/1'], sprintf('%s/%d', nm, j + 1), ...
            'autorouting', 'on');
    end
    rest = incomingSyn{k}(ndirect + 1:end);
    if ~isempty(rest)
        % chain SynSum junctions: first takes up to 8 synapse lines, each
        % next takes the previous junction's output + up to 7 more
        prevSrc = '';
        i = 1;
        while i <= numel(rest)
            nSS = nSS + 1;
            ss = sprintf('ss_%d_%s', nSS, nm);
            add_block('SNS_Library/SynSum', [mdl '/' ss], ...
                'Position', [x - 260 ny(k) + 20 + 26*nSS x - 216 ny(k) + 52 + 26*nSS]);
            set_param([mdl '/' ss], 'ShowName', 'off');
            pn = 0;
            if ~isempty(prevSrc)
                pn = pn + 1;
                add_line(mdl, [prevSrc '/1'], sprintf('%s/%d', ss, pn), ...
                    'autorouting', 'on');
            end
            while pn < 8 && i <= numel(rest)
                pn = pn + 1;
                add_line(mdl, [rest{i} '/1'], sprintf('%s/%d', ss, pn), ...
                    'autorouting', 'on');
                i = i + 1;
            end
            prevSrc = ss;
        end
        add_line(mdl, [prevSrc '/1'], sprintf('%s/%d', nm, NSYNPORT + 1), ...
            'autorouting', 'on');
    end
end

% ---- synapse Vpre wiring (one input per synapse) ---------------------------
for k = 1:numel(synapses)
    s = synapses(k);
    synbl = sprintf('syn_%d_%s_to_%s', k, s.src, s.dst);
    add_line(mdl, [s.src '/1'], [synbl '/1'], 'autorouting', 'on');  % V -> Vpre
end

% ---- MN S outputs -> Mux -> outport (actuator-id order) ------------------
if ~isempty(outputs)
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
end

% ---- model annotation + config -------------------------------------------
try
    if isfield(J, 'meta') && isfield(J.meta, 'source'), src = J.meta.source;
    else, src = jsonPath; end
    a = Simulink.Annotation(mdl, sprintf( ...
        ['SNS network (SNS_Library blocks, 2026-09-22 architecture: ' ...
         'synapses one-input onto neuron syn ports, summation inside ' ...
         'neurons).\nSource: %s'], src));
    a.Position = [40 200; 40 200];
catch
    % annotation API varies; model is complete without it
end
set_param(mdl, 'Description', sprintf( ...
    ['SNS network generated from %s. u = input currents nA; S = MN ' ...
     'drives in actuator order.'], jsonPath));
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002');   % production Euler semantics (chaos note)

save_system(mdl, slxDst);
fprintf('saved %s (%d SynSum junctions for >6-synapse neurons)\n', slxDst, nSS);

% ---- smoke run: compile + 0.1 s with zero input ---------------------------
set_param(mdl, 'StopTime', '0.1', 'SignalLogging', 'off');
sim(mdl);
fprintf('smoke run 0.1 s OK (zero input)\n');
close_system(mdl, 0);
end
