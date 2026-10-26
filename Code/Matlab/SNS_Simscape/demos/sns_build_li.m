%% sns_build_li.m - Li two-layer contact-driven CPG as a base-Simulink model
%
% Ben's order (2026-10-24): rebuild the Li model from the AnimatLab
% original into SNS Simscape. Source of truth: the GUID-level mine in
% Code\MuJoCo_SNS\spinal\connectome_templates.json key "li" (58 nodes /
% 98 edges; "walk tester rearranged.aproj", 2023), which carries a
% PER-SYNAPSE conductance (1-8 uS) on every edge - preserved here
% individually per Ben's 2026-10-24 per-target ruling.
%
% Fidelity notes (LOUD, per the no-silent-simplification rule):
%   - Wiring, signs, and per-synapse conductances: FAITHFUL to the mine.
%   - Stance-CPG dynamics: the mine documents "per-leg stance CPG
%     (6 nA tonic, tau_Ca 7.5 s burst termination) gated by CPG-inhibit"
%     and flags the original as spiking I&F. This v1 implements the
%     documented slow burst termination with the codebase's PROVEN graded
%     pattern (NonSpikingNeuron + slow Adp neuron, tau scaled to 7.5 s),
%     NOT a spike-event Ca integrator - an ADAPTATION, labeled as such;
%     swapping in the spiking cell is local to build_stance_cpg().
%   - SN-heel/SN-toe sensors are driven by a scripted contact script
%     (antiphase heel pulses) for this neural smoke run.
%
% Output: demos\SNS_Li_CPG.slx. Smoke runner: dev\run_li_smoke_20261024.m

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);

tpl = jsondecode(fileread(fullfile( ...
    root, '..', '..', 'MuJoCo_SNS', 'spinal', 'connectome_templates.json')));
li = tpl.li;
mdl = 'SNS_Li_CPG';
if bdIsLoaded(mdl), close_system(mdl, 0); end
slxPath = fullfile(here, [mdl '.slx']);
if exist(slxPath, 'file'), delete(slxPath); end
new_system(mdl);
set_param(mdl, 'UnconnectedInputMsg', 'none', 'ScreenColor', 'white');

%% node registry: canonical block name + type + position for every node
NODE = struct();   % scalar empty struct: allows NODE.(lbl) dynamic fields
labels = cell(numel(li.nodes), 1);
for k = 1:numel(li.nodes)
    nd = li.nodes(k);
    lbl = matlab.lang.makeValidName(nd.label);
    labels{k} = lbl;
    p = [60 60 140 140];
    if isfield(nd, 'x') && isfield(nd, 'y')
        p = [nd.x nd.y nd.x + 110 nd.y + 80];
    end
    switch nd.type
        case 'PORT-load',   NODE.(lbl).blk = ['in_' lbl];
        case 'HC-RG-E',     NODE.(lbl).blk = ['cpg_' lbl];
        case 'HC-RG-F',     NODE.(lbl).blk = ['cpg_' lbl];
        case 'MUSCLE',      NODE.(lbl).blk = ['m_' lbl];
        case {'SN-heel', 'SN-toe', 'SN-II'}, NODE.(lbl).blk = ['sn_' lbl];
        otherwise,         NODE.(lbl).blk = lbl;   % MN, IN-*, HC-PF-*
    end
    NODE.(lbl).type = nd.type;
    NODE.(lbl).pos = p;
end
fn = fieldnames(NODE);

%% create blocks
for k = 1:numel(fn)
    lbl = fn{k};
    b = NODE.(lbl).blk;
    p = NODE.(lbl).pos;
    switch NODE.(lbl).type
        case 'PORT-load'
            add_block('simulink/Sources/Constant', [mdl '/' b], ...
                'Value', '6', 'Position', p);       % 6 nA tonic (mine note)
        case {'HC-RG-E', 'HC-RG-F'}
            build_stance_cpg(mdl, b, p);
        case 'MN'
            add_block('SNS_Library/NonSpikingNeuron', [mdl '/' b], ...
                'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', ...
                'Slope', '1', 'Position', p);
        case 'MUSCLE'
            add_block('simulink/Sinks/Terminator', [mdl '/' b], ...
                'Position', p);
        case {'SN-heel', 'SN-toe'}
            % scripted contact sensor: 10 nA heel pulses (antiphase by side)
            ph = '0'; if contains(lbl, 'R'), ph = '1'; end %#ok<AGROW>
            add_block('SNS_Library/NonSpikingNeuron', [mdl '/' b], ...
                'Vrest', '-52', 'Gm', '0.12', 'Cm', '2', 'Thr', '-55', ...
                'Slope', '1', 'Position', p);
            add_block('simulink/Sources/Pulse Generator', ...
                [mdl '/contact_' b], 'PulseType', 'Time based', ...
                'Amplitude', '10', 'Period', '2', 'PulseWidth', '35', ...
                'PhaseDelay', ph, 'Position', [p(1) - 130 p(2) p(1) - 60 p(2) + 30]);
            add_line(mdl, ['contact_' b '/1'], [b '/1'], 'autorouting', 'on');
        otherwise
            add_block('SNS_Library/NonSpikingNeuron', [mdl '/' b], ...
                'Vrest', '-52', 'Gm', '0.12', 'Cm', '5', 'Thr', '-50', ...
                'Slope', '1', 'Position', p);
    end
end

%% per-target synapses: one block per edge, each with ITS OWN gain.
% Port allocator: neurons take synaptic inputs on ports 2..7 (syn1..syn6);
% stance CPGs on ports 2..4; MN->MUSCLE edges are activation taps (no syn).
NEXTPORT = struct();
SYN_X = 40;
EDGE_N = 0; SKIPPED = 0;
for k = 1:numel(li.edges)
    e = li.edges(k);
    src = matlab.lang.makeValidName(e.from);
    dst = matlab.lang.makeValidName(e.to);
    if ~isfield(NODE, src) || ~isfield(NODE, dst)
        SKIPPED = SKIPPED + 1;
        warning('edge %d skipped: endpoint missing (%s -> %s)', k, e.from, e.to);
        continue
    end
    sb_ = [NODE.(src).blk '/1'];                       % source V output
    if strcmp(NODE.(dst).type, 'MUSCLE')
        add_line(mdl, [NODE.(src).blk '/2'], [NODE.(dst).blk '/1'], ...
                 'autorouting', 'on');                % MN S(activation) tap
        EDGE_N = EDGE_N + 1;
        continue
    end
    % next free syn input port on the target
    if ~isfield(NEXTPORT, dst), NEXTPORT.(dst) = 2; end
    pnum = NEXTPORT.(dst);
    if strcmp(NODE.(dst).type, 'HC-RG-E') || strcmp(NODE.(dst).type, 'HC-RG-F')
        if pnum > 4, error('stance CPG %s out of syn ports at edge %d', dst, k); end
    elseif pnum > 7
        error('neuron %s exceeds 6 syn inputs (edge %d) - needs SynSum chain', dst, k);
    end
    NEXTPORT.(dst) = pnum + 1;
    EDGE_N = EDGE_N + 1;
    sb = sprintf('syn%03d', EDGE_N);
    esyn = '0';
    if strcmp(e.sign, 'inh'), esyn = '-72'; end
    add_block('SNS_Library/NonSpikingSynapse', [mdl '/' sb], ...
        'gmax', num2str(e.gain), 'Esyn', esyn, ...
        'ThrPre', '-50', 'SlopePre', '1', 'ShowName', 'off', ...
        'Position', [SYN_X + 160 * fix(EDGE_N / 34) 40 + 26 * EDGE_N ...
                     SYN_X + 205 + 160 * fix(EDGE_N / 34) 66 + 26 * EDGE_N]);
    add_line(mdl, sb_, [sb '/1'], 'autorouting', 'on');
    add_line(mdl, [sb '/1'], [NODE.(dst).blk '/' num2str(pnum)], ...
             'autorouting', 'on');
end
fprintf('SNS_Li_CPG: %d nodes, %d connections (%d skipped)\n', ...
    numel(fn), EDGE_N, SKIPPED);

set_param(mdl, 'SolverType', 'Variable-step', 'Solver', 'ode45', ...
          'MaxStep', '1e-3', 'StopTime', '30', 'UnconnectedInputMsg', 'none');
save_system(mdl, slxPath);
fprintf('built demos/SNS_Li_CPG.slx (wiring + per-synapse gains faithful to the mine)\n');

%% ---------------- local functions ----------------
function build_stance_cpg(mdl, b, p)
% Graded stance CPG v1 (ADAPTATION, flagged in the header): SNS
% NonSpikingNeuron + slow Adp neuron with tau scaled to the documented
% tau_Ca 7.5 s (Cm 750 nF / Gm 0.1 uS = 7500 ms). Adp S output feeds a
% LINEAR adaptive-inhibition synapse onto the cell (BPACPG-demo pattern;
% linear gain, not a saturating synapse, or the driving force collapses).
% Ports: 1 = Iapp (tonic drive / commissural currents), 2..4 = syn1..3
% external synaptic inputs, 1 out = cell V.
cp = [mdl '/' b];
add_block('simulink/Ports & Subsystems/Subsystem', cp, 'Position', p);
delete_line(cp, 'In1/1', 'Out1/1');
delete_block([cp '/In1']);
delete_block([cp '/Out1']);
add_block('simulink/Sources/In1', [cp '/Iapp'], 'Port', '1', ...
          'Position', [p(1) - 90 p(2) - 10 p(1) - 60 p(2) + 10]);
for sp = 1:3
    add_block('simulink/Sources/In1', [cp '/syn' num2str(sp)], ...
        'Port', num2str(sp + 1), ...
        'Position', [p(1) - 90 p(2) + 20 + 22 * sp p(1) - 60 p(2) + 40 + 22 * sp]);
end
add_block('simulink/Sinks/Out1', [cp '/V'], 'Port', '1', ...
          'Position', [p(1) + 190 p(2) p(1) + 220 p(2) + 20]);
add_block('SNS_Library/NonSpikingNeuron', [cp '/cell'], ...
    'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', 'Slope', '1', ...
    'Position', [p(1) + 20 p(2) - 10 p(1) + 120 p(2) + 70]);
add_block('SNS_Library/NonSpikingNeuron', [cp '/adp'], ...
    'Vrest', '-52', 'Gm', '0.1', 'Cm', '750', 'Thr', '-45', 'Slope', '1', ...
    'Position', [p(1) + 20 p(2) + 110 p(1) + 120 p(2) + 190]);
add_block('SNS_Library/NonSpikingSynapse', [cp '/s_adp'], ...
    'gmax', '0.5', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '1', ...
    'ShowName', 'off', 'Position', [p(1) + 150 p(2) + 130 p(1) + 185 p(2) + 160]);
add_line(cp, 'Iapp/1', 'cell/1', 'autorouting', 'on');
add_line(cp, 'cell/1', 'V/1', 'autorouting', 'on');
add_line(cp, 'cell/1', 'adp/1', 'autorouting', 'on');
add_line(cp, 'adp/2', 's_adp/1', 'autorouting', 'on');   % S output, linear
add_line(cp, 's_adp/1', 'cell/2', 'autorouting', 'on');  % cell syn1 = Adp
for sp = 1:3
    add_line(cp, ['syn' num2str(sp) '/1'], ...
             ['cell/' num2str(sp + 2)], 'autorouting', 'on');  % cell syn2..4
end
m = Simulink.Mask.create(cp);
m.Description = ['Li stance CPG v1 (ADAPTATION): graded cell + slow Adp ' ...
    'burst termination, tau = documented tau_Ca 7.5 s. Ports: 1 Iapp ' ...
    '(nA currents), 2-4 synaptic inputs, out = cell V (mV).'];
set_param(cp, 'MaskIconFrame', 'off');
end

