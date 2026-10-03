%% sns_build_spiking.m — ADD spiking blocks to SNS_Library.slx (2026-10-02)
%
% GOAL-2 addition (spiking campaign). ADDITIVE ONLY: this script OPENS the
% existing library and appends three NEW blocks; it never touches the
% committed blocks. Run AFTER sns_build_library.m + sns_build_actuators.m
% (same pattern as sns_build_actuators.m). Re-running it is idempotent: it
% deletes only ITS OWN blocks first.
%
% Doctrine: SNS-Toolbox 1.5.2 spiking semantics (Szczecinski et al. 2017 /
% 2020 generalized-LIF; backends.py SNS_Numpy.forward), the same lineage the
% library header cites. Hybrid pattern from toolbox Tutorials 3/5/9: spikes
% in the interneuron layers, ANALOG voltages into the motoneurons/plant.
%
% THE BLOCKS
%
% 1) SpikingNeuron — spiking LIF with INTERNAL synaptic summation (the
%    committed SpikingLIFNeuron has only an Iapp port; this one is the
%    toolbox-faithful drop-in with syn1..syn6 ports like NonSpikingNeuron).
%      Cm*dV/dt = Gm*(Vrest-V) + Iapp + sum(g*Esyn) - V*sum(g)   [between spikes]
%      V >= Vth  ->  spike output 1, V <- Vreset
%    (toolbox SpikingNeuron with threshold_proportionality_constant 0 =
%     FIXED threshold; the adaptive-threshold variant is future work).
%    Membrane STARTS at Vreset (Integrator reset/initial value, the
%    AnimatLab + toolbox default reset_potential = resting_potential case).
%    Ports: 1 = Iapp [nA], 2..7 = syn1..syn6 ([g; g*Esyn], auto-ground);
%    outputs: V [mV], spike event (0/1, one-step pulse).
%
% 2) SpikingSynapse — spiking chemical synapse, EVENT-wired. Input = the
%    presynaptic neuron's SPIKE line. Conductance (toolbox backends.py:131,178):
%      between spikes: dg/dt = -g/tau_syn  (exact exponential decay)
%      on each presynaptic spike (rising edge): g <- min(gmax, g + ginc)
%    Output = [g; g*Esyn], the SAME 2-wide synapse signal as
%    NonSpikingSynapse -> it plugs into ANY neuron's syn port (SpikingNeuron
%    OR NonSpikingNeuron). This one block is the toolbox SpikingSynapse,
%    conductance_increment included.
%
% 3) HybridSpikingSynapse — the SPIKING-TO-NON-SPIKING synapse of the
%    Tutorials-3/5/9 pattern. Identical conductance dynamics, but the input
%    is the presynaptic V [mV] line (wired exactly like a NonSpikingSynapse)
%    and the synapse detects the spike ITSELF: rising edge of
%    (Vpre >= ThrPre). Works off any neuron's V output (AnimatLab-style);
%    intended target = a NON-SPIKING neuron's syn port, so a hybrid circuit
%    keeps the same look as the non-spiking circuit.
%
% UNITS/convention mapping to sns_toolbox (for cross-platform ports):
%    C[uF] -> Cm = 1000*C [nF]; Gm [uS]; Vrest [mV];
%    toolbox reversal_potential is RELATIVE to the postsynaptic rest:
%    Esyn[block, absolute] = reversal_potential + Vrest_post;
%    Vth[block] = threshold_initial_value + Vrest_post;
%    time_constant [ms] -> tau_syn; conductance_increment [uS] -> ginc.
%    (Vrest_post = 0 in the units-test convention makes both identical.)
%
% Integrator port order (empirically verified on R2025a, goal2_probe3):
%    external IC + rising reset -> port 1 = xdot, port 2 = RESET,
%    port 3 = external INITIAL CONDITION.

lib = 'SNS_Library';
here = fileparts(mfilename('fullpath'));
if ~bdIsLoaded(lib)
    load_system(fullfile(here, [lib '.slx']));
end
set_param(lib, 'Lock', 'off');

% idempotent: clear only our own blocks
own = {'SpikingNeuron', 'SpikingSynapse', 'HybridSpikingSynapse'};
for k = 1:numel(own)
    b = [lib '/' own{k}];
    if ~isempty(find_system(lib, 'SearchDepth', 1, 'Name', own{k}))
        delete_block(b);
    end
end

NSPORT = 6;

%% ---------------- SpikingNeuron ----------------
blk = [lib '/SpikingNeuron'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 2020 130 2110]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Iapp'], 'Port', '1', 'Position', [25 48 55 62]);
add_block('simulink/Math Operations/Sum', [blk '/IappSum'], 'Inputs', '+', 'Position', [95 43 125 67]);
for k = 1:NSPORT
    add_block('simulink/Sources/In1', [blk '/syn' num2str(k)], 'Port', num2str(k+1), ...
        'PortDimensions', '2', 'Position', [25 78+30*(k-1) 55 92+30*(k-1)]);
end
for k = 1:NSPORT
    add_block('simulink/Signal Routing/Demux', [blk '/D' num2str(k)], 'Outputs', '2', ...
        'Position', [95 72+30*(k-1) 98 106+30*(k-1)]);
end
add_block('simulink/Math Operations/Sum', [blk '/sumG'], ...
    'Inputs', repmat('+', 1, NSPORT), 'Position', [170 90 200 90+26*NSPORT]);
add_block('simulink/Math Operations/Sum', [blk '/sumE'], ...
    'Inputs', repmat('+', 1, NSPORT), 'Position', [170 90+26*NSPORT+40 200 90+52*NSPORT+40]);
for k = 1:NSPORT
    add_line(blk, ['syn' num2str(k) '/1'], ['D' num2str(k) '/1'], 'autorouting', 'on');
    add_line(blk, ['D' num2str(k) '/1'], ['sumG/' num2str(k)], 'autorouting', 'on');
    add_line(blk, ['D' num2str(k) '/2'], ['sumE/' num2str(k)], 'autorouting', 'on');
end
% Isyn = sumE - sumG*V
add_block('simulink/Math Operations/Product', [blk '/gTimesV'], 'Inputs', '2', 'Position', [270 95 300 125]);
add_block('simulink/Math Operations/Gain', [blk '/negGv'], 'Gain', '-1', 'Position', [330 100 360 130]);
add_block('simulink/Math Operations/Sum', [blk '/Isyn'], 'Inputs', '++', 'Position', [400 110 430 140]);
% membrane (identical to NonSpikingNeuron) ...
add_block('simulink/Sources/Constant', [blk '/Vrest_c'], 'Value', 'Vrest', 'Position', [330 240 360 270]);
add_block('simulink/Math Operations/Sum', [blk '/dVm'], 'Inputs', '-+', 'Position', [400 245 430 275]);
add_block('simulink/Math Operations/Gain', [blk '/Gm'], 'Gain', 'Gm', 'Position', [460 250 490 280]);
add_block('simulink/Math Operations/Sum', [blk '/dV'], 'Inputs', '+++', 'Position', [520 140 550 170]);
add_block('simulink/Math Operations/Gain', [blk '/membrane'], 'Gain', '1000/Cm', 'Position', [580 140 610 170]);
% ... plus threshold + reset (port 2 = reset; IC param = Vreset).
% SPIKE DETECTION USES THE STATE PORT (outport 2 of membraneRC): the
% threshold->reset wiring is an algebraic loop through the reset port, and
% the state port is Simulink's documented way to break it (LIF pattern).
add_block('simulink/Continuous/Integrator', [blk '/membraneRC'], 'InitialCondition', 'Vreset', ...
    'ExternalReset', 'rising', 'Position', [640 138 670 172]);
set_param([blk '/membraneRC'], 'ShowStatePort', 'on');   % must be set_param, not add_block
add_block('simulink/Logic and Bit Operations/Relational Operator', [blk '/thr'], ...
    'Operator', '>=', 'Position', [700 150 730 180]);
add_block('simulink/Sources/Constant', [blk '/Vth_c'], 'Value', 'Vth', 'Position', [640 205 670 235]);
add_block('simulink/Sinks/Out1', [blk '/V_mV'], 'Port', '1', 'Position', [830 145 860 159]);
add_block('simulink/Sinks/Out1', [blk '/spike'], 'Port', '2', 'Position', [830 210 860 224]);
add_line(blk, 'Iapp/1', 'IappSum/1', 'autorouting', 'on');
add_line(blk, 'IappSum/1', 'dV/1', 'autorouting', 'on');
add_line(blk, 'sumG/1', 'gTimesV/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'gTimesV/2', 'autorouting', 'on');
add_line(blk, 'gTimesV/1', 'negGv/1', 'autorouting', 'on');
add_line(blk, 'sumE/1', 'Isyn/1', 'autorouting', 'on');
add_line(blk, 'negGv/1', 'Isyn/2', 'autorouting', 'on');
add_line(blk, 'Isyn/1', 'dV/2', 'autorouting', 'on');
add_line(blk, 'Vrest_c/1', 'dVm/2', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'dVm/1', 'autorouting', 'on');
add_line(blk, 'dVm/1', 'Gm/1', 'autorouting', 'on');
add_line(blk, 'Gm/1', 'dV/3', 'autorouting', 'on');
add_line(blk, 'dV/1', 'membrane/1', 'autorouting', 'on');
add_line(blk, 'membrane/1', 'membraneRC/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'V_mV/1', 'autorouting', 'on');
phm = get_param([blk '/membraneRC'], 'PortHandles');
pht = get_param([blk '/thr'], 'PortHandles');
add_line(blk, phm.State, pht.Inport(1));                      % STATE port -> detector
add_line(blk, 'Vth_c/1', 'thr/2', 'autorouting', 'on');
add_line(blk, 'thr/1', 'membraneRC/2', 'autorouting', 'on');   % reset (port 2)
add_line(blk, 'thr/1', 'spike/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS Spiking Neuron';
m.Description = ['Spiking leaky integrate-and-fire neuron with INTERNAL synaptic summation ' ...
    '(toolbox SpikingNeuron, fixed-threshold variant). Between spikes the membrane is the ' ...
    'same RC circuit as NonSpikingNeuron: Cm*dV/dt = Gm*(Vrest-V) + Iapp + sum(g*Esyn) - ' ...
    'V*sum(g). When V >= Vth the spike output goes 1 (one-step pulse) and V resets to ' ...
    'Vreset (membrane starts at Vreset). Ports: 1 = Iapp [nA]; syn1..syn6 = synapse lines ' ...
    '[g; g*Esyn] from NonSpikingSynapse / SpikingSynapse / HybridSpikingSynapse (unconnected ' ...
    'ports count as zero). Outputs: V [mV], spike event (0/1). Doctrine: Szczecinski et ' ...
    'al. 2017/2020; sns_toolbox 1.5.2. The adaptive-threshold variant (threshold ' ...
    'proportionality m, threshold increment) is future work.'];
m.Display = spikeNeuronIconCode();
finishMask(m, blk, {'Vrest','Resting potential Vrest (mV)','-52'; ...
                    'Vth','Spike threshold Vth (mV)','-45'; ...
                    'Vreset','Reset (and initial) potential Vreset (mV)','-60'; ...
                    'Gm','Membrane conductance Gm (uS)','0.1'; ...
                    'Cm','Membrane capacitance Cm (nF)','5'});

%% ---------------- SpikingSynapse ----------------
buildSpikeSynapse(lib, 'SpikingSynapse', [40 2160 100 2220], false);
%% ---------------- HybridSpikingSynapse ----------------
buildSpikeSynapse(lib, 'HybridSpikingSynapse', [40 2280 100 2340], true);

save_system(lib);
fprintf('SNS_Library.slx: added SpikingNeuron + SpikingSynapse + HybridSpikingSynapse (14 blocks total).\n');

%% ---------------- local functions ----------------
function buildSpikeSynapse(lib, name, pos, hybrid)
% One spike-event-driven synapse. hybrid=false: input = spike line.
% hybrid=true: input = Vpre [mV], spike detected at Vpre >= ThrPre.
blk = [lib '/' name];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', pos);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
if hybrid
    add_block('simulink/Sources/In1', [blk '/Vpre'], 'Port', '1', 'Position', [15 33 40 47]);
    add_block('simulink/Logic and Bit Operations/Relational Operator', [blk '/xing'], ...
        'Operator', '>=', 'Position', [55 25 85 55]);
    add_block('simulink/Sources/Constant', [blk '/ThrPre_c'], 'Value', 'ThrPre', 'Position', [15 63 45 93]);
    add_line(blk, 'Vpre/1', 'xing/1', 'autorouting', 'on');
    add_line(blk, 'ThrPre_c/1', 'xing/2', 'autorouting', 'on');
    ev = 'xing';     % event node source
else
    add_block('simulink/Sources/In1', [blk '/spike'], 'Port', '1', 'Position', [15 33 40 47]);
    ev = 'spike';
end
% conductance state g: xdot = -g*1000/tau_syn, reset (port 2) = event,
% external IC (port 3) = min(gmax, g + ginc*event).
% The IC branch reads the STATE PORT (gcore outport 2): feeding the IC port
% from the OUTPUT port would close an algebraic loop through the reset/IC
% path; the state port is Simulink's documented loop-breaker.
add_block('simulink/Math Operations/Gain', [blk '/leak'], 'Gain', '-1000/tau_syn', 'Position', [120 108 155 142]);
add_block('simulink/Continuous/Integrator', [blk '/gcore'], 'InitialConditionSource', 'external', ...
    'ExternalReset', 'rising', 'Position', [185 100 220 150]);
set_param([blk '/gcore'], 'ShowStatePort', 'on');
add_block('simulink/Math Operations/Gain', [blk '/gincG'], 'Gain', 'ginc', ...
    'OutDataTypeStr', 'double', 'Position', [120 165 155 195]);
add_block('simulink/Math Operations/Sum', [blk '/gPlus'], 'Inputs', '++', 'Position', [245 165 275 195]);
add_block('simulink/Discontinuities/Saturation', [blk '/gcap'], ...
    'UpperLimit', 'gmax', 'LowerLimit', '0', 'Position', [305 165 335 195]);
add_line(blk, 'gcore/1', 'leak/1', 'autorouting', 'on');
add_line(blk, 'leak/1', 'gcore/1', 'autorouting', 'on');          % xdot  (port 1)
add_line(blk, [ev '/1'], 'gcore/2', 'autorouting', 'on');         % reset (port 2)
add_line(blk, [ev '/1'], 'gincG/1', 'autorouting', 'on');
phg = get_param([blk '/gcore'], 'PortHandles');
phs = get_param([blk '/gPlus'], 'PortHandles');
add_line(blk, phg.State, phs.Inport(1));                         % STATE port -> IC feed
add_line(blk, 'gincG/1', 'gPlus/2', 'autorouting', 'on');
add_line(blk, 'gPlus/1', 'gcap/1', 'autorouting', 'on');
add_line(blk, 'gcap/1', 'gcore/3', 'autorouting', 'on');          % ext IC (port 3)
% output pair [g; g*Esyn]
add_block('simulink/Math Operations/Gain', [blk '/EsynG'], 'Gain', 'Esyn', 'Position', [305 60 340 90]);
add_block('simulink/Signal Routing/Mux', [blk '/pair'], 'Inputs', '2', 'Position', [390 90 393 140]);
add_block('simulink/Sinks/Out1', [blk '/syn_out'], 'Position', [430 107 460 121]);
add_line(blk, 'gcore/1', 'EsynG/1', 'autorouting', 'on');
add_line(blk, 'gcore/1', 'pair/1', 'autorouting', 'on');
add_line(blk, 'EsynG/1', 'pair/2', 'autorouting', 'on');
add_line(blk, 'pair/1', 'syn_out/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
if hybrid
    m.Type = 'SNS Hybrid Spiking Synapse';
    desc = ['HYBRID spiking-to-non-spiking synapse (toolbox Tutorials 3/5/9 pattern: spikes ' ...
        'above, analog voltages into the plant). Input: presynaptic membrane potential Vpre ' ...
        '[mV] from ANY neuron (spiking or not) — the synapse detects the spike itself: ' ...
        'rising edge of (Vpre >= ThrPre). Conductance: exponential decay with tau_syn; on ' ...
        'each detected spike g <- min(gmax, g + ginc). Output [g; g*Esyn] connects to the ' ...
        'syn1..syn6 port of a NON-SPIKING (or spiking) postsynaptic neuron. ' ...
        'CHOOSING ThrPre (AnimatLab semantics): it must sit BETWEEN the presynaptic ' ...
        'Vrest and Vth, comfortably BELOW Vth (a spiking presynaptic cell RESETS at ' ...
        'threshold, so its V output never renders values >= Vth — ThrPre = Vth fires ' ...
        'nothing; ThrPre ~ 0.9*Vth keys the g jump ~1 ms before the true spike).'];
else
    m.Type = 'SNS Spiking Synapse';
    desc = ['Spiking chemical synapse, event-wired: input = the SPIKE output line of a ' ...
        'SpikingNeuron / SpikingLIFNeuron. Conductance (sns_toolbox SpikingSynapse): ' ...
        'exponential decay with tau_syn; on each presynaptic spike (rising edge) ' ...
        'g <- min(gmax, g + ginc) — conductance_increment included, saturating at gmax. ' ...
        'Output [g; g*Esyn] connects to ANY neuron''s syn1..syn6 port (spiking or ' ...
        'non-spiking postsynaptic cell). Esyn >= 0 excitatory (icon: white triangle), ' ...
        '< 0 inhibitory (solid black circle).'];
end
m.Description = desc;
m.Display = spikeSynapseIconCode();
if hybrid
    finishMask(m, blk, {'gmax','Max conductance gmax (uS)','1'; ...
                        'ginc','Conductance increment per spike (uS)','1'; ...
                        'tau_syn','Decay time constant tau_syn (ms)','20'; ...
                        'Esyn','Reversal potential Esyn (mV). >=0 excit, <0 inhib','0'; ...
                        'ThrPre','Spike-detection threshold ThrPre (mV)','-45'});
else
    finishMask(m, blk, {'gmax','Max conductance gmax (uS)','1'; ...
                        'ginc','Conductance increment per spike (uS)','1'; ...
                        'tau_syn','Decay time constant tau_syn (ms)','20'; ...
                        'Esyn','Reversal potential Esyn (mV). >=0 excit, <0 inhib','0'});
end
end

function finishMask(m, blk, params)
    for k = 1:size(params, 1)
        m.addParameter('Name', params{k,1}, 'Type', 'edit', ...
            'Prompt', params{k,2}, 'Value', params{k,3});
    end
    set_param(blk, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
        'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');
end

function s = spikeNeuronIconCode()
    % Open circle (heavy outline) + SPIKE-TRAIN waveform (same family as the
    % committed SpikingLIFNeuron icon; Ben's 2026-09-22 spike-train ask).
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.95*cos(t_), 0.95*sin(t_), [0 0 0]);' ...
        'patch(0.84*cos(t_), 0.84*sin(t_), [0.99 0.96 0.81]);' ...
        'xs_ = [-0.58 -0.54 -0.50 -0.46 -0.36 -0.32 -0.28 -0.24 -0.14 -0.10 -0.06 -0.02 0.08 0.12 0.16 0.20 0.30 0.34 0.38 0.42 0.52 0.58];' ...
        'ys_ = [-0.22 -0.22 0.40 -0.22 -0.22 -0.22 0.40 -0.22 -0.22 -0.22 0.40 -0.22 -0.22 -0.22 0.40 -0.22 -0.22 -0.22 0.40 -0.22 -0.22 -0.22];' ...
        'color(''black'');' ...
        'plot(xs_, ys_);' ...
        'plot(xs_+0.012, ys_);' ...
        'plot(xs_-0.012, ys_);' ...
        }, newline);
end

function s = spikeSynapseIconCode()
    % Small axon bar carrying a SPIKE TRAIN, terminating in the postsynaptic
    % marker picked automatically from the sign of Esyn (same convention as
    % NonSpikingSynapse: white triangle = excitatory, black circle = inhibitory).
    s = strjoin({ ...
        'es_ = NaN;' ...
        'try' ...
        '    tmp_ = Esyn;' ...
        '    if ischar(tmp_) || isstring(tmp_), tmp_ = eval(tmp_); end' ...
        '    es_ = tmp_;' ...
        'catch' ...
        'end' ...
        'if isnan(es_)' ...
        '    patch([-1 1 1 -1], [-1 -1 1 1], [0.9 0.9 0.9]);' ...
        '    disp(''E?'');' ...
        'elseif es_ < 0' ...
        '    patch([-1 1 1 -1], [-1 -1 1 1], [0.97 0.86 0.85]);' ...
        '    patch([-1 0.58 0.58 -1], [0.42 0.42 0.58 0.58], [0 0 0]);' ...
        '    t_ = linspace(0, 2*pi, 73);' ...
        '    patch(0.18*cos(t_) + 0.72, 0.50*sin(t_), [0 0 0]);' ...
        '    color(''black'');' ...
        '    plot([-1 -0.75 -0.65 -0.55 -0.35 -0.25 -0.15 0.05 0.15 0.25 0.45 0.55 0.58], [0.50 0.50 0.72 0.50 0.50 0.72 0.50 0.50 0.72 0.50 0.50 0.50 0.50]);' ...
        'else' ...
        '    patch([-1 1 1 -1], [-1 -1 1 1], [0.85 0.94 0.86]);' ...
        '    patch([-1 0.30 0.30 -1], [0.42 0.42 0.58 0.58], [0 0 0]);' ...
        '    patch([0.30 0.90 0.90], [0.50 0.80 -0.20], [0 0 0]);' ...
        '    patch([0.44 0.76 0.76], [0.50 0.70 -0.10], [1 1 1]);' ...
        '    color(''black'');' ...
        '    plot([-1 -0.75 -0.65 -0.55 -0.35 -0.25 -0.15 0.05 0.15 0.25 0.30], [0.50 0.50 0.72 0.50 0.50 0.72 0.50 0.50 0.72 0.50 0.50]);' ...
        'end' ...
        }, newline);
end
