%% sns_build_library.m — build SNS_Library.slx
% Synthetic Nervous System block library for the Bipedal_Robot Simscape/Simulink work.
%
% ARCHITECTURE (Ben's 2026-09-22 redesign):
%   * NonSpikingSynapse = ONE INPUT (Vpre) -> ONE OUTPUT. The synapse outputs a
%     2-wide synaptic signal [g; g*Esyn] (uS and uS*mV); the POSTSYNAPTIC
%     neuron evaluates the driving force itself. No Vpost sense wire, no
%     postsynaptic input port on the synapse.
%   * NonSpikingNeuron sums its synaptic inputs INTERNALLY: port 1 = injected
%     current Iapp [nA] (descending drive / afferent currents; a vector input
%     is summed element-wise), ports 2..7 = synapse inputs [g; g*Esyn].
%     Unconnected ports auto-ground to zero (model diagnostic
%     UnconnectedInputMsg = 'none', the default) — NO Sum block in front of
%     the neuron, and a neuron with fewer than 6 synapses needs no dummy
%     wiring. Isyn = sum(g*Esyn) - V*sum(g) — algebraically identical to
%     sum_k gmax_k*Sat(Vpre_k)*(Esyn_k - V), so all previous tuning and the
%     numpy units reference remain valid.
%   * SynSum = junction for generated big models (e.g. SNS_SpinalNetwork):
%     sums up to EIGHT synapse signals into one [sum(g); sum(g*Esyn)] line so
%     a neuron receiving >6 synapses needs only one wire per SynSum.
%
% Neuron model = non-spiking leaky integrator, literally an RC membrane:
%   Cm dV/dt = Gm*(Vrest - V) + Iapp + Isyn
%   tau_m = Cm/Gm  (nF/uS = ms). Same convention as Animatlab non-spiking
%   neurons and sns_toolbox NonSpikingNeuron (Szczecinski et al. 2017).
% Synapse model (Animatlab + SNS toolbox non-spiking chemical synapse):
%   g(Vpre) = gmax * Sat(Vpre),  Sat = clip((Vpre - ThrPre)/SlopePre, 0, 1)
%   Isyn    = g*(Esyn - Vpost)  [uS*mV = nA]  (evaluated in the neuron)
%   Esyn > Vrest  -> excitatory (depolarizing), Esyn < Vrest -> inhibitory.
%
% --- APPEARANCE CONVENTIONS (Ben, 2026-09-09 + 2026-09-22) ---------------------
%   * NON-SPIKING neuron  = open circle + GRADED-POTENTIAL WAVEFORM (smooth hump)
%   * SPIKING neuron      = open circle + SPIKE-TRAIN WAVEFORM
%   * Ia muscle spindle   = SPINDLE-SHAPED capsule (fusiform, tapered ends)
%   * Ib Golgi tendon     = capsule with braided collagen strands, "Ib"
%   * muscle / BPA blocks = fusiform WITH STRIATIONS across the belly
%   * MuscleActivation    = PENTAGON (pink)
%   * synapse             = small pass-through axon bar terminating in the
%     postsynaptic marker at its output edge:
%       Esyn <  0 -> inhibitory  -> SOLID BLACK CIRCLE
%       Esyn >= 0 -> excitatory  -> WHITE TRIANGLE, heavy black edge,
%                                   tip pointing back toward the presynaptic side
%   Shape is the primary code; tints are the redundant CVD-safe cue.
% THICK STROKES: Simulink mask icons cannot set LineWidth, so outlines are
% drawn as filled shapes (black outer shape + smaller filled inner shape).
% Keep CIRCLE-icon blocks SQUARE (width == height) so circles stay round.

lib = 'SNS_Library';
if bdIsLoaded(lib), close_system(lib, 0); end
if exist([lib '.slx'], 'file'), delete([lib '.slx']); end
new_system(lib, 'Library');
load_system(lib);
NSPORT = 6;   % synapse input ports on NonSpikingNeuron (ports 2..NSPORT+1)

%% ---------------- NonSpikingNeuron (internal synaptic summation) ----------------
blk = [lib '/NonSpikingNeuron'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 40 130 130]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Iapp'], 'Port', '1', 'Position', [25 48 55 62]);
add_block('simulink/Math Operations/Sum', [blk '/IappSum'], 'Inputs', '+', 'Position', [95 43 125 67]);  % scalar or vector -> scalar
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
% Isyn = sumE - sumG*V  (membrane V feeds back from the integrator)
add_block('simulink/Math Operations/Product', [blk '/gTimesV'], 'Inputs', '2', 'Position', [270 95 300 125]);
add_block('simulink/Math Operations/Gain', [blk '/negGv'], 'Gain', '-1', 'Position', [330 100 360 130]);
add_block('simulink/Math Operations/Sum', [blk '/Isyn'], 'Inputs', '++', 'Position', [400 110 430 140]);
% membrane: dV = Gm*(Vrest - V) + Iapp + Isyn, V' = 1000/Cm * dV
add_block('simulink/Sources/Constant', [blk '/Vrest_c'], 'Value', 'Vrest', 'Position', [330 240 360 270]);
add_block('simulink/Math Operations/Sum', [blk '/dVm'], 'Inputs', '-+', 'Position', [400 245 430 275]);
add_block('simulink/Math Operations/Gain', [blk '/Gm'], 'Gain', 'Gm', 'Position', [460 250 490 280]);
add_block('simulink/Math Operations/Sum', [blk '/dV'], 'Inputs', '+++', 'Position', [520 140 550 170]);
add_block('simulink/Math Operations/Gain', [blk '/membrane'], 'Gain', '1000/Cm', 'Position', [580 140 610 170]);
add_block('simulink/Continuous/Integrator', [blk '/membraneRC'], 'InitialCondition', 'Vrest', 'Position', [640 138 670 172]);
% S(V) drive output
add_block('simulink/Math Operations/Sum', [blk '/thrSub'], 'Inputs', '+-', 'Position', [640 220 670 250]);
add_block('simulink/Sources/Constant', [blk '/Thr_c'], 'Value', 'Thr', 'Position', [560 255 590 285]);
add_block('simulink/Math Operations/Gain', [blk '/invSlope'], 'Gain', '1/Slope', 'Position', [700 220 730 250]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', '1', 'LowerLimit', '0', 'Position', [760 220 790 250]);
add_block('simulink/Sinks/Out1', [blk '/V_mV'], 'Port', '1', 'Position', [830 145 860 159]);
add_block('simulink/Sinks/Out1', [blk '/S_drive'], 'Port', '2', 'Position', [830 225 860 239]);
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
add_line(blk, 'membraneRC/1', 'thrSub/1', 'autorouting', 'on');
add_line(blk, 'Thr_c/1', 'thrSub/2', 'autorouting', 'on');
add_line(blk, 'thrSub/1', 'invSlope/1', 'autorouting', 'on');
add_line(blk, 'invSlope/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'S_drive/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS NonSpiking Neuron';
m.Description = ['Non-spiking leaky integrate-and-fire neuron (RC membrane) with INTERNAL ' ...
    'synaptic summation. Port 1 (Iapp): injected current [nA] — descending drive or afferent ' ...
    'currents (a vector input is summed element-wise). Ports syn1..syn6: synaptic inputs from ' ...
    'NonSpikingSynapse blocks ([g; g*Esyn]); each synapse connects to the NEXT FREE syn port, ' ...
    'unconnected ports count as zero (leave them empty). ' ...
    'Cm*dV/dt = Gm*(Vrest-V) + Iapp + sum(g*Esyn) - V*sum(g). ' ...
    'Outputs: membrane potential V [mV] and normalized drive S(V) [0..1]. ' ...
    'Conventions: Animatlab / sns_toolbox (Szczecinski et al. 2017).'];
m.Display = neuronIconCode();
finishMask(m, blk, {'Vrest','Resting potential Vrest (mV)','-52'; ...
                    'Gm','Membrane conductance Gm (uS)','0.1'; ...
                    'Cm','Membrane capacitance Cm (nF)','5'; ...
                    'Thr','Saturation threshold Thr (mV)','-55'; ...
                    'Slope','Saturation slope (1/mV)','1'});

%% ---------------- NonSpikingSynapse (1 input -> 1 output) ----------------
blk = [lib '/NonSpikingSynapse'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 220 100 280]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Vpre'], 'Port', '1', 'Position', [15 53 40 67]);
add_block('simulink/Math Operations/Sum', [blk '/thrSub'], 'Inputs', '+-', 'Position', [60 48 90 78]);
add_block('simulink/Sources/Constant', [blk '/ThrPre_c'], 'Value', 'ThrPre', 'Position', [15 95 45 125]);
add_block('simulink/Math Operations/Gain', [blk '/invSlope'], 'Gain', '1/SlopePre', 'Position', [110 48 140 78]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', '1', 'LowerLimit', '0', 'Position', [160 48 190 78]);
add_block('simulink/Math Operations/Gain', [blk '/gmax_g'], 'Gain', 'gmax', 'Position', [210 48 240 78]);
% g*Esyn branch
add_block('simulink/Math Operations/Gain', [blk '/gEsyn_g'], 'Gain', 'gmax*Esyn', 'Position', [210 108 250 138]);
add_block('simulink/Signal Routing/Mux', [blk '/pair'], 'Inputs', '2', 'Position', [290 68 293 118]);
add_block('simulink/Sinks/Out1', [blk '/syn_out'], 'Position', [340 85 370 99]);
add_line(blk, 'Vpre/1', 'thrSub/1', 'autorouting', 'on');
add_line(blk, 'ThrPre_c/1', 'thrSub/2', 'autorouting', 'on');
add_line(blk, 'thrSub/1', 'invSlope/1', 'autorouting', 'on');
add_line(blk, 'invSlope/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'gmax_g/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'gEsyn_g/1', 'autorouting', 'on');
add_line(blk, 'gmax_g/1', 'pair/1', 'autorouting', 'on');
add_line(blk, 'gEsyn_g/1', 'pair/2', 'autorouting', 'on');
add_line(blk, 'pair/1', 'syn_out/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS NonSpiking Synapse';
m.Description = ['Non-spiking chemical synapse, ONE input -> ONE output. Input: presynaptic ' ...
    'membrane potential Vpre [mV]. Output: 2-wide synaptic signal [g; g*Esyn] that connects ' ...
    'to a NonSpikingSynapse input port (syn1..syn6) of the POSTSYNAPTIC NonSpikingNeuron — ' ...
    'place the synapse next to the neuron it synapses onto. ' ...
    'g = gmax*Sat(Vpre), Sat = clip((Vpre-ThrPre)/SlopePre, 0, 1); the neuron forms ' ...
    'Isyn = sum(g*Esyn) - Vpost*sum(g). Esyn >= 0 mV -> EXCITATORY (icon: white triangle); ' ...
    'Esyn < 0 mV -> INHIBITORY (icon: solid black circle). Diagram language: ' ...
    'Szczecinski et al. 2017 Fig. 2.'];
m.Display = synapseIconCode();
finishMask(m, blk, {'gmax','Max conductance gmax (uS)','1'; ...
                    'Esyn','Reversal potential Esyn (mV). >=0 excit, <0 inhib','0'; ...
                    'ThrPre','Pre saturation threshold ThrPre (mV)','-55'; ...
                    'SlopePre','Pre saturation slope SlopePre (1/mV)','1'});

%% ---------------- SynSum (junction for >6 synapses onto one neuron) -----------
blk = [lib '/SynSum'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 340 96 396]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
NSYN = 8;
for k = 1:NSYN
    add_block('simulink/Sources/In1', [blk '/in' num2str(k)], 'Port', num2str(k), ...
        'PortDimensions', '2', 'Position', [20 28+26*(k-1) 45 42+26*(k-1)]);
    add_block('simulink/Signal Routing/Demux', [blk '/Dm' num2str(k)], 'Outputs', '2', ...
        'Position', [80 22+26*(k-1) 83 56+26*(k-1)]);
    add_line(blk, ['in' num2str(k) '/1'], ['Dm' num2str(k) '/1'], 'autorouting', 'on');
end
add_block('simulink/Math Operations/Sum', [blk '/sumG'], 'Inputs', repmat('+', 1, NSYN), ...
    'Position', [150 60 180 60+24*NSYN]);
add_block('simulink/Math Operations/Sum', [blk '/sumE'], 'Inputs', repmat('+', 1, NSYN), ...
    'Position', [150 60+24*NSYN+40 180 60+48*NSYN+40]);
for k = 1:NSYN
    add_line(blk, ['Dm' num2str(k) '/1'], ['sumG/' num2str(k)], 'autorouting', 'on');
    add_line(blk, ['Dm' num2str(k) '/2'], ['sumE/' num2str(k)], 'autorouting', 'on');
end
add_block('simulink/Signal Routing/Mux', [blk '/pair'], 'Inputs', '2', 'Position', [230 90 233 140]);
add_block('simulink/Sinks/Out1', [blk '/syn_out'], 'Position', [270 107 300 121]);
add_line(blk, 'sumG/1', 'pair/1', 'autorouting', 'on');
add_line(blk, 'sumE/1', 'pair/2', 'autorouting', 'on');
add_line(blk, 'pair/1', 'syn_out/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS Synapse Sum Junction';
m.Description = ['Sums up to eight NonSpikingSynapse outputs ([g; g*Esyn] lines) into ONE ' ...
    'synaptic signal for a NonSpikingNeuron syn port. Only needed when a neuron receives ' ...
    'MORE than 6 synapses (the neuron has 6 syn ports): chain SynSum blocks for even more ' ...
    '(a SynSum output can feed another SynSum input). Unconnected inputs count as zero.'];
m.Display = synsumIconCode();
finishMask(m, blk, {});

%% ---------------- SpikingLIFNeuron ----------------
blk = [lib '/SpikingLIFNeuron'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 480 130 570]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Iapp'], 'Port', '1', 'Position', [25 78 55 92]);
add_block('simulink/Math Operations/Sum', [blk '/IappSum'], 'Inputs', '+', 'Position', [95 73 125 97]);
add_block('simulink/Sources/Constant', [blk '/Vrest_c'], 'Value', 'Vrest', 'Position', [90 150 120 180]);
add_block('simulink/Math Operations/Sum', [blk '/sumI'], 'Inputs', '++', 'Position', [250 70 280 100]);
add_block('simulink/Math Operations/Sum', [blk '/leak'], 'Inputs', '-+', 'Position', [150 108 180 138]);
add_block('simulink/Math Operations/Gain', [blk '/Rm_g'], 'Gain', 'Rm', 'Position', [200 113 230 143]);
add_block('simulink/Math Operations/Gain', [blk '/tauInv'], 'Gain', '1000/tau', 'Position', [310 70 340 100]);
add_block('simulink/Continuous/Integrator', [blk '/Vmem'], 'InitialCondition', 'Vrest', 'ExternalReset', 'rising', 'Position', [370 63 400 107]);
add_block('simulink/Logic and Bit Operations/Relational Operator', [blk '/thr'], 'Operator', '>=', 'Position', [430 130 460 160]);
add_block('simulink/Sources/Constant', [blk '/Vth_c'], 'Value', 'Vth', 'Position', [370 175 400 205]);
add_block('simulink/Sinks/Out1', [blk '/V_mV'], 'Port', '1', 'Position', [560 78 590 92]);
add_block('simulink/Sinks/Out1', [blk '/spike'], 'Port', '2', 'Position', [560 138 590 152]);
add_line(blk, 'Iapp/1', 'IappSum/1', 'autorouting', 'on');
add_line(blk, 'IappSum/1', 'sumI/2', 'autorouting', 'on');
add_line(blk, 'Vrest_c/1', 'leak/2', 'autorouting', 'on');
add_line(blk, 'leak/1', 'Rm_g/1', 'autorouting', 'on');
add_line(blk, 'Rm_g/1', 'sumI/1', 'autorouting', 'on');
add_line(blk, 'sumI/1', 'tauInv/1', 'autorouting', 'on');
add_line(blk, 'tauInv/1', 'Vmem/1', 'autorouting', 'on');
add_line(blk, 'Vmem/1', 'leak/1', 'autorouting', 'on');
add_line(blk, 'Vmem/1', 'thr/1', 'autorouting', 'on');
add_line(blk, 'Vth_c/1', 'thr/2', 'autorouting', 'on');
add_line(blk, 'Vmem/1', 'V_mV/1', 'autorouting', 'on');
add_line(blk, 'thr/1', 'Vmem/2', 'autorouting', 'on');
add_line(blk, 'thr/1', 'spike/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS Spiking LIF Neuron';
m.Description = ['Spiking leaky integrate-and-fire neuron: tau*dV/dt = (Vrest-V) + Rm*Iapp; ' ...
    'V>=Vth -> output spike flag and membrane resets to Vrest. Input: injected current ' ...
    'Iapp [nA] (vector inputs are summed element-wise).'];
m.Display = lifIconCode();
finishMask(m, blk, {'Vrest','Resting potential Vrest (mV)','-60'; ...
                    'Vth','Spike threshold Vth (mV)','-50'; ...
                    'Rm','Membrane resistance Rm (MOhm)','10'; ...
                    'tau','Membrane time constant tau (ms)','20'});

%% ---------------- IaMuscleSpindle ----------------
blk = [lib '/IaMuscleSpindle'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 640 130 730]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/stretch'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Sources/In1', [blk '/vel'], 'Port', '2', 'Position', [25 138 55 152]);
add_block('simulink/Math Operations/Gain', [blk '/Wl_g'], 'Gain', 'Wl', 'Position', [100 55 130 85]);
add_block('simulink/Math Operations/Gain', [blk '/Wv_g'], 'Gain', 'Wv', 'Position', [100 130 130 160]);
add_block('simulink/Math Operations/Sum', [blk '/sumIa'], 'Inputs', '++', 'Position', [170 90 200 120]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', 'Imax', 'LowerLimit', '0', 'Position', [240 90 270 120]);
add_block('simulink/Sinks/Out1', [blk '/Isyn_nA'], 'Position', [310 98 340 112]);
add_line(blk, 'stretch/1', 'Wl_g/1', 'autorouting', 'on');
add_line(blk, 'vel/1', 'Wv_g/1', 'autorouting', 'on');
add_line(blk, 'Wl_g/1', 'sumIa/1', 'autorouting', 'on');
add_line(blk, 'Wv_g/1', 'sumIa/2', 'autorouting', 'on');
add_line(blk, 'sumIa/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'Isyn_nA/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS Ia Muscle Spindle Afferent';
m.Description = 'Ia spindle afferent: Isyn = clip(Wl*stretch + Wv*velocity, 0, Imax) [nA]. Inputs are normalized (0..1) length and velocity signals. Icon: spindle-shaped capsule.';
m.Display = spindleIconCode('Ia');
finishMask(m, blk, {'Imax','Peak current Imax (nA)','10'; ...
                    'Wl','Length weight Wl (nA)','6'; ...
                    'Wv','Velocity weight Wv (nA)','8'});

%% ---------------- IbGolgiTendon ----------------
blk = [lib '/IbGolgiTendon'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 800 130 890]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/force'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Math Operations/Gain', [blk '/Kf_g'], 'Gain', 'Kf', 'Position', [110 55 140 85]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', 'Imax', 'LowerLimit', '0', 'Position', [180 55 210 85]);
add_block('simulink/Sinks/Out1', [blk '/Isyn_nA'], 'Position', [250 63 280 77]);
add_line(blk, 'force/1', 'Kf_g/1', 'autorouting', 'on');
add_line(blk, 'Kf_g/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'Isyn_nA/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS Ib Golgi Tendon Afferent';
m.Description = 'Ib Golgi tendon organ afferent: Isyn = clip(Kf*forceNorm, 0, Imax) [nA]. Input is normalized force (0..1). Icon: capsule with braided collagen strands.';
m.Display = spindleIconCode('Ib');
finishMask(m, blk, {'Imax','Peak current Imax (nA)','10'; ...
                    'Kf','Force weight Kf (nA)','10'});

%% ---------------- MuscleActivation (PENTAGON, pink) ----------------
blk = [lib '/MuscleActivation'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 960 130 1050]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/S'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Math Operations/Sum', [blk '/dA'], 'Inputs', '-+', 'Position', [90 55 120 85]);
add_block('simulink/Math Operations/Gain', [blk '/tauInv'], 'Gain', '1000/tauAct', 'Position', [150 55 180 85]);
add_block('simulink/Continuous/Integrator', [blk '/Aint'], 'InitialCondition', 'A0', 'Position', [210 55 240 85]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', '1', 'LowerLimit', '0', 'Position', [270 55 300 85]);
add_block('simulink/Sinks/Out1', [blk '/activation'], 'Position', [340 63 370 77]);
add_line(blk, 'S/1', 'dA/2', 'autorouting', 'on');
add_line(blk, 'Aint/1', 'dA/1', 'autorouting', 'on');
add_line(blk, 'dA/1', 'tauInv/1', 'autorouting', 'on');
add_line(blk, 'tauInv/1', 'Aint/1', 'autorouting', 'on');
add_line(blk, 'Aint/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'activation/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS Muscle Activation';
m.Description = ['First-order activation dynamics: tauAct*dA/dt = (S - A). Maps motoneuron ' ...
    'drive S(V) [0..1] to muscle activation A [0..1]. A0 = initial activation (set it to ' ...
    'the steady hold activation so a model starts settled, with no startup transient).'];
m.Display = pentagonIconCode();
finishMask(m, blk, {'tauAct','Activation time constant tauAct (ms)','50'; ...
                    'A0','Initial activation A(0) [0..1]','0'});

%% ---------------- BPAForce ----------------
blk = [lib '/BPAForce'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 1120 130 1210]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/activation'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Sources/In1', [blk '/strain'], 'Port', '2', 'Position', [25 138 55 152]);
add_block('simulink/Math Operations/Gain', [blk '/Fmax_g'], 'Gain', 'Fmax', 'Position', [100 55 130 85]);
add_block('simulink/Math Operations/Sum', [blk '/Lterm'], 'Inputs', '+-', 'Position', [100 130 130 160]);
add_block('simulink/Sources/Constant', [blk '/epsMax_c'], 'Value', 'epsMax', 'Position', [30 175 60 205]);
add_block('simulink/Math Operations/Product', [blk '/Fprod'], 'Inputs', '2', 'Position', [200 85 230 115]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', 'Fmax', 'LowerLimit', '0', 'Position', [260 85 290 115]);
add_block('simulink/Sinks/Out1', [blk '/force_N'], 'Position', [330 93 360 107]);
% Fmax_g already forms Fmax*A; apply activation exactly once.
add_line(blk, 'activation/1', 'Fmax_g/1', 'autorouting', 'on');
add_line(blk, 'epsMax_c/1', 'Lterm/1', 'autorouting', 'on');
add_line(blk, 'strain/1', 'Lterm/2', 'autorouting', 'on');
add_line(blk, 'Fmax_g/1', 'Fprod/1', 'autorouting', 'on');
add_line(blk, 'Lterm/1', 'Fprod/2', 'autorouting', 'on');
add_line(blk, 'Fprod/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'force_N/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS BPA Force Element';
m.Description = 'Simple BPA (pneumatic artificial muscle) force: F = Fmax*A*max(0, epsMax - strain). Replace with MonoPam/Xi-corrected prediction when available.';
m.Display = muscleIconCode('BPA');
finishMask(m, blk, {'Fmax','Max isometric force Fmax (N)','500'; ...
                    'epsMax','Max contraction strain epsMax (0..1)','0.25'});

save_system(lib);
fprintf('SNS_Library.slx saved with 8 library blocks (2026-09-22 architecture).\n');

%% ---------------- local functions ----------------
% NOTE on mask icon code: Simulink mask drawing commands (plot/patch/text) are
% a RESTRICTED subset — numeric arguments only, no LineSpec strings ('k-'), no
% name-value pairs ('LineWidth',...), and color() takes a color NAME only.
% Use patch(...,[r g b]) for fills, color('black') before plot for edges,
% and disp() for centered text.
function finishMask(m, blk, params)
    for k = 1:size(params, 1)
        m.addParameter('Name', params{k,1}, 'Type', 'edit', ...
            'Prompt', params{k,2}, 'Value', params{k,3});
    end
    % Square blocks + autoscaled mask icons -> circles stay circular.
    % Frame off = the drawn shape IS the block, like a circuit diagram.
    set_param(blk, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
        'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');
end

function s = neuronIconCode()
    % Open circle (heavy outline) + GRADED-POTENTIAL waveform: a non-spiking
    % neuron signals with smoothly varying membrane potential, so its icon
    % shows a smooth depolarizing hump (Ben 2026-09-22).
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.95*cos(t_), 0.95*sin(t_), [0 0 0]);' ...
        'patch(0.84*cos(t_), 0.84*sin(t_), [0.99 0.96 0.81]);' ...
        'x_ = linspace(-0.55, 0.55, 61);' ...
        'y_ = -0.20 + 0.40*exp(-((x_+0.08)/0.26).^2);' ...
        'color(''black'');' ...
        'plot(x_, y_);' ...
        'plot([-0.55 0.55], [-0.20 -0.20]);' ...
        }, newline);
end

function s = lifIconCode()
    % Open circle (heavy outline) + SPIKE-TRAIN waveform (4 action potentials).
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

function s = spindleIconCode(lbl)
    % SPINDLE-SHAPED capsule (fusiform, tapered at both ends) — the block
    % itself looks like a muscle spindle. Green-gray afferent tint; wavy
    % intrafusal fiber strand below the label.
    s = strjoin({ ...
        'x_ = linspace(-1, 1, 41);' ...
        'w_ = 0.58*(1 - abs(x_).^1.5);' ...
        'patch([x_ fliplr(x_)], [w_ -fliplr(w_)], [0 0 0]);' ...
        'x2_ = 0.93*x_;' ...
        'w2_ = 0.58*(1 - abs(x2_).^1.5);' ...
        'patch([x2_ fliplr(x2_)], [w2_ -fliplr(w2_)], [0.85 0.93 0.85]);' ...
        'color(''black'');' ...
        'plot(linspace(-0.55, 0.55, 27), 0.05*sin(linspace(0, 5*pi, 27)) - 0.24);' ...
        'plot(linspace(-0.55, 0.55, 27), 0.05*sin(linspace(0, 5*pi, 27)) - 0.24 + 0.06);' ...
        ['disp(''' lbl ''');'] ...
        }, newline);
end

function s = pentagonIconCode()
    % PENTAGON, pink — the muscle-activation stage (Ben 2026-09-22: "make
    % them a pentagon shape, or maybe the gain triangles but pink").
    a_ = [0, 1.2566, 2.5133, 3.7699, 5.0265];
    s = strjoin({ ...
        ['patch(0.95*cos([' num2str(a_) ']), 0.95*sin([' num2str(a_) ']), [0 0 0]);'] ...
        ['patch(0.82*cos([' num2str(a_) ']), 0.82*sin([' num2str(a_) ']), [0.96 0.72 0.78]);'] ...
        }, newline);
end

function s = muscleIconCode(lbl)
    % Fusiform muscle with heavy outline AND STRIATIONS across the belly
    % (Ben 2026-09-22: "muscle blocks should have striations").
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.95*cos(t_), 0.50*sin(t_), [0 0 0]);' ...
        'patch(0.84*cos(t_), 0.41*sin(t_), [0.94 0.76 0.74]);' ...
        'color(''black'');' ...
        'plot([-0.45 -0.55], [-0.24 0.24]);' ...
        'plot([-0.18 -0.28], [-0.30 0.30]);' ...
        'plot([0.09 -0.01], [-0.30 0.30]);' ...
        'plot([0.36 0.26], [-0.30 0.30]);' ...
        'plot([0.55 0.47], [-0.22 0.22]);' ...
        ['disp(''' lbl ''');'] ...
        }, newline);
end

function s = synapseIconCode()
    % Small ONE-IN/ONE-OUT connection icon (2026-09-22, per Ben): a compact
    % pass-through axon bar terminating in the postsynaptic marker at the
    % OUTPUT (right) edge — no Vpost sense wire anymore:
    %   Esyn <  0 -> inhibitory  -> SOLID BLACK CIRCLE terminal
    %   Esyn >= 0 -> excitatory  -> WHITE TRIANGLE terminal, heavy black
    %                               edge, tip pointing back toward the
    %                               presynaptic side (Ben 2026-09-09)
    %   unknown Esyn expression -> gray backdrop + "E?"
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
        '    patch([-1 0.70 0.70 -1], [0.34 0.34 0.48 0.48], [0 0 0]);' ...
        '    t_ = linspace(0, 2*pi, 73);' ...
        '    patch(0.20*cos(t_) + 0.76, 0.41*sin(t_), [0 0 0]);' ...
        'else' ...
        '    patch([-1 1 1 -1], [-1 -1 1 1], [0.85 0.94 0.86]);' ...
        '    patch([-1 0.42 0.42 -1], [0.34 0.34 0.48 0.48], [0 0 0]);' ...
        '    patch([0.42 0.98 0.98], [0.41 0.72 -0.30], [0 0 0]);' ...
        '    patch([0.55 0.86 0.86], [0.41 0.62 -0.18], [1 1 1]);' ...
        'end' ...
        }, newline);
end

function s = synsumIconCode()
    % Compact junction glyph: two synaptic lines merging into one.
    s = strjoin({ ...
        'patch([-1 1 1 -1], [-1 -1 1 1], [0.93 0.93 0.93]);' ...
        'color(''black'');' ...
        'plot([-1 -0.1], [0.45 0.10]);' ...
        'plot([-1 -0.1], [0.10 0.10]);' ...
        'plot([-1 -0.1], [-0.45 0.10]);' ...
        'plot([-0.1 1], [0.10 0.10]);' ...
        'plot([-0.1 -0.1], [0.16 0.04]);' ...
        }, newline);
end
