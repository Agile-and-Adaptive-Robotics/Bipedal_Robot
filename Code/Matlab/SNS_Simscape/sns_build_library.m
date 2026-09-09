%% sns_build_library.m — build SNS_Library.slx
% Synthetic Nervous System block library for the Bipedal_Robot Simscape/Simulink work.
%
% Neuron model = non-spiking leaky integrator, literally an RC membrane:
%   Cm dV/dt = Gm*(Vrest - V) + sum(Isyn)
%   tau_m = Cm/Gm  (nF/uS = ms). Same convention as Animatlab non-spiking neurons
%   and sns_toolbox NonSpikingNeuron (Szczecinski et al. 2017).
% Synapse model (Animatlab + SNS toolbox non-spiking chemical synapse):
%   Isyn = gmax * Sat(Vpre) * (Esyn - Vpost)   [uS*mV = nA]
%   Sat(Vpre) = clip((Vpre - ThrPre)/SlopePre, 0, 1)
%   Esyn > Vrest  -> excitatory (depolarizing), Esyn < Vrest -> inhibitory.
%
% Parameter conventions copied from SNS toolbox defaults / Animatlab:
%   Vrest = -52 mV, tau_m = 50 ms, synapse Thr = -55 mV, slope = 1 /mV,
%   excitatory Esyn = 0 mV, inhibitory Esyn = -72 mV.
%
% --- APPEARANCE CONVENTIONS (Ben, 2026-09-09) ---------------------------------
% Diagram language follows Rybak/Shevtsova, Animatlab, and Szczecinski et al.
% 2017 "functional subnetwork" papers:
%   * neurons                = open circles (black edge, white fill)
%   * sensory afferents      = open circles labeled Ia / Ib
%   * muscles                = fusiform ellipses (light green tint)
%   * EXCITATORY connection  = WHITE TRIANGLE with black edges (INVERTED per
%     Ben 2026-09-09: tip points back toward the presynaptic side, flat base
%     at the postsynaptic side)
%   * INHIBITORY connection  = SOLID BLACK CIRCLE
% The E/I marker on NonSpikingSynapse is chosen AUTOMATICALLY from the sign of
% Esyn, so the icon always tells the truth about the connection. Shape is the
% primary code; block tints use the colorblind-safe Okabe-Ito palette as a
% redundant cue (excitatory = light orange, inhibitory = light blue; orange/blue
% reads correctly under deuteranopia/protanopia/tritanopia).
% Keep every masked block SQUARE (width == height) — mask icons autoscale to the
% block rectangle, so a non-square block turns circles into ellipses.

lib = 'SNS_Library';
if bdIsLoaded(lib), close_system(lib, 0); end
if exist([lib '.slx'], 'file'), delete([lib '.slx']); end
new_system(lib, 'Library');
load_system(lib);

%% ---------------- NonSpikingNeuron ----------------
blk = [lib '/NonSpikingNeuron'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 40 160 160]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Isyn'], 'Port', '1', 'Position', [25 78 55 92]);
add_block('simulink/Sources/Constant', [blk '/Vrest_c'], 'Value', 'Vrest', 'Position', [90 150 120 180]);
add_block('simulink/Math Operations/Sum', [blk '/dV'], 'Inputs', '++', 'Position', [250 70 280 100]);
add_block('simulink/Math Operations/Sum', [blk '/dVm'], 'Inputs', '-+', 'Position', [160 108 190 138]);
add_block('simulink/Math Operations/Gain', [blk '/Gm'], 'Gain', 'Gm', 'Position', [210 113 240 143]);
add_block('simulink/Math Operations/Gain', [blk '/membrane'], 'Gain', '1000/Cm', 'Position', [310 70 340 100]);
add_block('simulink/Continuous/Integrator', [blk '/membraneRC'], 'InitialCondition', 'Vrest', 'Position', [370 68 400 102]);
add_block('simulink/Math Operations/Sum', [blk '/thrSub'], 'Inputs', '+-', 'Position', [400 130 430 160]);
add_block('simulink/Sources/Constant', [blk '/Thr_c'], 'Value', 'Thr', 'Position', [340 175 370 205]);
add_block('simulink/Math Operations/Gain', [blk '/invSlope'], 'Gain', '1/Slope', 'Position', [460 130 490 160]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', '1', 'LowerLimit', '0', 'Position', [520 130 550 160]);
add_block('simulink/Sinks/Out1', [blk '/V_mV'], 'Port', '1', 'Position', [590 78 620 92]);
add_block('simulink/Sinks/Out1', [blk '/S_drive'], 'Port', '2', 'Position', [590 138 620 152]);
add_line(blk, 'Isyn/1', 'dV/2', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'thrSub/1', 'autorouting', 'on');
add_line(blk, 'Thr_c/1', 'thrSub/2', 'autorouting', 'on');
add_line(blk, 'thrSub/1', 'invSlope/1', 'autorouting', 'on');
add_line(blk, 'Vrest_c/1', 'dVm/2', 'autorouting', 'on');
add_line(blk, 'dVm/1', 'Gm/1', 'autorouting', 'on');
add_line(blk, 'Gm/1', 'dV/1', 'autorouting', 'on');
add_line(blk, 'dV/1', 'membrane/1', 'autorouting', 'on');
add_line(blk, 'membrane/1', 'membraneRC/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'V_mV/1', 'autorouting', 'on');
add_line(blk, 'membraneRC/1', 'dVm/1', 'autorouting', 'on');
add_line(blk, 'invSlope/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'S_drive/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS NonSpiking Neuron';
m.Description = ['Non-spiking leaky integrate-and-fire neuron (RC membrane): ' ...
    'Cm*dV/dt = Gm*(Vrest-V) + Isyn. Outputs membrane potential V [mV] and normalized drive S(V) [0..1]. ' ...
    'Parameter conventions: Animatlab / sns_toolbox (Szczecinski et al. 2017).'];
m.Display = neuronIconCode('NS');
finishMask(m, blk, {'Vrest','Resting potential Vrest (mV)','-52'; ...
                    'Gm','Membrane conductance Gm (uS)','0.1'; ...
                    'Cm','Membrane capacitance Cm (nF)','5'; ...
                    'Thr','Saturation threshold Thr (mV)','-55'; ...
                    'Slope','Saturation slope (1/mV)','1'});

%% ---------------- NonSpikingSynapse ----------------
blk = [lib '/NonSpikingSynapse'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 220 140 320]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Vpre'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Sources/In1', [blk '/Vpost'], 'Port', '2', 'Position', [25 178 55 192]);
add_block('simulink/Math Operations/Sum', [blk '/thrSub'], 'Inputs', '+-', 'Position', [60 55 90 85]);
add_block('simulink/Sources/Constant', [blk '/ThrPre_c'], 'Value', 'ThrPre', 'Position', [20 110 50 140]);
add_block('simulink/Math Operations/Gain', [blk '/invSlope'], 'Gain', '1/SlopePre', 'Position', [110 55 140 85]);
add_block('simulink/Discontinuities/Saturation', [blk '/Sat'], 'UpperLimit', '1', 'LowerLimit', '0', 'Position', [150 55 180 85]);
add_block('simulink/Math Operations/Gain', [blk '/gmax_g'], 'Gain', 'gmax', 'Position', [210 55 240 85]);
add_block('simulink/Math Operations/Sum', [blk '/drive'], 'Inputs', '+-', 'Position', [150 170 180 200]);
add_block('simulink/Sources/Constant', [blk '/Esyn_c'], 'Value', 'Esyn', 'Position', [90 195 120 225]);
add_block('simulink/Math Operations/Product', [blk '/Iout'], 'Inputs', '2', 'Position', [290 100 320 130]);
add_block('simulink/Sinks/Out1', [blk '/Isyn_nA'], 'Position', [360 108 390 122]);
add_line(blk, 'Vpre/1', 'thrSub/1', 'autorouting', 'on');
add_line(blk, 'ThrPre_c/1', 'thrSub/2', 'autorouting', 'on');
add_line(blk, 'thrSub/1', 'invSlope/1', 'autorouting', 'on');
add_line(blk, 'invSlope/1', 'Sat/1', 'autorouting', 'on');
add_line(blk, 'Sat/1', 'gmax_g/1', 'autorouting', 'on');
add_line(blk, 'gmax_g/1', 'Iout/1', 'autorouting', 'on');
add_line(blk, 'Esyn_c/1', 'drive/1', 'autorouting', 'on');
add_line(blk, 'Vpost/1', 'drive/2', 'autorouting', 'on');
add_line(blk, 'drive/1', 'Iout/2', 'autorouting', 'on');
add_line(blk, 'Iout/1', 'Isyn_nA/1', 'autorouting', 'on');
m = Simulink.Mask.create(blk);
m.Type = 'SNS NonSpiking Synapse';
m.Description = ['Non-spiking chemical synapse: Isyn = gmax*Sat(Vpre)*(Esyn - Vpost) [nA]. ' ...
    'Esyn >= 0 mV -> EXCITATORY (icon: white triangle, black edges); ' ...
    'Esyn < 0 mV -> INHIBITORY (icon: solid black circle). ' ...
    'ThrPre/SlopePre presynaptic. Diagram language: Szczecinski et al. 2017 Fig. 2.'];
m.Display = synapseIconCode();
finishMask(m, blk, {'gmax','Max conductance gmax (uS)','1'; ...
                    'Esyn','Reversal potential Esyn (mV). >=0 excit, <0 inhib','0'; ...
                    'ThrPre','Pre saturation threshold ThrPre (mV)','-55'; ...
                    'SlopePre','Pre saturation slope SlopePre (1/mV)','1'});

%% ---------------- SpikingLIFNeuron ----------------
blk = [lib '/SpikingLIFNeuron'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 380 160 500]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/Isyn'], 'Port', '1', 'Position', [25 78 55 92]);
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
add_line(blk, 'Isyn/1', 'sumI/2', 'autorouting', 'on');
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
m.Description = ['Spiking leaky integrate-and-fire neuron: tau*dV/dt = (Vrest-V) + Rm*Isyn; ' ...
    'V>=Vth -> output spike flag and membrane resets to Vrest.'];
m.Display = lifIconCode();
finishMask(m, blk, {'Vrest','Resting potential Vrest (mV)','-60'; ...
                    'Vth','Spike threshold Vth (mV)','-50'; ...
                    'Rm','Membrane resistance Rm (MOhm)','10'; ...
                    'tau','Membrane time constant tau (ms)','20'});

%% ---------------- IaMuscleSpindle ----------------
blk = [lib '/IaMuscleSpindle'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 560 160 680]);
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
m.Description = 'Ia spindle afferent: Isyn = clip(Wl*stretch + Wv*velocity, 0, Imax) [nA]. Inputs are normalized (0..1) length and velocity signals.';
m.Display = afferentIconCode('Ia');
finishMask(m, blk, {'Imax','Peak current Imax (nA)','10'; ...
                    'Wl','Length weight Wl (nA)','6'; ...
                    'Wv','Velocity weight Wv (nA)','8'});

%% ---------------- IbGolgiTendon ----------------
blk = [lib '/IbGolgiTendon'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 740 160 860]);
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
m.Description = 'Ib Golgi tendon organ afferent: Isyn = clip(Kf*forceNorm, 0, Imax) [nA]. Input is normalized force (0..1).';
m.Display = afferentIconCode('Ib');
finishMask(m, blk, {'Imax','Peak current Imax (nA)','10'; ...
                    'Kf','Force weight Kf (nA)','10'});

%% ---------------- MuscleActivation ----------------
blk = [lib '/MuscleActivation'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 920 160 1040]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
add_block('simulink/Sources/In1', [blk '/S'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Math Operations/Sum', [blk '/dA'], 'Inputs', '-+', 'Position', [90 55 120 85]);
add_block('simulink/Math Operations/Gain', [blk '/tauInv'], 'Gain', '1000/tauAct', 'Position', [150 55 180 85]);
add_block('simulink/Continuous/Integrator', [blk '/Aint'], 'InitialCondition', '0', 'Position', [210 55 240 85]);
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
m.Description = 'First-order activation dynamics: tauAct*dA/dt = (S - A). Maps motoneuron drive S(V) [0..1] to muscle activation A [0..1].';
m.Display = muscleIconCode('MA');
finishMask(m, blk, {'tauAct','Activation time constant tauAct (ms)','50'});

%% ---------------- BPAForce ----------------
blk = [lib '/BPAForce'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 1100 160 1220]);
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
m.Description = 'Simple BPA (pneumatic muscle) force: F = Fmax*A*max(0, epsMax - strain). Replace with MonoPam/Xi-corrected prediction when available.';
m.Display = muscleIconCode('BPA');
finishMask(m, blk, {'Fmax','Max isometric force Fmax (N)','500'; ...
                    'epsMax','Max contraction strain epsMax (0..1)','0.25'});

save_system(lib);
fprintf('SNS_Library.slx saved with 7 library blocks (journal diagram icons).\n');

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

function s = neuronIconCode(lbl)
    % Open circle, black edge, white fill, small label inside.
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.92*cos(t_), 0.92*sin(t_), [1 1 1]);' ...
        'color(''black'');' ...
        'plot(0.92*cos(t_), 0.92*sin(t_));' ...
        ['disp(''' lbl ''');'] ...
        }, newline);
end

function s = lifIconCode()
    % Open circle with a spike glyph inside (distinguishes spiking cell).
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.92*cos(t_), 0.92*sin(t_), [1 1 1]);' ...
        'color(''black'');' ...
        'plot(0.92*cos(t_), 0.92*sin(t_));' ...
        'plot([-0.55 -0.25 0.0 0.25 0.55], [-0.15 -0.15 0.45 -0.15 -0.15]);' ...
        }, newline);
end

function s = afferentIconCode(lbl)
    % Open circle, light-gray tint (sensory), afferent label inside.
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.92*cos(t_), 0.92*sin(t_), [0.93 0.93 0.93]);' ...
        'color(''black'');' ...
        'plot(0.92*cos(t_), 0.92*sin(t_));' ...
        ['disp(''' lbl ''');'] ...
        }, newline);
end

function s = muscleIconCode(lbl)
    % Fusiform (spindle-shaped) muscle, light green tint (Okabe-Ito green).
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.92*cos(t_), 0.42*sin(t_), [0.82 0.92 0.87]);' ...
        'color(''black'');' ...
        'plot(0.92*cos(t_), 0.42*sin(t_));' ...
        ['disp(''' lbl ''');'] ...
        }, newline);
end

function s = synapseIconCode()
    % E/I marker chosen automatically from the SIGN of Esyn:
    %   Esyn <  0 -> inhibitory  -> SOLID BLACK CIRCLE on light-blue backdrop
    %   Esyn >= 0 -> excitatory  -> WHITE TRIANGLE, black edges, on light-orange
    %   unknown expression   -> gray backdrop + "E?" text
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
        '    patch([-1 1 1 -1], [-1 -1 1 1], [0.87 0.93 0.97]);' ...
        '    t_ = linspace(0, 2*pi, 73);' ...
        '    patch(0.78*cos(t_), 0.78*sin(t_), [0 0 0]);' ...
        'else' ...
        '    patch([-1 1 1 -1], [-1 -1 1 1], [1 0.94 0.85]);' ...
        '    patch([0 -0.85 0.85], [-0.8 0.62 0.62], [1 1 1]);' ...
        '    color(''black'');' ...
        '    plot([0 -0.85 0.85 0], [-0.8 0.62 0.62 0]);' ...
        'end' ...
        }, newline);
end
