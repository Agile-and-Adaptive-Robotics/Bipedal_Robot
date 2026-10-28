%% sns_build_demo.m — build KneeReflexDemo.slx (SNS circuit + 1-DOF knee model)
%
% Circuit topology (all E/I connections are NonSpikingSynapse blocks; sign from Esyn).
% 2026-09-22 architecture: synapses are ONE-INPUT (Vpre) -> ONE-OUTPUT ([g; g*Esyn])
% and connect to the syn1..syn6 ports of the neuron they synapse onto — synaptic
% summation happens INSIDE the neuron, so there are no Sum blocks in the circuit.
%
%   Descending drive  ------------------------------> MN_ext, MN_flex (Iapp port)
%   Ia spindle (ext) -> SN_Ia_ext --Exc(E=0) ----->  MN_ext syn1   (stretch reflex)
%   Ib GTO (ext)     -> SN_Ib_ext --Inh(E=-72) --->  MN_ext syn2   (autogenic inhibition)
%   Ia spindle (flex)-> SN_Ia_flex -Inh(E=-72) --->  MN_ext syn3   (reciprocal inhibition)
%   Ia spindle (flex)-> SN_Ia_flex -Exc(E=0) ----->  MN_flex syn1
%   Ib GTO (flex)    -> SN_Ib_flex -Inh(E=-72) -->  MN_flex syn2
%   Ia spindle (ext) -> SN_Ia_ext --Inh(E=-72) --->  MN_flex syn3   (reciprocal inhibition)
%
% Sensory neurons (SN_*) convert afferent currents into graded presynaptic voltage.
% Synapse saturation: ThrPre=-45 mV (above rest -52), SlopePre=0.5/mV -> graded, off at rest.
%
% APPEARANCE (Ben 2026-09-09 + 2026-09-22):
%   neurons = circles with GRADED-POTENTIAL waveform; afferents = spindle capsules
%   labeled Ia/Ib; muscles = striated fusiforms; activation = pink pentagon;
%   synapse = small pass-through axon with E-triangle / I-ball terminal, placed
%   CLOSE to the neuron it synapses onto. Tints Okabe-Ito (colorblind-safe).
%
% Knee model (reduced-order 1-DOF, stand-in for the Simscape Multibody import):
%   I*thdd = Tflex - Text + Tload - b*thd - K*(th - th0)
%   theta = 0 deg full extension, +90 deg full flexion.
% Muscle strain over ROM: eps_ext = 0.15*th/ROM, eps_flex = 0.15*(1 - th/ROM).

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder
mdl = 'KneeReflexDemo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '5', 'ScreenColor', 'white', ...
    'UnconnectedInputMsg', 'none');   % unused neuron syn ports auto-ground

%% ---- parameters stored in model PreLoadFcn (so the model opens runnable) ----
paramCmd = strjoin({ ...
    'I_knee = 0.06;      % kg*m^2 below-knee inertia (rig estimate)' ...
    'b_knee = 0.5;        % N*m*s/rad joint damping' ...
    'K_knee = 0.5;       % N*m/rad return spring' ...
    'th0    = 0.26;      % rad spring rest angle (~15 deg)' ...
    'ROM    = 1.5708;    % rad (90 deg) normalizing range' ...
    'r_arm  = 0.035;     % m BPA moment arm' ...
    'Tload  = 0.5;       % N*m constant flexion load torque' ...
    'desc_ext = 4.2;     % nA descending drive to extensor MN' ...
    'desc_flex = 2.5;    % nA descending drive to flexor MN' ...
    'Fmax_ext = 500;     % N' ...
    'Fmax_flex = 450;    % N' ...
    'epsScale = 0.15;    % muscle strain span over full joint ROM' ...
    }, newline);
set_param(mdl, 'PreLoadFcn', paramCmd);
eval(paramCmd);

%% ---- knee model subsystem: [T_flex, T_ext] -> [theta, thd, normalized signals]
knee = [mdl '/KneeModel'];
add_block('simulink/Ports & Subsystems/Subsystem', knee, 'Position', [250 620 430 760]);
delete_line(knee, 'In1/1', 'Out1/1');
delete_block([knee '/In1']);
delete_block([knee '/Out1']);
KB = @(n, p, v) add_block(p, [knee '/' n], 'Position', v);
KB('T_flx',   'simulink/Sources/In1',                  [25 63 55 77]);
KB('T_ext',   'simulink/Sources/In1',                  [25 108 55 122]);
KB('muscleT', 'simulink/Math Operations/Sum',          [95 78 125 112]);   % Tflex - Text
KB('Tload_c', 'simulink/Sources/Constant',             [95 130 125 160]);  % Tload
KB('netT',    'simulink/Math Operations/Sum',          [180 95 210 145]);  % muscle + Tload - damp - spring
KB('dampS',   'simulink/Math Operations/Sum',          [180 190 210 230]); % b*thd + K*(th-th0)
KB('invI',    'simulink/Math Operations/Gain',         [250 105 280 135]);
KB('thd_i',   'simulink/Continuous/Integrator',        [310 105 340 135]);
KB('th_i',    'simulink/Continuous/Integrator',        [370 105 400 135]);
KB('dampG',   'simulink/Math Operations/Gain',         [250 195 280 225]);
KB('sprdfl',  'simulink/Math Operations/Sum',          [180 255 210 285]); % th - th0
KB('th0_c',   'simulink/Sources/Constant',             [95 262 125 292]);
KB('springG', 'simulink/Math Operations/Gain',         [250 255 280 285]);
KB('one_c',   'simulink/Sources/Constant',             [250 340 280 370]);
KB('thn_g',   'simulink/Math Operations/Gain',         [310 335 340 365]); % 1/ROM
KB('thdn_g',  'simulink/Math Operations/Gain',         [310 375 340 405]);
KB('flxstr',  'simulink/Math Operations/Sum',          [370 340 400 370]); % 1 - th/ROM
KB('velflp',  'simulink/Math Operations/Gain',         [370 380 400 410]); % -thd/ROM
KB('strx_g',  'simulink/Math Operations/Gain',         [450 340 480 370]); % epsScale*th/ROM
KB('strf_g',  'simulink/Math Operations/Gain',         [450 385 480 415]); % epsScale*(1-th/ROM)
set_param([knee '/muscleT'], 'Inputs', '+-');
set_param([knee '/netT'],    'Inputs', '++-');
set_param([knee '/dampS'],   'Inputs', '++');
set_param([knee '/sprdfl'],  'Inputs', '+-');
set_param([knee '/flxstr'],  'Inputs', '+-');
set_param([knee '/invI'],    'Gain', '1/I_knee');
set_param([knee '/dampG'],   'Gain', 'b_knee');
set_param([knee '/springG'], 'Gain', 'K_knee');
set_param([knee '/thn_g'],   'Gain', '1/ROM');
set_param([knee '/thdn_g'],  'Gain', '1/ROM');
set_param([knee '/velflp'],  'Gain', '-1/ROM');
set_param([knee '/strx_g'],  'Gain', 'epsScale');
set_param([knee '/strf_g'],  'Gain', 'epsScale');
set_param([knee '/th0_c'],   'Value', 'th0');
set_param([knee '/Tload_c'], 'Value', 'Tload');
set_param([knee '/one_c'],   'Value', '1');
set_param([knee '/th_i'],    'InitialCondition', '0.26');
outs = {'O_th',1,'theta (rad)'; 'O_thd',2,'thd (rad/s)'; 'O_thn',3,'th/ROM'; ...
        'O_thdn',4,'thd/ROM'; 'O_flxstr',5,'flexor stretch (norm)'; ...
        'O_flxvel',6,'flexor velocity (norm)'; 'O_strx',7,'extensor strain'; ...
        'O_strf',8,'flexor strain'};
for k = 1:size(outs, 1)
    add_block('simulink/Sinks/Out1', [knee '/' outs{k,1}], 'Port', num2str(outs{k,2}), ...
        'Position', [540 33+46*(k-1) 570 47+46*(k-1)]);
end
kl = @(a, b) add_line(knee, a, b, 'autorouting', 'on');
kl('T_flx/1', 'muscleT/1');
kl('T_ext/1', 'muscleT/2');
kl('muscleT/1', 'netT/1');
kl('Tload_c/1', 'netT/2');
kl('dampS/1', 'netT/3');   % enters with - sign
kl('netT/1', 'invI/1');
kl('invI/1', 'thd_i/1');
kl('thd_i/1', 'th_i/1');
kl('th_i/1', 'sprdfl/1');
kl('th0_c/1', 'sprdfl/2');
kl('sprdfl/1', 'springG/1');
kl('thd_i/1', 'dampG/1');
kl('dampG/1', 'dampS/1');
kl('springG/1', 'dampS/2');
kl('th_i/1', 'O_th/1');
kl('thd_i/1', 'O_thd/1');
kl('th_i/1', 'thn_g/1');
kl('thd_i/1', 'thdn_g/1');
kl('thn_g/1', 'O_thn/1');
kl('thdn_g/1', 'O_thdn/1');
kl('one_c/1', 'flxstr/1');
kl('thn_g/1', 'flxstr/2');
kl('flxstr/1', 'O_flxstr/1');
kl('thdn_g/1', 'velflp/1');
kl('velflp/1', 'O_flxvel/1');
kl('thn_g/1', 'strx_g/1');
kl('flxstr/1', 'strf_g/1');
kl('strx_g/1', 'O_strx/1');
kl('strf_g/1', 'O_strf/1');
m = Simulink.Mask.create(knee);
m.Type = 'SNS Knee Model (1-DOF)';
m.Description = ['Reduced-order 1-DOF knee: I*thdd = Tflex - Text + Tload - b*thd - ' ...
    'K*(th-th0). theta = 0 deg full extension, +90 deg full flexion. Outputs theta, ' ...
    'thd, and the normalized stretch/velocity/strain signals for the afferents.'];
m.Display = kneeIconCode();
set_param(knee, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
    'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');

%% ---- afferents (spindle capsules / GTO capsules) ----
snsInst(mdl, 'SNS_Library/IaMuscleSpindle', 'Ia_ext',  [40 40 100 100], 'Imax', '10', 'Wl', '6', 'Wv', '8');
snsInst(mdl, 'SNS_Library/IaMuscleSpindle', 'Ia_flex', [40 190 100 250], 'Imax', '10', 'Wl', '6', 'Wv', '8');
snsInst(mdl, 'SNS_Library/IbGolgiTendon', 'Ib_ext',    [40 340 100 400], 'Imax', '10', 'Kf', '0.02');
snsInst(mdl, 'SNS_Library/IbGolgiTendon', 'Ib_flex',   [40 440 100 500], 'Imax', '10', 'Kf', '0.0222');
add_line(mdl, 'KneeModel/3', 'Ia_ext/1', 'autorouting', 'on');
add_line(mdl, 'KneeModel/4', 'Ia_ext/2', 'autorouting', 'on');
add_line(mdl, 'KneeModel/5', 'Ia_flex/1', 'autorouting', 'on');
add_line(mdl, 'KneeModel/6', 'Ia_flex/2', 'autorouting', 'on');

%% ---- sensory neurons: afferent current -> graded presynaptic voltage ----
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ia_ext',  [180 40 260 120], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ia_flex', [180 190 260 270], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ib_ext',  [180 340 260 420], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ib_flex', [180 440 260 520], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
add_line(mdl, 'Ia_ext/1', 'SN_Ia_ext/1', 'autorouting', 'on');
add_line(mdl, 'Ia_flex/1', 'SN_Ia_flex/1', 'autorouting', 'on');
add_line(mdl, 'Ib_ext/1', 'SN_Ib_ext/1', 'autorouting', 'on');
add_line(mdl, 'Ib_flex/1', 'SN_Ib_flex/1', 'autorouting', 'on');

%% ---- motoneurons (graded-waveform non-spiking neurons) ----
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'MN_ext',  [660 40 760 200], 'Vrest', '-52', 'Gm', '0.5', 'Cm', '2.5', 'Thr', '-45', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'MN_flex', [660 320 760 480], 'Vrest', '-52', 'Gm', '0.5', 'Cm', '2.5', 'Thr', '-45', 'Slope', '1');
add_block('simulink/Sources/Constant', [mdl '/desc_ext_c'], 'Value', 'desc_ext', 'Position', [560 -10 600 14]);
add_block('simulink/Sources/Constant', [mdl '/desc_flex_c'], 'Value', 'desc_flex', 'Position', [560 270 600 294]);
add_line(mdl, 'desc_ext_c/1', 'MN_ext/1', 'autorouting', 'on');
add_line(mdl, 'desc_flex_c/1', 'MN_flex/1', 'autorouting', 'on');

%% ---- synapses: small, ONE input (Vpre), placed against their postsynaptic MN ----
% MN_ext syn ports (block ports 2,3,4 = syn1..syn3)
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaext_exc',         [580 58 620 90],  'gmax', '0.18', 'Esyn', '0',   'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Ibext_inh',         [580 100 620 132], 'gmax', '0.20', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaflex_inh_on_ext', [580 142 620 174], 'gmax', '0.15', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% MN_flex syn ports
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaflex_exc',        [580 338 620 370], 'gmax', '0.18', 'Esyn', '0',   'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Ibflex_inh',        [580 380 620 412], 'gmax', '0.20', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaext_inh_on_flex', [580 422 620 454], 'gmax', '0.15', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% presynaptic: sensory neuron V -> synapse Vpre
add_line(mdl, 'SN_Ia_ext/1', 'syn_Iaext_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib_ext/1', 'syn_Ibext_inh/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_flex/1', 'syn_Iaflex_inh_on_ext/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_flex/1', 'syn_Iaflex_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib_flex/1', 'syn_Ibflex_inh/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_ext/1', 'syn_Iaext_inh_on_flex/1', 'autorouting', 'on');
% synapse output -> postsynaptic neuron syn port (next free syn1..syn3)
add_line(mdl, 'syn_Iaext_exc/1', 'MN_ext/2', 'autorouting', 'on');
add_line(mdl, 'syn_Ibext_inh/1', 'MN_ext/3', 'autorouting', 'on');
add_line(mdl, 'syn_Iaflex_inh_on_ext/1', 'MN_ext/4', 'autorouting', 'on');
add_line(mdl, 'syn_Iaflex_exc/1', 'MN_flex/2', 'autorouting', 'on');
add_line(mdl, 'syn_Ibflex_inh/1', 'MN_flex/3', 'autorouting', 'on');
add_line(mdl, 'syn_Iaext_inh_on_flex/1', 'MN_flex/4', 'autorouting', 'on');

%% ---- muscles: MN S -> pentagon activation -> striated BPA -> torque ----
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_ext',  [820 100 880 160], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_flex', [820 380 880 440], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_ext',  [930 100 1000 160], 'Fmax', 'Fmax_ext', 'epsMax', '0.25');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_flex', [930 380 1000 440], 'Fmax', 'Fmax_flex', 'epsMax', '0.25');
add_block('simulink/Math Operations/Gain', [mdl '/r_ext'], 'Gain', 'r_arm', 'Position', [1040 115 1070 145]);
add_block('simulink/Math Operations/Gain', [mdl '/r_flex'], 'Gain', 'r_arm', 'Position', [1040 395 1070 425]);
add_line(mdl, 'MN_ext/2', 'Act_ext/1', 'autorouting', 'on');
add_line(mdl, 'MN_flex/2', 'Act_flex/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'BPA_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'BPA_flex/1', 'autorouting', 'on');
add_line(mdl, 'KneeModel/7', 'BPA_ext/2', 'autorouting', 'on');
add_line(mdl, 'KneeModel/8', 'BPA_flex/2', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'Ib_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'Ib_flex/1', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'r_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'r_flex/1', 'autorouting', 'on');
add_line(mdl, 'r_flex/1', 'KneeModel/1', 'autorouting', 'on');
add_line(mdl, 'r_ext/1', 'KneeModel/2', 'autorouting', 'on');

%% ---- logging ----
logNames = {'th', 'thd', 'A_ext', 'A_flex', 'V_MN_ext', 'V_MN_flex', 'F_ext', 'F_flex', 'Ia_ext_c', 'Ib_ext_c'};
for k = 1:numel(logNames)
    y = 60 + 55*(k-1);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logNames{k}], ...
        'VariableName', ['log_' logNames{k}], 'SaveFormat', 'Timeseries', ...
        'Position', [1160 y 1230 y+30]);
end
add_line(mdl, 'KneeModel/1', 'log_th/1', 'autorouting', 'on');
add_line(mdl, 'KneeModel/2', 'log_thd/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'log_A_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'log_A_flex/1', 'autorouting', 'on');
add_line(mdl, 'MN_ext/1', 'log_V_MN_ext/1', 'autorouting', 'on');
add_line(mdl, 'MN_flex/1', 'log_V_MN_flex/1', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'log_F_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'log_F_flex/1', 'autorouting', 'on');
add_line(mdl, 'Ia_ext/1', 'log_Ia_ext_c/1', 'autorouting', 'on');
add_line(mdl, 'Ib_ext/1', 'log_Ib_ext_c/1', 'autorouting', 'on');

%% ---- annotations: title + diagram-language legend ----
try
    anno = Simulink.Annotation(mdl, ...
        'SNS knee reflex demo: non-spiking RC neurons + E/I synapses + antagonist BPAs on a 1-DOF knee model');
    anno.Position = [40 -75 900 -35];
catch
end
try
    legend_txt = ['Diagram language: circle + graded waveform = non-spiking neuron; spindle = Ia afferent; ' ...
        'capsule "Ib" = Golgi tendon; pentagon = activation; striated fusiform = muscle/BPA.  ' ...
        'Connection markers: solid black circle = INHIBITORY (Esyn < 0), ' ...
        'open triangle = EXCITATORY (Esyn >= 0). Tints: Okabe-Ito (colorblind-safe).'];
    anno2 = Simulink.Annotation(mdl, legend_txt);
    anno2.Position = [40 -25 1150 5];
    anno2.FontSize = 10;
catch
end
% synapse names off (a synapse is a connection, not a labelled component)
for nm = {'syn_Iaext_exc', 'syn_Ibext_inh', 'syn_Iaflex_inh_on_ext', ...
          'syn_Iaflex_exc', 'syn_Ibflex_inh', 'syn_Iaext_inh_on_flex'}
    try
        set_param([mdl '/' nm{1}], 'ShowName', 'off', 'NamePlacement', 'alternate');
    catch
    end
end

save_system(mdl);
fprintf('KneeReflexDemo.slx built (2026-09-22 architecture).\n');

function h = snsInst(mdl, libPath, name, pos, varargin)
    h = add_block(libPath, [mdl '/' name], 'Position', pos, varargin{:});
end

function s = kneeIconCode()
    % Femur + shank glyph with the joint circle.
    s = strjoin({ ...
        'patch([-1 1 1 -1], [-1 -1 1 1], [0 0 0]);' ...
        'patch([-0.92 0.92 0.92 -0.92], [-0.84 -0.84 0.84 0.84], [0.93 0.93 0.93]);' ...
        'color(''black'');' ...
        'plot([-0.55 0.05], [0.75 -0.05]);' ...
        'plot([0.05 0.60], [-0.05 -0.75]);' ...
        't_ = linspace(0, 2*pi, 25);' ...
        'patch(0.10*cos(t_) + 0.05, 0.10*sin(t_) - 0.05, [0 0 0]);' ...
        'patch(0.055*cos(t_) + 0.05, 0.055*sin(t_) - 0.05, [1 1 1]);' ...
        }, newline);
end
