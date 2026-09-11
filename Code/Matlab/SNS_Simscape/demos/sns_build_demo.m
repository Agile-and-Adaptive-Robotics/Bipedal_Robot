%% sns_build_demo.m — build KneeReflexDemo.slx (SNS circuit + reduced-order knee plant)
%
% Circuit topology (all E/I connections are NonSpikingSynapse blocks; sign from Esyn):
%   Descending drive  ------------------------------> MN_ext, MN_flex (bias below threshold)
%   Ia spindle (ext) -> SN_Ia_ext --Exc(E=0) -----> (+) MN_ext     (stretch reflex)
%   Ib GTO (ext)     -> SN_Ib_ext --Inh(E=-72) ---> (-) MN_ext     (autogenic inhibition)
%   Ia spindle (flex)-> SN_Ia_flex -Inh(E=-72) ---> (-) MN_ext     (reciprocal inhibition)
%   Ia spindle (flex)-> SN_Ia_flex -Exc(E=0) -----> (+) MN_flex
%   Ib GTO (flex)    -> SN_Ib_flex -Inh(E=-72) -->  (-) MN_flex
%   Ia spindle (ext) -> SN_Ia_ext --Inh(E=-72) ---> (-) MN_flex     (reciprocal inhibition)
%
% Sensory neurons (SN_*) convert afferent currents into graded presynaptic voltage.
% Synapse saturation: ThrPre=-45 mV (above rest -52), SlopePre=0.5/mV -> graded, off at rest.
%
% APPEARANCE (Szczecinski 2017 / Rybak / Animatlab diagram language):
%   neurons = circles, afferents labeled Ia/Ib, muscles = ellipses,
%   synapse icon shows SOLID BLACK CIRCLE = inhibitory (Esyn<0),
%   WHITE TRIANGLE w/ black edges = excitatory (Esyn>=0), auto from Esyn sign.
%   Tints are Okabe-Ito (colorblind-safe); shape is the primary code.
%   All neural blocks are SQUARE so mask-icon circles render round.
%
% Plant (reduced-order 1-DOF knee, stand-in for the Simscape Multibody import):
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
set_param(mdl, 'Solver', 'ode45', 'StopTime', '5', 'ScreenColor', 'white');

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

%% ---- plant (bottom center): I*thdd = Tflex - Text + Tload - b thd - K(th-th0) ----
add_block('simulink/Math Operations/Sum', [mdl '/netTorque'], 'Inputs', '+++', 'Position', [640 560 670 620]);
add_block('simulink/Math Operations/Sum', [mdl '/netTorque2'], 'Inputs', '--+', 'Position', [560 565 590 615]);
add_block('simulink/Math Operations/Gain', [mdl '/invI'], 'Gain', '1/I_knee', 'Position', [700 575 730 605]);
add_block('simulink/Continuous/Integrator', [mdl '/thd_int'], 'InitialCondition', '0', 'Position', [760 575 790 605]);
add_block('simulink/Continuous/Integrator', [mdl '/th_int'], 'InitialCondition', '0.26', 'Position', [820 575 850 605]);
add_block('simulink/Math Operations/Gain', [mdl '/damp'], 'Gain', 'b_knee', 'Position', [760 650 790 680]);
add_block('simulink/Math Operations/Gain', [mdl '/spring'], 'Gain', 'K_knee', 'Position', [760 700 790 730]);
add_block('simulink/Math Operations/Sum', [mdl '/spr_defl'], 'Inputs', '+-', 'Position', [700 705 730 735]);
add_block('simulink/Sources/Constant', [mdl '/th0_c'], 'Value', 'th0', 'Position', [640 715 670 745]);
add_block('simulink/Sources/Constant', [mdl '/Tload_c'], 'Value', 'Tload', 'Position', [490 585 520 615]);
% netTorque: [+,+,+] = [netMuscle, -damp-spring, Tload]; netTorque2: [-,-,+]= [Text?, ...]
add_line(mdl, 'netTorque/1', 'invI/1', 'autorouting', 'on');
add_line(mdl, 'invI/1', 'thd_int/1', 'autorouting', 'on');
add_line(mdl, 'thd_int/1', 'th_int/1', 'autorouting', 'on');
add_line(mdl, 'thd_int/1', 'damp/1', 'autorouting', 'on');
add_line(mdl, 'th_int/1', 'spr_defl/1', 'autorouting', 'on');
add_line(mdl, 'th0_c/1', 'spr_defl/2', 'autorouting', 'on');
add_line(mdl, 'spr_defl/1', 'spring/1', 'autorouting', 'on');
% muscle torque net: Tflex - Text
add_block('simulink/Math Operations/Sum', [mdl '/muscleT'], 'Inputs', '+-', 'Position', [400 570 430 600]);
add_line(mdl, 'muscleT/1', 'netTorque/1', 'autorouting', 'on');
add_line(mdl, 'netTorque2/1', 'netTorque/2', 'autorouting', 'on');   % -(damp+spring)
add_line(mdl, 'damp/1', 'netTorque2/1', 'autorouting', 'on');
add_line(mdl, 'spring/1', 'netTorque2/2', 'autorouting', 'on');
add_line(mdl, 'Tload_c/1', 'netTorque2/3', 'autorouting', 'on');

%% ---- kinematic normalizations ----
add_block('simulink/Math Operations/Gain', [mdl '/th_norm'], 'Gain', '1/ROM', 'Position', [900 575 930 605]);
add_block('simulink/Math Operations/Gain', [mdl '/thd_norm'], 'Gain', '1/ROM', 'Position', [900 640 930 670]);
add_block('simulink/Math Operations/Sum', [mdl '/flex_stretch'], 'Inputs', '+-', 'Position', [900 700 930 730]);
add_block('simulink/Sources/Constant', [mdl '/one_c'], 'Value', '1', 'Position', [840 710 870 740]);
add_block('simulink/Math Operations/Gain', [mdl '/vel_flip'], 'Gain', '-1', 'Position', [900 770 930 800]);
add_block('simulink/Math Operations/Gain', [mdl '/strain_ext_g'], 'Gain', 'epsScale', 'Position', [900 830 930 860]);
add_block('simulink/Math Operations/Gain', [mdl '/strain_flex_g'], 'Gain', 'epsScale', 'Position', [900 890 930 920]);
add_line(mdl, 'th_int/1', 'th_norm/1', 'autorouting', 'on');
add_line(mdl, 'thd_int/1', 'thd_norm/1', 'autorouting', 'on');
add_line(mdl, 'one_c/1', 'flex_stretch/1', 'autorouting', 'on');
add_line(mdl, 'th_norm/1', 'flex_stretch/2', 'autorouting', 'on');
add_line(mdl, 'thd_norm/1', 'vel_flip/1', 'autorouting', 'on');
add_line(mdl, 'th_norm/1', 'strain_ext_g/1', 'autorouting', 'on');
add_line(mdl, 'flex_stretch/1', 'strain_flex_g/1', 'autorouting', 'on');

%% ---- afferents (sensory currents, nA) — square blocks so circle icons stay round ----
snsInst(mdl, 'SNS_Library/IaMuscleSpindle', 'Ia_ext', [120 60 200 140], 'Imax', '10', 'Wl', '6', 'Wv', '8');
snsInst(mdl, 'SNS_Library/IaMuscleSpindle', 'Ia_flex', [120 200 200 280], 'Imax', '10', 'Wl', '6', 'Wv', '8');
snsInst(mdl, 'SNS_Library/IbGolgiTendon', 'Ib_ext', [120 360 200 440], 'Imax', '10', 'Kf', '0.02');
snsInst(mdl, 'SNS_Library/IbGolgiTendon', 'Ib_flex', [120 460 200 540], 'Imax', '10', 'Kf', '0.0222');
add_line(mdl, 'th_norm/1', 'Ia_ext/1', 'autorouting', 'on');
add_line(mdl, 'thd_norm/1', 'Ia_ext/2', 'autorouting', 'on');
add_line(mdl, 'flex_stretch/1', 'Ia_flex/1', 'autorouting', 'on');
add_line(mdl, 'vel_flip/1', 'Ia_flex/2', 'autorouting', 'on');

%% ---- sensory neurons: afferent current -> graded presynaptic voltage ----
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ia_ext', [250 40 350 140], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ia_flex', [250 190 350 290], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ib_ext', [250 350 350 450], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'SN_Ib_flex', [250 490 350 590], 'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
add_line(mdl, 'Ia_ext/1', 'SN_Ia_ext/1', 'autorouting', 'on');
add_line(mdl, 'Ia_flex/1', 'SN_Ia_flex/1', 'autorouting', 'on');
add_line(mdl, 'Ib_ext/1', 'SN_Ib_ext/1', 'autorouting', 'on');
add_line(mdl, 'Ib_flex/1', 'SN_Ib_flex/1', 'autorouting', 'on');

%% ---- motoneurons (non-spiking LIF / RC) ----
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'MN_ext', [560 80 670 190], 'Vrest', '-52', 'Gm', '0.5', 'Cm', '2.5', 'Thr', '-45', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'MN_flex', [560 260 670 370], 'Vrest', '-52', 'Gm', '0.5', 'Cm', '2.5', 'Thr', '-45', 'Slope', '1');

%% ---- synapse summing nodes ----
add_block('simulink/Math Operations/Sum', [mdl '/sumMN_ext'], 'Inputs', '++++', 'Position', [480 95 510 155]);
add_block('simulink/Math Operations/Sum', [mdl '/sumMN_flex'], 'Inputs', '++++', 'Position', [480 255 510 315]);
add_block('simulink/Sources/Constant', [mdl '/desc_ext_c'], 'Value', 'desc_ext', 'Position', [470 165 500 195]);
add_block('simulink/Sources/Constant', [mdl '/desc_flex_c'], 'Value', 'desc_flex', 'Position', [470 325 500 355]);
add_line(mdl, 'desc_ext_c/1', 'sumMN_ext/1', 'autorouting', 'on');
add_line(mdl, 'desc_flex_c/1', 'sumMN_flex/1', 'autorouting', 'on');
add_line(mdl, 'sumMN_ext/1', 'MN_ext/1', 'autorouting', 'on');
add_line(mdl, 'sumMN_flex/1', 'MN_flex/1', 'autorouting', 'on');

%% ---- synapses: [pre V, post V] -> current (Esyn sets E vs I; icon auto-draws marker) ----
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaext_exc', [380 60 450 130], 'gmax', '0.18', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Ibext_inh', [380 150 450 220], 'gmax', '0.20', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaflex_exc', [380 240 450 310], 'gmax', '0.18', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaflex_inh_on_ext', [380 330 450 400], 'gmax', '0.15', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Ibflex_inh', [380 420 450 490], 'gmax', '0.20', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'syn_Iaext_inh_on_flex', [380 510 450 580], 'gmax', '0.15', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% afferent neuron V -> synapse pre
add_line(mdl, 'SN_Ia_ext/1', 'syn_Iaext_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib_ext/1', 'syn_Ibext_inh/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_flex/1', 'syn_Iaflex_inh_on_ext/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_flex/1', 'syn_Iaflex_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib_flex/1', 'syn_Ibflex_inh/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_ext/1', 'syn_Iaext_inh_on_flex/1', 'autorouting', 'on');
% MN membrane V -> synapse post
add_line(mdl, 'MN_ext/1', 'syn_Iaext_exc/2', 'autorouting', 'on');
add_line(mdl, 'MN_ext/1', 'syn_Ibext_inh/2', 'autorouting', 'on');
add_line(mdl, 'MN_ext/1', 'syn_Iaflex_inh_on_ext/2', 'autorouting', 'on');
add_line(mdl, 'MN_flex/1', 'syn_Iaflex_exc/2', 'autorouting', 'on');
add_line(mdl, 'MN_flex/1', 'syn_Ibflex_inh/2', 'autorouting', 'on');
add_line(mdl, 'MN_flex/1', 'syn_Iaext_inh_on_flex/2', 'autorouting', 'on');
% synapse outputs -> MN current sums
add_line(mdl, 'syn_Iaext_exc/1', 'sumMN_ext/2', 'autorouting', 'on');
add_line(mdl, 'syn_Ibext_inh/1', 'sumMN_ext/3', 'autorouting', 'on');
add_line(mdl, 'syn_Iaflex_inh_on_ext/1', 'sumMN_ext/4', 'autorouting', 'on');
add_line(mdl, 'syn_Iaflex_exc/1', 'sumMN_flex/2', 'autorouting', 'on');
add_line(mdl, 'syn_Ibflex_inh/1', 'sumMN_flex/3', 'autorouting', 'on');
add_line(mdl, 'syn_Iaext_inh_on_flex/1', 'sumMN_flex/4', 'autorouting', 'on');

%% ---- muscles ----
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_ext', [720 90 800 140], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_flex', [720 250 800 300], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_ext', [840 80 920 140], 'Fmax', 'Fmax_ext', 'epsMax', '0.25');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_flex', [840 240 920 300], 'Fmax', 'Fmax_flex', 'epsMax', '0.25');
add_line(mdl, 'MN_ext/2', 'Act_ext/1', 'autorouting', 'on');
add_line(mdl, 'MN_flex/2', 'Act_flex/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'BPA_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'BPA_flex/1', 'autorouting', 'on');
add_line(mdl, 'strain_ext_g/1', 'BPA_ext/2', 'autorouting', 'on');
add_line(mdl, 'strain_flex_g/1', 'BPA_flex/2', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'Ib_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'Ib_flex/1', 'autorouting', 'on');

%% ---- torque arms into plant: muscleT = Tflex - Text ----
add_block('simulink/Math Operations/Gain', [mdl '/r_ext'], 'Gain', 'r_arm', 'Position', [960 90 990 120]);
add_block('simulink/Math Operations/Gain', [mdl '/r_flex'], 'Gain', 'r_arm', 'Position', [960 250 990 280]);
add_line(mdl, 'BPA_ext/1', 'r_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'r_flex/1', 'autorouting', 'on');
add_line(mdl, 'r_flex/1', 'muscleT/1', 'autorouting', 'on');
add_line(mdl, 'r_ext/1', 'muscleT/2', 'autorouting', 'on');

%% ---- logging ----
logNames = {'th', 'thd', 'A_ext', 'A_flex', 'V_MN_ext', 'V_MN_flex', 'F_ext', 'F_flex', 'Ia_ext_c', 'Ib_ext_c'};
for k = 1:numel(logNames)
    y = 80 + 80*(k-1);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logNames{k}], 'VariableName', ['log_' logNames{k}], 'SaveFormat', 'Timeseries', 'Position', [1080 y 1150 y+30]);
end
add_line(mdl, 'th_int/1', 'log_th/1', 'autorouting', 'on');
add_line(mdl, 'thd_int/1', 'log_thd/1', 'autorouting', 'on');
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
        'SNS knee reflex demo: non-spiking RC neurons + E/I synapses + antagonist BPAs on 1-DOF knee');
    anno.Position = [40 -75 900 -35];
catch
end
try
    legend_txt = ['Diagram language (Szczecinski 2017 / Rybak / Animatlab): ' ...
        'open circle = neuron; circle "Ia"/"Ib" = afferent; ellipse = muscle.  ' ...
        'Connection markers: solid black circle = INHIBITORY (Esyn < 0), ' ...
        'open triangle = EXCITATORY (Esyn >= 0). Tints: Okabe-Ito (colorblind-safe).'];
    anno2 = Simulink.Annotation(mdl, legend_txt);
    anno2.Position = [40 -25 1150 5];
    anno2.FontSize = 10;
catch
end
% keep flagged block labels out of autorouted wires
for nm = {'muscleT', 'netTorque2', 'spr_defl', 'strain_ext_g', 'sumMN_ext', 'sumMN_flex'}
    try
        set_param([mdl '/' nm{1}], 'NamePlacement', 'alternate');
    catch
    end
end

save_system(mdl);
fprintf('KneeReflexDemo.slx built.\n');

function h = snsInst(mdl, libPath, name, pos, varargin)
    h = add_block(libPath, [mdl '/' name], 'Position', pos, varargin{:});
end
