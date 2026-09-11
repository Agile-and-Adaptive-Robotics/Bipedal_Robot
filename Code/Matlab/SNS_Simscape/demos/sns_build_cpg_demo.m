%% sns_build_cpg_demo.m — build BPACPGLegDemo.slx (half-center CPG + BPA antagonist pair + 1-DOF knee)
%
% Built 2026-09-10. A "simple CPG" in the SNS library language:
%
%   tonic drive -----> RG_ext ------Exc-------------------> Adp_ext (slow)
%        |               |  <------Inh---- (adaptive) ------/   |
%        |               |                                      |
%        |               +--> S(t) -> Act_ext -> BPA_ext -> T_ext
%        |                                                                      mutual
%   tonic drive -----> RG_flex ------Exc------------------> Adp_flex (slow)   inhibition
%                        |  <------Inh---- (adaptive) ------/    |            between
%                        +--> S(t) -> Act_flex -> BPA_flex -> T_flex          RG_ext/RG_flex
%
% Half-center motif (classic SNS/CPG construction, cf. Rybak & Shevtsova,
% SNS-toolbox CPG examples): RG_ext and RG_flex mutually inhibit (Esyn=-72).
% Each RG slowly recruits its own Adp neuron via a slow excitatory synapse;
% the Adp neuron feeds adaptive inhibition BACK onto its RG. When the
% adaptive inhibition accumulates past the tonic drive, the active RG shuts
% off, its partner is released from inhibition, and the pair alternates.
% Period is set mainly by the Adp membrane tau (Cm/Gm) and the loop gains.
%
% Plant: masked 1-DOF knee subsystem (I*thdd = Tflex - Text + Tload
% - b*thd - K*(th-th0)); muscle strain over ROM as in KneeReflexDemo.
% APPEARANCE: same diagram language as the restyled SNS_Library (heavy
% outlines, Animatlab colors); the plant is a single masked block so the
% sheet reads like a circuit, not Simulink math.

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder
mdl = 'BPACPGLegDemo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '10', 'ScreenColor', 'white');

%% ---- parameters stored in model PreLoadFcn (model opens runnable) ----
paramCmd = strjoin({ ...
    '% --- plant ---' ...
    'I_knee = 0.06;      % kg*m^2 below-knee inertia' ...
    'b_knee = 0.5;       % N*m*s/rad joint damping' ...
    'K_knee = 0.5;       % N*m/rad return spring' ...
    'th0    = 0.26;      % rad spring rest angle' ...
    'ROM    = 1.5708;    % rad normalizing range' ...
    'r_arm  = 0.035;     % m BPA moment arm' ...
    'Tload  = 0.3;       % N*m constant flexion load' ...
    'Fmax_ext = 500;     % N' ...
    'Fmax_flex = 450;    % N' ...
    'epsScale = 0.15;    % muscle strain span over ROM' ...
    '% --- CPG (params tuned in tune_cpg2.m ODE prototype, 2026-09-10) ---' ...
    'drive_ext  = 4.0;   % nA tonic drive to RG_ext (free V = -32 mV)' ...
    'drive_flex = 3.6;   % nA tonic drive to RG_flex (asymmetry breaks the latch)' ...
    'g_adp = 15;         % adaptive current at S_adp = 1 (nA)' ...
    'tauAct = 50;        % ms muscle activation' ...
    }, newline);
set_param(mdl, 'PreLoadFcn', paramCmd);
eval(paramCmd);

%% ---- masked plant subsystem: [T_flex, T_ext] -> [theta, eps_ext, eps_flex] ----
plant = [mdl '/KneePlant'];
add_block('simulink/Ports & Subsystems/Subsystem', plant, 'Position', [900 330 1060 470]);
delete_line(plant, 'In1/1', 'Out1/1');
delete_block([plant '/In1']);
delete_block([plant '/Out1']);
add_block('simulink/Sources/In1', [plant '/T_flex'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Sources/In1', [plant '/T_ext'], 'Port', '2', 'Position', [25 178 55 192]);
add_block('simulink/Math Operations/Sum', [plant '/netT'], 'Inputs', '+++', 'Position', [120 110 150 170]);
add_block('simulink/Math Operations/Sum', [plant '/muscleT'], 'Inputs', '+-', 'Position', [60 100 90 140]);
add_block('simulink/Math Operations/Sum', [plant '/negT'], 'Inputs', '--', 'Position', [250 205 280 245]);
add_block('simulink/Math Operations/Gain', [plant '/invI'], 'Gain', '1/I_knee', 'Position', [190 125 220 155]);
add_block('simulink/Continuous/Integrator', [plant '/thd_int'], 'InitialCondition', '0', 'Position', [310 125 340 155]);
add_block('simulink/Continuous/Integrator', [plant '/th_int'], 'InitialCondition', '0.26', 'Position', [370 125 400 155]);
add_block('simulink/Math Operations/Gain', [plant '/damp'], 'Gain', 'b_knee', 'Position', [310 210 340 240]);
add_block('simulink/Math Operations/Gain', [plant '/spring'], 'Gain', 'K_knee', 'Position', [430 210 460 240]);
add_block('simulink/Math Operations/Sum', [plant '/spr_defl'], 'Inputs', '+-', 'Position', [370 215 400 245]);
add_block('simulink/Sources/Constant', [plant '/th0_c'], 'Value', 'th0', 'Position', [310 260 340 290]);
add_block('simulink/Sources/Constant', [plant '/Tload_c'], 'Value', 'Tload', 'Position', [60 160 90 190]);
add_block('simulink/Math Operations/Gain', [plant '/th_norm'], 'Gain', '1/ROM', 'Position', [460 125 490 155]);
add_block('simulink/Math Operations/Sum', [plant '/flex_stretch'], 'Inputs', '+-', 'Position', [520 190 550 220]);
add_block('simulink/Sources/Constant', [plant '/one_c'], 'Value', '1', 'Position', [460 250 490 280]);
add_block('simulink/Math Operations/Gain', [plant '/strain_ext_g'], 'Gain', 'epsScale', 'Position', [580 120 610 150]);
add_block('simulink/Math Operations/Gain', [plant '/strain_flex_g'], 'Gain', 'epsScale', 'Position', [580 190 610 220]);
add_block('simulink/Sinks/Out1', [plant '/theta'], 'Port', '1', 'Position', [650 48 680 62]);
add_block('simulink/Sinks/Out1', [plant '/eps_ext'], 'Port', '2', 'Position', [650 128 680 142]);
add_block('simulink/Sinks/Out1', [plant '/eps_flex'], 'Port', '3', 'Position', [650 198 680 212]);
pl = @(a, b) add_line(plant, a, b, 'autorouting', 'on');
% netT = [T_flex - T_ext] + [-(b*thd + K*(th-th0))] + Tload
pl('T_flex/1', 'muscleT/1');
pl('T_ext/1', 'muscleT/2');
pl('muscleT/1', 'netT/1');
pl('Tload_c/1', 'netT/3');
pl('netT/1', 'invI/1');
pl('invI/1', 'thd_int/1');
pl('thd_int/1', 'th_int/1');
pl('thd_int/1', 'damp/1');
pl('th_int/1', 'spr_defl/1');
pl('th0_c/1', 'spr_defl/2');
pl('spr_defl/1', 'spring/1');
pl('damp/1', 'negT/1');
pl('spring/1', 'negT/2');
pl('negT/1', 'netT/2');
pl('th_int/1', 'theta/1');
pl('th_int/1', 'th_norm/1');
pl('th_norm/1', 'strain_ext_g/1');
pl('th_norm/1', 'flex_stretch/2');
pl('one_c/1', 'flex_stretch/1');
pl('flex_stretch/1', 'strain_flex_g/1');
pl('strain_ext_g/1', 'eps_ext/1');
pl('strain_flex_g/1', 'eps_flex/1');
m = Simulink.Mask.create(plant);
m.Type = 'SNS Knee Plant (1-DOF)';
m.Description = ['Reduced-order 1-DOF knee: I*thdd = T_flex - T_ext + Tload - b*thd - K*(th-th0). ' ...
    'Outputs theta [rad] and normalized muscle strains eps_ext/eps_flex [0..1]. ' ...
    'Stand-in for the Simscape Multibody import of 09_BA_003.'];
m.Display = plantIconCode();
m.addParameter('Name', 'Tload', 'Type', 'edit', 'Prompt', 'Constant load torque (N*m, flexion +)', 'Value', '0.3');
set_param(plant, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
    'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');
% plant initial angle: start slightly flexed like the reflex demo
set_param([plant '/th_int'], 'InitialCondition', '0.26');

%% ---- CPG core: half-center RG pair with adaptive inhibition ----
% RG neurons: fast RC membrane; S spans V in [Thr, Thr+Slope^-1].
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'RG_ext',  [200 60  300 160], ...
    'Vrest', '-52', 'Gm', '0.2', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'RG_flex', [200 330 300 430], ...
    'Vrest', '-52', 'Gm', '0.2', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
% Adaption neurons: SLOW membrane (tau ~ 600 ms) integrating their RG's output.
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'Adp_ext',  [420 60  520 160], ...
    'Vrest', '-52', 'Gm', '0.05', 'Cm', '30', 'Thr', '-45', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'Adp_flex', [420 330 520 430], ...
    'Vrest', '-52', 'Gm', '0.05', 'Cm', '30', 'Thr', '-45', 'Slope', '1');
% mutual inhibition between the RGs (half-center, WTA-strength)
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'inh_ext_on_flex', [330 240 400 300], ...
    'gmax', '0.8', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'inh_flex_on_ext', [330 160 400 220], ...
    'gmax', '0.8', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% slow excitatory RG -> own Adapting cell
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'exc_ext_adp', [340 70 410 120], ...
    'gmax', '0.12', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'exc_flex_adp', [340 340 410 390], ...
    'gmax', '0.12', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
% adaptive current: LINEAR gain on the Adp neuron's S output (tune_cpg2.m v2
% design; a synapse here made the winner latch because its driving force
% (Esyn - V) collapsed at high RG voltage)
add_block('simulink/Math Operations/Gain', [mdl '/gadp_ext'], 'Gain', '-g_adp', 'Position', [300 15 340 45]);
add_block('simulink/Math Operations/Gain', [mdl '/gadp_flex'], 'Gain', '-g_adp', 'Position', [300 445 340 475]);

% input current sums (drive + adaptive + mutual inhibition) -> RG neurons
add_block('simulink/Math Operations/Sum', [mdl '/sumRG_ext'], 'Inputs', '+++', 'Position', [120 85 150 135]);
add_block('simulink/Math Operations/Sum', [mdl '/sumRG_flex'], 'Inputs', '+++', 'Position', [120 355 150 405]);
add_block('simulink/Sources/Constant', [mdl '/drive_ext_c'], 'Value', 'drive_ext', 'Position', [40 105 70 135]);
add_block('simulink/Sources/Constant', [mdl '/drive_flex_c'], 'Value', 'drive_flex', 'Position', [40 375 70 405]);
add_line(mdl, 'drive_ext_c/1', 'sumRG_ext/1', 'autorouting', 'on');
add_line(mdl, 'drive_flex_c/1', 'sumRG_flex/1', 'autorouting', 'on');
add_line(mdl, 'sumRG_ext/1', 'RG_ext/1', 'autorouting', 'on');
add_line(mdl, 'sumRG_flex/1', 'RG_flex/1', 'autorouting', 'on');
% RG V -> synapses (pre sides)
add_line(mdl, 'RG_ext/1', 'inh_ext_on_flex/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'inh_flex_on_ext/1', 'autorouting', 'on');
add_line(mdl, 'RG_ext/1', 'exc_ext_adp/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'exc_flex_adp/1', 'autorouting', 'on');
% Adp S output -> adaptive current gain -> own RG's sum (port 2)
add_line(mdl, 'Adp_ext/2', 'gadp_ext/1', 'autorouting', 'on');
add_line(mdl, 'gadp_ext/1', 'sumRG_ext/2', 'autorouting', 'on');
add_line(mdl, 'Adp_flex/2', 'gadp_flex/1', 'autorouting', 'on');
add_line(mdl, 'gadp_flex/1', 'sumRG_flex/2', 'autorouting', 'on');
% post sides. NonSpikingSynapse ports: 1 = Isyn output (current), 2 = Vpost
% (voltage feedback input). Mutual inhibition: current into the PARTNER's
% sum block (port 3), Vpost feedback from the partner's membrane V.
add_line(mdl, 'inh_ext_on_flex/1', 'sumRG_flex/3', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'inh_ext_on_flex/2', 'autorouting', 'on');
add_line(mdl, 'inh_flex_on_ext/1', 'sumRG_ext/3', 'autorouting', 'on');
add_line(mdl, 'RG_ext/1', 'inh_flex_on_ext/2', 'autorouting', 'on');
% slow excitation onto Adp neurons: Isyn output -> Adp input; Vpost = Adp V
add_line(mdl, 'exc_ext_adp/1', 'Adp_ext/1', 'autorouting', 'on');
add_line(mdl, 'Adp_ext/1', 'exc_ext_adp/2', 'autorouting', 'on');
add_line(mdl, 'exc_flex_adp/1', 'Adp_flex/1', 'autorouting', 'on');
add_line(mdl, 'Adp_flex/1', 'exc_flex_adp/2', 'autorouting', 'on');

%% ---- drive to muscles: RG S -> activation -> BPA force -> plant ----
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_ext',  [600 70 680 120], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_flex', [600 340 680 390], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_ext',  [740 60 820 120], 'Fmax', 'Fmax_ext', 'epsMax', '0.25');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_flex', [740 330 820 390], 'Fmax', 'Fmax_flex', 'epsMax', '0.25');
add_line(mdl, 'RG_ext/2', 'Act_ext/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/2', 'Act_flex/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'BPA_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'BPA_flex/1', 'autorouting', 'on');
% BPA force -> torque -> plant; strain feedback from plant to BPA strain inputs
add_block('simulink/Math Operations/Gain', [mdl '/r_ext'], 'Gain', 'r_arm', 'Position', [840 75 870 105]);
add_block('simulink/Math Operations/Gain', [mdl '/r_flex'], 'Gain', 'r_arm', 'Position', [840 345 870 375]);
add_line(mdl, 'BPA_flex/1', 'r_flex/1', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'r_ext/1', 'autorouting', 'on');
add_line(mdl, 'r_flex/1', 'KneePlant/1', 'autorouting', 'on');
add_line(mdl, 'r_ext/1', 'KneePlant/2', 'autorouting', 'on');
add_line(mdl, 'KneePlant/2', 'BPA_ext/2', 'autorouting', 'on');
add_line(mdl, 'KneePlant/3', 'BPA_flex/2', 'autorouting', 'on');

%% ---- logging ----
logNames = {'th', 'A_ext', 'A_flex', 'V_RG_ext', 'V_RG_flex', 'F_ext', 'F_flex'};
for k = 1:numel(logNames)
    y = 60 + 90*(k-1);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logNames{k}], ...
        'VariableName', ['log_' logNames{k}], 'SaveFormat', 'Timeseries', ...
        'Position', [1150 y 1220 y+30]);
end
add_line(mdl, 'KneePlant/1', 'log_th/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'log_A_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'log_A_flex/1', 'autorouting', 'on');
add_line(mdl, 'RG_ext/1', 'log_V_RG_ext/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'log_V_RG_flex/1', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'log_F_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'log_F_flex/1', 'autorouting', 'on');

%% ---- annotations ----
try
    anno = Simulink.Annotation(mdl, ...
        'SNS CPG leg demo: half-center RG (mutual inhibition + adaptive inhibition) driving antagonist BPAs on a 1-DOF knee');
    anno.Position = [40 -90 1100 -50];
catch
end

save_system(mdl);
fprintf('BPACPGLegDemo.slx built.\n');

%% ---------------- local functions ----------------
function h = snsInst(mdl, libPath, name, pos, varargin)
    h = add_block(libPath, [mdl '/' name], 'Position', pos, varargin{:});
end

function s = plantIconCode()
    % Rounded box, heavy outline (black outer box + slightly smaller light
    % inner box), knee-joint glyph (two links + pivot circle) and label.
    s = strjoin({ ...
        'patch([-1 1 1 -1], [-1 -1 1 1], [0 0 0]);' ...
        'patch([-0.92 0.92 0.92 -0.92], [-0.84 -0.84 0.84 0.84], [0.93 0.93 0.93]);' ...
        'color(''black'');' ...
        'plot([-0.55 -0.15], [0.45 0.05]);' ...
        'plot([-0.15 0.30], [0.05 0.50]);' ...
        'plot([-0.15 -0.15], [0.05 -0.45]);' ...
        't_ = linspace(0, 2*pi, 49);' ...
        'patch(-0.15 + 0.09*cos(t_), 0.05 + 0.09*sin(t_), [1 1 1]);' ...
        'patch(-0.15 + 0.055*cos(t_), 0.05 + 0.055*sin(t_), [0 0 0]);' ...
        'disp(''knee'');' ...
        }, newline);
end
