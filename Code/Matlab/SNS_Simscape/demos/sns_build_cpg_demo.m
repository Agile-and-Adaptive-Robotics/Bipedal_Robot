%% sns_build_cpg_demo.m — build BPACPGLegDemo.slx (half-center CPG + BPA antagonist pair + 1-DOF knee)
%
% Built 2026-09-10; ported to the 2026-09-22 library architecture (one-input
% synapses onto neuron syn ports, summation inside the neurons, shared
% KneeModel subsystem copied from KneeReflexDemo).
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
% the Adp neuron feeds adaptive inhibition BACK onto its RG (LINEAR gain on
% the Adp S output — a synapse there latches the winner because its driving
% force collapses at high RG voltage). The drive and the adaptive current
% both land on the RG's Iapp port (the neuron sums a vector Iapp internally).
%
% Model: masked 1-DOF KneeModel subsystem (same block as KneeReflexDemo).

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder
mdl = 'BPACPGLegDemo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '10', 'ScreenColor', 'white', ...
    'UnconnectedInputMsg', 'none');
load_system(fullfile(cdto, 'KneeReflexDemo'));   % source of the shared KneeModel

%% ---- parameters stored in model PreLoadFcn (model opens runnable) ----
paramCmd = strjoin({ ...
    '% --- knee model ---' ...
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

%% ---- knee model subsystem: 1:1 copy of KneeReflexDemo's KneeModel ----
add_block('KneeReflexDemo/KneeModel', [mdl '/KneeModel'], 'Position', [880 330 1060 470]);
set_param([mdl '/KneeModel'], 'Name', 'KneeModel');

%% ---- CPG core: half-center RG pair with adaptive inhibition ----
% RG neurons: fast RC membrane; S spans V in [Thr, Thr+Slope^-1].
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'RG_ext',  [240 60  340 180], ...
    'Vrest', '-52', 'Gm', '0.2', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'RG_flex', [240 330 340 450], ...
    'Vrest', '-52', 'Gm', '0.2', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
% Adaption neurons: SLOW membrane (tau ~ 600 ms) integrating their RG's output.
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'Adp_ext',  [520 60  620 180], ...
    'Vrest', '-52', 'Gm', '0.05', 'Cm', '30', 'Thr', '-45', 'Slope', '1');
snsInst(mdl, 'SNS_Library/NonSpikingNeuron', 'Adp_flex', [520 330 620 450], ...
    'Vrest', '-52', 'Gm', '0.05', 'Cm', '30', 'Thr', '-45', 'Slope', '1');
% mutual inhibition between the RGs (half-center, WTA-strength) — placed
% against the RG each one synapses ONTO
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'inh_ext_on_flex', [170 355 210 387], ...
    'gmax', '0.8', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'inh_flex_on_ext', [170 85 210 117], ...
    'gmax', '0.8', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% slow excitatory RG -> own Adapting cell (against the Adp neuron)
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'exc_ext_adp', [450 80 490 112], ...
    'gmax', '0.12', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst(mdl, 'SNS_Library/NonSpikingSynapse', 'exc_flex_adp', [450 350 490 382], ...
    'gmax', '0.12', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
% adaptive current: LINEAR gain on the Adp neuron's S output (tune_cpg2.m v2
% design; a synapse here made the winner latch because its driving force
% (Esyn - V) collapsed at high RG voltage)
add_block('simulink/Math Operations/Gain', [mdl '/gadp_ext'], 'Gain', '-g_adp', 'Position', [380 25 420 55]);
add_block('simulink/Math Operations/Gain', [mdl '/gadp_flex'], 'Gain', '-g_adp', 'Position', [380 465 420 495]);
% drive + adaptive current bundled onto the RG Iapp port (neuron sums vectors)
add_block('simulink/Signal Routing/Mux', [mdl '/bundle_ext'], 'Inputs', '2', 'Position', [150 95 153 145]);
add_block('simulink/Signal Routing/Mux', [mdl '/bundle_flex'], 'Inputs', '2', 'Position', [150 365 153 415]);
add_block('simulink/Sources/Constant', [mdl '/drive_ext_c'], 'Value', 'drive_ext', 'Position', [60 100 100 130]);
add_block('simulink/Sources/Constant', [mdl '/drive_flex_c'], 'Value', 'drive_flex', 'Position', [60 370 100 400]);
for nm = {'bundle_ext', 'bundle_flex', 'gadp_ext', 'gadp_flex'}
    try, set_param([mdl '/' nm{1}], 'ShowName', 'off'); catch, end
end
add_line(mdl, 'drive_ext_c/1', 'bundle_ext/1', 'autorouting', 'on');
add_line(mdl, 'gadp_ext/1', 'bundle_ext/2', 'autorouting', 'on');
add_line(mdl, 'drive_flex_c/1', 'bundle_flex/1', 'autorouting', 'on');
add_line(mdl, 'gadp_flex/1', 'bundle_flex/2', 'autorouting', 'on');
add_line(mdl, 'bundle_ext/1', 'RG_ext/1', 'autorouting', 'on');
add_line(mdl, 'bundle_flex/1', 'RG_flex/1', 'autorouting', 'on');
% mutual inhibition: pre = partner RG V, output -> own syn1 port
add_line(mdl, 'RG_ext/1', 'inh_ext_on_flex/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'inh_flex_on_ext/1', 'autorouting', 'on');
add_line(mdl, 'inh_ext_on_flex/1', 'RG_flex/2', 'autorouting', 'on');
add_line(mdl, 'inh_flex_on_ext/1', 'RG_ext/2', 'autorouting', 'on');
% slow excitation onto Adp neurons: pre = own RG V, output -> Adp syn1
add_line(mdl, 'RG_ext/1', 'exc_ext_adp/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'exc_flex_adp/1', 'autorouting', 'on');
add_line(mdl, 'exc_ext_adp/1', 'Adp_ext/2', 'autorouting', 'on');
add_line(mdl, 'exc_flex_adp/1', 'Adp_flex/2', 'autorouting', 'on');
% Adp S output -> adaptive current gain -> own RG's Iapp bundle
add_line(mdl, 'Adp_ext/2', 'gadp_ext/1', 'autorouting', 'on');
add_line(mdl, 'Adp_flex/2', 'gadp_flex/1', 'autorouting', 'on');

%% ---- drive to muscles: RG S -> activation -> BPA force -> knee model ----
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_ext',  [680 90 740 150], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/MuscleActivation', 'Act_flex', [680 360 740 420], 'tauAct', '50');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_ext',  [780 90 850 150], 'Fmax', 'Fmax_ext', 'epsMax', '0.25');
snsInst(mdl, 'SNS_Library/BPAForce', 'BPA_flex', [780 360 850 420], 'Fmax', 'Fmax_flex', 'epsMax', '0.25');
add_line(mdl, 'RG_ext/2', 'Act_ext/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/2', 'Act_flex/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'BPA_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'BPA_flex/1', 'autorouting', 'on');
add_block('simulink/Math Operations/Gain', [mdl '/r_ext'], 'Gain', 'r_arm', 'Position', [870 105 900 135]);
add_block('simulink/Math Operations/Gain', [mdl '/r_flex'], 'Gain', 'r_arm', 'Position', [870 375 900 405]);
add_line(mdl, 'BPA_flex/1', 'r_flex/1', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'r_ext/1', 'autorouting', 'on');
add_line(mdl, 'r_flex/1', 'KneeModel/1', 'autorouting', 'on');
add_line(mdl, 'r_ext/1', 'KneeModel/2', 'autorouting', 'on');
add_line(mdl, 'KneeModel/7', 'BPA_ext/2', 'autorouting', 'on');
add_line(mdl, 'KneeModel/8', 'BPA_flex/2', 'autorouting', 'on');

%% ---- logging ----
logNames = {'th', 'A_ext', 'A_flex', 'V_RG_ext', 'V_RG_flex', 'F_ext', 'F_flex'};
for k = 1:numel(logNames)
    y = 60 + 90*(k-1);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logNames{k}], ...
        'VariableName', ['log_' logNames{k}], 'SaveFormat', 'Timeseries', ...
        'Position', [1150 y 1220 y+30]);
end
add_line(mdl, 'KneeModel/1', 'log_th/1', 'autorouting', 'on');
add_line(mdl, 'Act_ext/1', 'log_A_ext/1', 'autorouting', 'on');
add_line(mdl, 'Act_flex/1', 'log_A_flex/1', 'autorouting', 'on');
add_line(mdl, 'RG_ext/1', 'log_V_RG_ext/1', 'autorouting', 'on');
add_line(mdl, 'RG_flex/1', 'log_V_RG_flex/1', 'autorouting', 'on');
add_line(mdl, 'BPA_ext/1', 'log_F_ext/1', 'autorouting', 'on');
add_line(mdl, 'BPA_flex/1', 'log_F_flex/1', 'autorouting', 'on');

%% ---- annotations ----
try
    anno = Simulink.Annotation(mdl, ...
        'SNS CPG leg demo: half-center RG (mutual inhibition + adaptive inhibition) driving antagonist BPAs on a 1-DOF knee model');
    anno.Position = [40 -90 1100 -50];
catch
end

save_system(mdl);
close_system('KneeReflexDemo', 0);
fprintf('BPACPGLegDemo.slx built (2026-09-22 architecture).\n');

%% ---------------- local functions ----------------
function h = snsInst(mdl, libPath, name, pos, varargin)
    h = add_block(libPath, [mdl '/' name], 'Position', pos, varargin{:});
end
