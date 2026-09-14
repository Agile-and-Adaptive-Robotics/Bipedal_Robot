%% sns_build_deng_cpg.m - Deng 2019 two-layer CPG as separate Simulink files
%
% From Nourse 2023 Fig 6A + Tables A4-A7 (text: Code\MuJoCo_SNS\spinal\
% _nourse2023.txt; lineage Deng 2019 Biomimetics 4(1):21 - Ben's Animatlab
% biped port is the same architecture). Builds:
%
%   SNS_Deng_Library.slx     persistent-Na half-center neuron "HCNeuron"
%                            (Table A4/A5: Cm 5 nF, Gm 1 uS, Vrest -60 mV,
%                            GNa 1.5 uS, ENa 50 mV, m: S 0.2 / E -40 / K 1,
%                            h: S -0.6 / E -60 / K 0.5, tau_h FIXED 350 ms)
%   demos\SNS_Deng_RG.slx    RG layer: HC_ext/HC_flx + IN-laminated mutual
%                            inhibition (HC -2.749-> IN -2.749-> HC, Esyn
%                            -40/-70 mV, window -60..-25); I_stim inport;
%                            V_ext/V_flx outports
%   demos\SNS_Deng_PF.slx    PF layer: hip pair + knee/ankle pair, same
%                            construction; RG->PF weak exc 0.1 uS (window
%                            -60..-40); V_RG_ext/V_RG_flx inports; 4 V outs
%   demos\SNS_Deng_CPGDemo.slx   ONE 10 nA 1 ms pulse (t=0.1 s) -> 20 s
%                            continuous run + hip MNs (PF->MN 2.565/3.632
%                            uS, Esyn -10, window -60..-50, MN Vrest -100)
%                            + Deng Fig 6B activation map
%
% tau_h note: sns_toolbox's tau_h(V) collapses at depolarized V and
% QUENCHES this circuit (verified deng_cpg_ode.py); Animatlab's port uses
% tau_h.max FIXED. Fixed tau -> RG self-runs at 1.94 s period from one
% pulse (numpy reference results\deng_cpg_ref.mat).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);                       % SNS_Simscape
addpath(root);

DENG_HC_PARAMS = { ...
    'Cm',   'Membrane capacitance Cm (nF)', '5'; ...
    'Gm',   'Membrane conductance Gm (uS)', '1'; ...
    'Vrest','Resting potential Vrest (mV)', '-60'; ...
    'GNa',  'Persistent Na conductance GNa (uS)', '1.5'; ...
    'ENa',  'Na reversal ENa (mV)', '50'; ...
    'Sm',   'm-gate slope Sm', '0.2'; ...
    'Em',   'm-gate midpoint Em (mV)', '-40'; ...
    'Km',   'm-gate constant Km', '1'; ...
    'Sh',   'h-gate slope Sh', '-0.6'; ...
    'Eh',   'h-gate midpoint Eh (mV)', '-60'; ...
    'Kh',   'h-gate exponent Kh', '0.5'; ...
    'tauH', 'h-gate time constant tauH (ms), FIXED', '350'; ...
    'h0',   'initial h (h_inf(Vrest) = 0.6667)', '0.6666667'};

%% ============ 1) SNS_Deng_Library: persistent-Na HC neuron =============
lib = 'SNS_Deng_Library';
if bdIsLoaded(lib), close_system(lib, 0); end
if exist(fullfile(root, [lib '.slx']), 'file')
    delete(fullfile(root, [lib '.slx']));
end
new_system(lib, 'Library');
load_system(lib);

blk = [lib '/HCNeuron'];
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 40 160 160]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
% --- membrane ---
add_block('simulink/Sources/In1', [blk '/Iin'], 'Port', '1', 'Position', [25 63 55 77]);
add_block('simulink/Sources/Constant', [blk '/Vrest_c'], 'Value', 'Vrest', 'Position', [90 170 120 200]);
add_block('simulink/Math Operations/Sum', [blk '/leak'], 'Inputs', '-+', 'Position', [160 110 190 140]);
add_block('simulink/Math Operations/Gain', [blk '/Gm_g'], 'Gain', 'Gm', 'Position', [220 110 250 140]);
add_block('simulink/Math Operations/Sum', [blk '/dVsum'], 'Inputs', '+++', 'Position', [290 60 320 120]);
add_block('simulink/Math Operations/Gain', [blk '/memb'], 'Gain', '1000/Cm', 'Position', [350 60 380 100]);
add_block('simulink/Continuous/Integrator', [blk '/Vint'], 'InitialCondition', 'Vrest', 'Position', [410 56 440 104]);
add_block('simulink/Sinks/Out1', [blk '/V_mV'], 'Port', '1', 'Position', [720 66 750 94]);
% --- m gate (instantaneous): m = 1/(1 + Km*exp(Sm*(Em-V))) ---
add_block('simulink/Sources/Constant', [blk '/Em_c'], 'Value', 'Em', 'Position', [90 230 120 260]);
add_block('simulink/Math Operations/Sum', [blk '/mArg'], 'Inputs', '+-', 'Position', [160 210 190 240]);
add_block('simulink/Math Operations/Math Function', [blk '/mExp'], 'Operator', 'exp', 'Position', [210 210 240 240]);
add_block('simulink/Math Operations/Gain', [blk '/Sm_g'], 'Gain', 'Sm', 'Position', [185 208 205 232]);
add_block('simulink/Math Operations/Gain', [blk '/mK'], 'Gain', 'Km', 'Position', [260 210 290 240]);
add_block('simulink/Sources/Constant', [blk '/one_c'], 'Value', '1', 'Position', [260 180 290 200]);
add_block('simulink/Math Operations/Sum', [blk '/mDen'], 'Inputs', '++', 'Position', [310 210 340 240]);
add_block('simulink/Math Operations/Math Function', [blk '/mInv'], 'Operator', 'reciprocal', 'Position', [360 210 390 240]);
% --- h gate: h_inf = 1/(1 + Kh*exp(Sh*(Eh-V)));  dh/dt = (h_inf-h)/tauH ---
add_block('simulink/Sources/Constant', [blk '/Eh_c'], 'Value', 'Eh', 'Position', [90 350 120 380]);
add_block('simulink/Math Operations/Sum', [blk '/hArg'], 'Inputs', '+-', 'Position', [160 330 190 360]);
add_block('simulink/Math Operations/Math Function', [blk '/hExp'], 'Operator', 'exp', 'Position', [210 330 240 360]);
add_block('simulink/Math Operations/Gain', [blk '/Sh_g'], 'Gain', 'Sh', 'Position', [185 328 205 352]);
add_block('simulink/Math Operations/Gain', [blk '/hK'], 'Gain', 'Kh', 'Position', [260 330 290 360]);
add_block('simulink/Sources/Constant', [blk '/one2_c'], 'Value', '1', 'Position', [260 300 290 320]);
add_block('simulink/Math Operations/Sum', [blk '/hDen'], 'Inputs', '++', 'Position', [310 330 340 360]);
add_block('simulink/Math Operations/Math Function', [blk '/hInf'], 'Operator', 'reciprocal', 'Position', [360 330 390 360]);
add_block('simulink/Math Operations/Sum', [blk '/dh'], 'Inputs', '+-', 'Position', [420 330 450 360]);
add_block('simulink/Math Operations/Gain', [blk '/invTauH'], 'Gain', '1000/tauH', 'Position', [470 330 500 360]);
add_block('simulink/Continuous/Integrator', [blk '/hint'], 'InitialCondition', 'h0', 'Position', [530 326 560 364]);
% --- I_Na = GNa * m * h^K * (ENa - V) ---
add_block('simulink/Sources/Constant', [blk '/ENa_c'], 'Value', 'ENa', 'Position', [400 180 430 210]);
add_block('simulink/Math Operations/Sum', [blk '/naDrive'], 'Inputs', '+-', 'Position', [460 150 490 180]);
add_block('simulink/Math Operations/Math Function', [blk '/hPow'], 'Operator', 'pow', 'Position', [520 330 550 360]);
add_block('simulink/Sources/Constant', [blk '/Kh_c'], 'Value', 'Kh', 'Position', [520 390 550 420]);
add_block('simulink/Math Operations/Product', [blk '/naProd'], 'Inputs', '3', 'Position', [600 140 630 200]);
add_block('simulink/Math Operations/Gain', [blk '/GNa_g'], 'Gain', 'GNa', 'Position', [650 150 680 180]);

% wiring
add_line(blk, 'Iin/1', 'dVsum/1');
add_line(blk, 'Vint/1', 'leak/1');
add_line(blk, 'Vrest_c/1', 'leak/2');
add_line(blk, 'leak/1', 'Gm_g/1');
add_line(blk, 'Gm_g/1', 'dVsum/2');
add_line(blk, 'GNa_g/1', 'dVsum/3');
add_line(blk, 'dVsum/1', 'memb/1');
add_line(blk, 'memb/1', 'Vint/1');
add_line(blk, 'Vint/1', 'V_mV/1');
add_line(blk, 'Vint/1', 'mArg/2');
add_line(blk, 'Vint/1', 'hArg/2');
add_line(blk, 'Vint/1', 'naDrive/2');
add_line(blk, 'Em_c/1', 'mArg/1');
add_line(blk, 'mArg/1', 'Sm_g/1');
add_line(blk, 'Sm_g/1', 'mExp/1');
add_line(blk, 'mExp/1', 'mK/1');
add_line(blk, 'mK/1', 'mDen/2');
add_line(blk, 'one_c/1', 'mDen/1');
add_line(blk, 'mDen/1', 'mInv/1');
add_line(blk, 'Eh_c/1', 'hArg/1');
add_line(blk, 'hArg/1', 'Sh_g/1');
add_line(blk, 'Sh_g/1', 'hExp/1');
add_line(blk, 'hExp/1', 'hK/1');
add_line(blk, 'hK/1', 'hDen/2');
add_line(blk, 'one2_c/1', 'hDen/1');
add_line(blk, 'hDen/1', 'hInf/1');
add_line(blk, 'hInf/1', 'dh/1');
add_line(blk, 'hint/1', 'dh/2');
add_line(blk, 'dh/1', 'invTauH/1');
add_line(blk, 'invTauH/1', 'hint/1');
add_line(blk, 'ENa_c/1', 'naDrive/1');
add_line(blk, 'mInv/1', 'naProd/1');
add_line(blk, 'hint/1', 'hPow/1');
add_line(blk, 'Kh_c/1', 'hPow/2');
add_line(blk, 'hPow/1', 'naProd/2');
add_line(blk, 'naDrive/1', 'naProd/3');
add_line(blk, 'naProd/1', 'GNa_g/1');

m = Simulink.Mask.create(blk);
m.Description = ['Deng/Nourse persistent-Na half-center neuron (Tables A4/A5).' newline ...
    'Cm dV/dt = Gm(Vrest-V) + Iin + GNa*m*h^K*(ENa-V); m instantaneous,' newline ...
    'dh/dt = (h_inf - h)/tauH with tauH FIXED (Animatlab interpretation).' newline ...
    'Nourse 2023 defaults reproduce the rat hindlimb RG half-center.'];
for k = 1:size(DENG_HC_PARAMS, 1)
    m.addParameter('Name', DENG_HC_PARAMS{k, 1}, 'Type', 'edit', ...
        'Prompt', DENG_HC_PARAMS{k, 2}, 'Value', DENG_HC_PARAMS{k, 3});
end
set_param(blk, 'MaskIconFrame', 'off');
save_system(lib);
fprintf('built SNS_Deng_Library.slx (HCNeuron)\n');

%% ============ 2) RG layer ==============================================
build_hc_layer('SNS_Deng_RG', fullfile(here, 'SNS_Deng_RG.slx'), ...
               {{'', 'stim'}});                 % one pair, stimulus input
fprintf('built demos/SNS_Deng_RG.slx\n');

%% ============ 3) PF layer ==============================================
build_hc_layer('SNS_Deng_PF', fullfile(here, 'SNS_Deng_PF.slx'), ...
               {{'hip'}, {'ka'}});              % two pairs, RG-driven
fprintf('built demos/SNS_Deng_PF.slx\n');

%% ============ 4) demo: one pulse -> continuous run =====================
mdl = 'SNS_Deng_CPGDemo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(fullfile(here, [mdl '.slx']), 'file')
    delete(fullfile(here, [mdl '.slx']));
end
new_system(mdl);
load_system('SNS_Deng_RG');
load_system('SNS_Deng_PF');
add_block('simulink/Ports & Subsystems/Model', [mdl '/RG'], ...
          'ModelFile', 'SNS_Deng_RG.slx', 'Position', [300 100 400 200]);
add_block('simulink/Ports & Subsystems/Model', [mdl '/PF'], ...
          'ModelFile', 'SNS_Deng_PF.slx', 'Position', [560 80 680 240]);
add_block('simulink/Sources/Pulse Generator', [mdl '/kick'], ...
          'PulseType', 'Time based', 'Amplitude', '10', 'Period', '40', ...
          'PulseWidth', '0.05', 'PhaseDelay', '0.1', ...
          'Position', [80 130 140 170]);   % 10 nA, 20 ms kick at t = 0.1 s
add_line(mdl, 'kick/1', 'RG/1');
add_line(mdl, 'RG/1', 'PF/1');
add_line(mdl, 'RG/2', 'PF/2');
% hip MNs + Deng Fig 6B activation map
mn_specs = {'MN_hip_ext', 2.565, 1; 'MN_hip_flx', 3.632, 2};
for c = 1:size(mn_specs, 1)
    nm = mn_specs{c, 1};
    y0 = 40 + 120 * (c - 1);
    add_block('SNS_Library/NonSpikingSynapse', [mdl '/syn_' nm], ...
        'gmax', num2str(mn_specs{c, 2}), 'Esyn', '-10', ...
        'ThrPre', '-60', 'SlopePre', '10', ...
        'Position', [700 y0 760 y0 + 60]);
    add_block('SNS_Library/NonSpikingNeuron', [mdl '/' nm], ...
        'Vrest', '-100', 'Gm', '1', 'Cm', '5', ...
        'Thr', '-60', 'Slope', '10', ...
        'Position', [790 y0 870 y0 + 80]);
    add_block('simulink/User-Defined Functions/Fcn', [mdl '/act_' nm], ...
        'Expression', '1/(1+exp(0.1532*(-70-u)))-0.01', ...
        'Position', [920 y0 1010 y0 + 40]);
    add_line(mdl, sprintf('PF/%d', mn_specs{c, 3}), sprintf('syn_%s/1', nm));
    add_line(mdl, [nm '/1'], sprintf('syn_%s/2', nm));
    add_line(mdl, sprintf('syn_%s/1', nm), [nm '/1']);
    add_line(mdl, [nm '/1'], sprintf('act_%s/1', nm));
end
% logging
logit = {'RG', 1, 'V_RG_ext'; 'RG', 2, 'V_RG_flx'; ...
         'PF', 1, 'V_PF_hip_e'; 'PF', 2, 'V_PF_hip_f'; ...
         'PF', 3, 'V_PF_ka_e'; 'PF', 4, 'V_PF_ka_f'; ...
         'MN_hip_ext', 1, 'V_MN_ext'; 'MN_hip_flx', 1, 'V_MN_flx'; ...
         'act_MN_hip_ext', 1, 'act_ext'; 'act_MN_hip_flx', 1, 'act_flx'};
for k = 1:size(logit, 1)
    ph = get_param([mdl '/' logit{k, 1}], 'PortHandles');
    set_param(ph.Outport(logit{k, 2}), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', 'DataLoggingName', logit{k, 3});
end
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '1e-4', 'StopTime', '20');
save_system(mdl, fullfile(here, [mdl '.slx']));
fprintf('built demos/SNS_Deng_CPGDemo.slx (one 10 nA/1 ms pulse, 20 s)\n');

%% ============ local functions ==========================================
function build_hc_layer(name, slxPath, jointSpecs)
% One model with a half-center pair per joint spec.
% jointSpecs: {{'', 'stim'}}  -> RG layer (names HC_*/IN_*, I_stim inport)
%             {{'hip'}, {'ka'}} -> PF layer (names hip_HC_* etc., RG inports)
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
if bdIsLoaded(name), close_system(name, 0); end
if exist(slxPath, 'file'), delete(slxPath); end
new_system(name);
load_system('SNS_Deng_Library');
load_system('SNS_Library');

isRG = numel(jointSpecs{1}) > 1 && strcmp(jointSpecs{1}{2}, 'stim');
WIN = {'-60', '35'};                          % HC<->IN window -60..-25

for s = 1:numel(jointSpecs)
    j = jointSpecs{s}{1};
    y = 60 + 340 * (s - 1);                   % row per pair
    pre = j;                                  % 'hip' or '' (RG)
    add_block('SNS_Deng_Library/HCNeuron', ...
              [name '/pair' num2str(s) '_HC_ext'], 'Position', [640 y 740 y + 100]);
    add_block('SNS_Deng_Library/HCNeuron', ...
              [name '/pair' num2str(s) '_HC_flx'], 'Position', [640 y + 150 740 y + 250]);
    add_block('SNS_Library/NonSpikingNeuron', ...
              [name '/pair' num2str(s) '_IN_ext'], ...
        'Vrest', '-60', 'Gm', '1', 'Cm', '5', 'Thr', WIN{1}, 'Slope', WIN{2}, ...
        'Position', [440 y 540 y + 100]);
    add_block('SNS_Library/NonSpikingNeuron', ...
              [name '/pair' num2str(s) '_IN_flx'], ...
        'Vrest', '-60', 'Gm', '1', 'Cm', '5', 'Thr', WIN{1}, 'Slope', WIN{2}, ...
        'Position', [440 y + 150 540 y + 250]);
    % HC->IN exc and IN->HC inh
    synspecs = {'s1', 'HC_ext', 'IN_ext', '-40'; 's2', 'IN_ext', 'HC_flx', '-70'; ...
                's3', 'HC_flx', 'IN_flx', '-40'; 's4', 'IN_flx', 'HC_ext', '-70'};
    for q = 1:4
        sb = sprintf('pair%d_%s', s, synspecs{q, 1});
        add_block('SNS_Library/NonSpikingSynapse', [name '/' sb], ...
            'gmax', '2.749', 'Esyn', synspecs{q, 4}, ...
            'ThrPre', WIN{1}, 'SlopePre', WIN{2}, ...
            'Position', [300 y + (q - 1) * 55 360 y + (q - 1) * 55 + 40]);
        add_line(name, ['pair' num2str(s) '_' synspecs{q, 2} '/1'], [sb '/1']);
        add_line(name, ['pair' num2str(s) '_' synspecs{q, 3} '/1'], [sb '/2']);
        % s2/s4 (IN->HC) route through the HC current Sums, not directly
        if q == 1 || q == 3
            add_line(name, [sb '/1'], ['pair' num2str(s) '_' synspecs{q, 3} '/1']);
        end
    end
    % merge all currents into each HC (Sum -> Iin)
    add_block('simulink/Math Operations/Sum', ...
              [name '/pair' num2str(s) '_sum_ext'], 'Inputs', '++', ...
              'Position', [590 y + 260 615 y + 290]);
    add_block('simulink/Math Operations/Sum', ...
              [name '/pair' num2str(s) '_sum_flx'], 'Inputs', '++', ...
              'Position', [590 y + 300 615 y + 330]);
    add_line(name, ['pair' num2str(s) '_s4/1'], ['pair' num2str(s) '_sum_ext/1']);
    add_line(name, ['pair' num2str(s) '_s2/1'], ['pair' num2str(s) '_sum_flx/1']);
    add_line(name, ['pair' num2str(s) '_sum_ext/1'], ['pair' num2str(s) '_HC_ext/1']);
    add_line(name, ['pair' num2str(s) '_sum_flx/1'], ['pair' num2str(s) '_HC_flx/1']);
    % stimulus (RG) or RG-drive synapses (PF)
    if isRG
        add_block('simulink/Sources/In1', [name '/I_stim'], 'Port', '1', ...
                  'Position', [500 y + 300 530 y + 330]);
        add_line(name, 'I_stim/1', ['pair' num2str(s) '_sum_ext/2']);
        extName = 'HC_ext'; flxName = 'HC_flx';
    else
        if s == 1     % RG-drive inports exist once for the whole PF layer
            add_block('simulink/Sources/In1', [name '/V_RG_ext'], 'Port', '1', ...
                      'Position', [20 60 50 90]);
            add_block('simulink/Sources/In1', [name '/V_RG_flx'], 'Port', '2', ...
                      'Position', [20 360 50 390]);
        end
        for half = 1:2
            if half == 1, rh = 'ext'; srcp = 'V_RG_ext'; else, rh = 'flx'; srcp = 'V_RG_flx'; end
            sb = ['pair' num2str(s) '_s_RG_' rh];
            add_block('SNS_Library/NonSpikingSynapse', [name '/' sb], ...
                'gmax', '0.1', 'Esyn', '-40', 'ThrPre', '-60', 'SlopePre', '20', ...
                'Position', [100 y + 250 * (half - 1) 160 y + 310 * (half - 1) + 60]);
            add_line(name, [srcp '/1'], [sb '/1']);
            add_line(name, ['pair' num2str(s) '_HC_' rh '/1'], [sb '/2']);
            add_line(name, [sb '/1'], ...
                     ['pair' num2str(s) '_sum_' rh '/2']);
        end
        extName = [j '_HC_ext']; flxName = [j '_HC_flx'];
    end
    % outports (per pair)
    add_block('simulink/Sinks/Out1', [name '/V_' extName], ...
              'Port', num2str(2 * s - 1), 'Position', [800 y 830 y + 30]);
    add_block('simulink/Sinks/Out1', [name '/V_' flxName], ...
              'Port', num2str(2 * s), 'Position', [800 y + 150 830 y + 180]);
    add_line(name, ['pair' num2str(s) '_HC_ext/1'], ['V_' extName '/1']);
    add_line(name, ['pair' num2str(s) '_HC_flx/1'], ['V_' flxName '/1']);
end
set_param(name, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '1e-4');
save_system(name, slxPath);
close_system(name, 0);
end
