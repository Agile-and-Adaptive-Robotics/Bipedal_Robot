function build_sns_w2l_cpg_20260930()
% SNS_W2L_CPG - REPRESENTATIVE CORE of the W2L contact-driven bilateral CPG
% built from SNS_Library blocks (laptop R2025b, 2026-09-30 campaign).
%
% *** THIS IS THE REPRESENTATIVE CORE, NOT THE FULL NET. *** The full
% transcribed net (Code\MuJoCo_SNS\spinal\w2l_cpg\build_w2l_net.py) has 95
% neurons / 208 synapses. This model realizes, per side:
%   - RG E/F half-centers as persistent-Na neurons (NaP current with FIXED
%     tau_h = 250 ms, the AnimatLab/Deng semantics the numpy backend
%     SNS_NumpyFixedTau implements) - custom masked block 'NapHC'
%   - laminated mutual inhibition (E->InE exc, InE->F inh, F->InF exc,
%     InF->E inh, g 4.0) + weak direct HC<->HC escape excitation (g 0.5)
%   - crossed commissural inhibition: c1 (RG_flx->c1 exc 1.0, c1->contra
%     RG_flx inh 3.0) and V3 (RG_ext->V3 exc 1.0, V3->contra InE exc 0.08)
%   - PF E/F half-centers (g: RG->PF 2.4, PF cross-IN laminated 4.0)
%   - 2 representative MNs (PF->MN g 2.0)
%   - DRIVE: tonic 2 nA (E) / 3 nA (F) + the .aproj 10 nA/10 ms antiphase
%     kickoff pair (L RG ext, R RG flx)
%   - heel-contact input channel (graded contact neuron tau 40 ms + g 1.0
%     contact synapses onto the ipsilateral extensor layers RG/PF/MN),
%     alternating 1 Hz trains at 20 nA, L phase 0 / R phase 0.5 s.
% MISSING vs the full net: Renshaw cells, all Ia/II/Ib afferent chains,
% the 2nd PF pair per side (hip+knee in the full net; 1 here), toe contact,
% 8 of the 12 MN outputs, IaIN reciprocal machinery.
%
% Units / equations mirror sns_toolbox 1.5.2 as used by the numpy reference:
%   Cm dV/dt = Gm(Vrest-V) + Iapp + sum(g*Esyn) - V*sum(g) + GNa*m*h*(ENa-V)
%   m = 1/(1+kM*exp(sM*(eM-V))) instantaneous; dh/dt = (hInf-h)/tauH, FIXED
%   tauH; hInf = 1/(1+kH*exp(sH*(eH-V))). NaP params = spinal/params.py NAP
%   verbatim (g 12 uS, E 8 mV, m: k1/s0.8/e2; h: k1/s-2/e3.5; tau 0.25 s).
%   Membrane tau: RG/IN/C/V3 50 ms, PF 80, PF-IN 80, MN 30, contact 40.
%
% VALIDATION (this file runs it): 12 s fixed-step Euler ode1 @ 2 ms (= the
% numpy dt), tonic + kickoff + alternating heel trains, i.e. the smoke_w2l.py
% protocol. PASS = both RG-E burst periods within 20% of the numpy reference
% (w2l_numpy_ref.json, full 95-neuron net, same protocol) AND L/R RG-E
% antiphase Pearson r < -0.5 in both.
%
% PROVEN 2026-10-01: the full net does NOT oscillate under tonic+kickoff
% only (R RG-E latches at -1.154 mV max over 2-12 s, zero bursts) - the W2L
% architecture is contact-driven, so the heel trains are part of the drive.

here = fileparts(mfilename('fullpath'));      % ...\SNS_Simscape\dev
root = fileparts(here);                       % ...\SNS_Simscape
addpath(root);
load_system('SNS_Library');

% ---- W2L gains (build_w2l_net.py W2L_GAINS, verbatim) ----------------------
EXC = 8; INH = -5;                            % E_REV_EXC / E_REV_INH (mV)
g_lam = 4.0;      g_dir = 0.5;                % rg_laminate / rg_direct_exc
g_c1pre = 1.0;    g_c1 = 3.0;                 % comm_c1 / c1_inh
g_v3pre = 1.0;    g_v3 = 0.08;                % comm_v3 / v3_weak
g_pfdrive = 2.4;  g_pfcross = 4.0;            % pf_drive / pf_cross
g_pfmn = 2.0;     g_contact = 1.0;            % pf_to_mn / contact
tonicE = 2; tonicF = 3; kickAmp = 10; kickT = 0.01;   % nA, s (.aproj pair)
heelI = 20; heelPeriod = 1.0; heelWidth = 30;         % nA, s, % of period
T_END = 12; SKIP = 2;                                  % s (numpy ref = same)

mdl = 'SNS_W2L_CPG';
if bdIsLoaded(mdl), close_system(mdl, 0); end
slxDst = fullfile(root, 'results', [mdl '.slx']);
if exist(slxDst, 'file'), delete(slxDst); end
new_system(mdl);
set_param(mdl, 'UnconnectedInputMsg', 'none');

%% ---------------- 1) NapHC prototype (persistent-Na half-center) -----------
proto = [mdl '/NapHC_proto'];
add_block('simulink/Ports & Subsystems/Subsystem', proto, 'Position', [40 40 140 140]);
delete_line(proto, 'In1/1', 'Out1/1');
delete_block([proto '/In1']); delete_block([proto '/Out1']);
NS = 4;                                       % syn1..syn4 ports
add_block('simulink/Sources/In1', [proto '/Iapp'], 'Port', '1', 'Position', [25 48 55 62]);
add_block('simulink/Math Operations/Sum', [proto '/IappSum'], 'Inputs', '+', 'Position', [95 43 125 67]);
for k = 1:NS
    add_block('simulink/Sources/In1', [proto '/syn' num2str(k)], 'Port', num2str(k+1), ...
        'PortDimensions', '2', 'Position', [25 78+30*(k-1) 55 92+30*(k-1)]);
    add_block('simulink/Signal Routing/Demux', [proto '/D' num2str(k)], 'Outputs', '2', ...
        'Position', [95 72+30*(k-1) 98 106+30*(k-1)]);
    add_line(proto, ['syn' num2str(k) '/1'], ['D' num2str(k) '/1'], 'autorouting', 'on');
end
add_block('simulink/Math Operations/Sum', [proto '/sumG'], ...
    'Inputs', repmat('+', 1, NS), 'Position', [170 90 200 90+26*NS]);
add_block('simulink/Math Operations/Sum', [proto '/sumE'], ...
    'Inputs', repmat('+', 1, NS), 'Position', [170 90+26*NS+40 200 90+52*NS+40]);
for k = 1:NS
    add_line(proto, ['D' num2str(k) '/1'], ['sumG/' num2str(k)], 'autorouting', 'on');
    add_line(proto, ['D' num2str(k) '/2'], ['sumE/' num2str(k)], 'autorouting', 'on');
end
% Isyn = sumE - sumG*V
add_block('simulink/Math Operations/Product', [proto '/gTimesV'], 'Inputs', '2', 'Position', [270 95 300 125]);
add_block('simulink/Math Operations/Gain', [proto '/negGv'], 'Gain', '-1', 'Position', [330 100 360 130]);
add_block('simulink/Math Operations/Sum', [proto '/Isyn'], 'Inputs', '++', 'Position', [400 110 430 140]);
add_line(proto, 'sumG/1', 'gTimesV/1', 'autorouting', 'on');
add_line(proto, 'gTimesV/1', 'negGv/1', 'autorouting', 'on');
add_line(proto, 'sumE/1', 'Isyn/1', 'autorouting', 'on');
add_line(proto, 'negGv/1', 'Isyn/2', 'autorouting', 'on');
% leak Gm*(Vrest-V)
add_block('simulink/Sources/Constant', [proto '/Vrest_c'], 'Value', 'Vrest', 'Position', [330 240 360 270]);
add_block('simulink/Math Operations/Sum', [proto '/dVm'], 'Inputs', '-+', 'Position', [400 245 430 275]);
add_block('simulink/Math Operations/Gain', [proto '/Gm'], 'Gain', 'Gm', 'Position', [460 250 490 280]);
add_line(proto, 'Vrest_c/1', 'dVm/2', 'autorouting', 'on');
add_line(proto, 'dVm/1', 'Gm/1', 'autorouting', 'on');
% m gate: mInf = 1/(1+kM*exp(sM*(eM-V)))
add_block('simulink/Sources/Constant', [proto '/eM_c'], 'Value', 'eM', 'Position', [90 210 120 240]);
add_block('simulink/Math Operations/Sum', [proto '/mArg'], 'Inputs', '+-', 'Position', [160 210 190 240]);
add_block('simulink/Math Operations/Gain', [proto '/sM_g'], 'Gain', 'sM', 'Position', [210 208 230 232]);
add_block('simulink/Math Operations/Math Function', [proto '/mExp'], 'Operator', 'exp', 'Position', [250 210 280 240]);
add_block('simulink/Math Operations/Gain', [proto '/mK'], 'Gain', 'kM', 'Position', [300 210 330 240]);
add_block('simulink/Sources/Constant', [proto '/one_m'], 'Value', '1', 'Position', [300 180 320 200]);
add_block('simulink/Math Operations/Sum', [proto '/mDen'], 'Inputs', '++', 'Position', [360 200 390 240]);
add_block('simulink/Math Operations/Math Function', [proto '/mInf'], 'Operator', 'reciprocal', 'Position', [420 205 450 235]);
add_line(proto, 'eM_c/1', 'mArg/1', 'autorouting', 'on');
add_line(proto, 'mArg/1', 'sM_g/1', 'autorouting', 'on');
add_line(proto, 'sM_g/1', 'mExp/1', 'autorouting', 'on');
add_line(proto, 'mExp/1', 'mK/1', 'autorouting', 'on');
add_line(proto, 'one_m/1', 'mDen/1', 'autorouting', 'on');
add_line(proto, 'mK/1', 'mDen/2', 'autorouting', 'on');
add_line(proto, 'mDen/1', 'mInf/1', 'autorouting', 'on');
% h gate: hInf = 1/(1+kH*exp(sH*(eH-V))); dh = (hInf-h)/tauH, tauH FIXED
add_block('simulink/Sources/Constant', [proto '/eH_c'], 'Value', 'eH', 'Position', [90 330 120 360]);
add_block('simulink/Math Operations/Sum', [proto '/hArg'], 'Inputs', '+-', 'Position', [160 330 190 360]);
add_block('simulink/Math Operations/Gain', [proto '/sH_g'], 'Gain', 'sH', 'Position', [210 328 230 352]);
add_block('simulink/Math Operations/Math Function', [proto '/hExp'], 'Operator', 'exp', 'Position', [250 330 280 360]);
add_block('simulink/Math Operations/Gain', [proto '/hK'], 'Gain', 'kH', 'Position', [300 330 330 360]);
add_block('simulink/Sources/Constant', [proto '/one_h'], 'Value', '1', 'Position', [300 300 320 320]);
add_block('simulink/Math Operations/Sum', [proto '/hDen'], 'Inputs', '++', 'Position', [360 320 390 360]);
add_block('simulink/Math Operations/Math Function', [proto '/hInf'], 'Operator', 'reciprocal', 'Position', [420 325 450 355]);
add_block('simulink/Math Operations/Sum', [proto '/dh'], 'Inputs', '+-', 'Position', [500 330 530 360]);
add_block('simulink/Math Operations/Gain', [proto '/invTauH'], 'Gain', '1000/tauH', 'Position', [560 330 590 360]);
add_block('simulink/Continuous/Integrator', [proto '/hint'], 'InitialCondition', 'h0', 'Position', [620 326 650 364]);
add_line(proto, 'eH_c/1', 'hArg/1', 'autorouting', 'on');
add_line(proto, 'hArg/1', 'sH_g/1', 'autorouting', 'on');
add_line(proto, 'sH_g/1', 'hExp/1', 'autorouting', 'on');
add_line(proto, 'hExp/1', 'hK/1', 'autorouting', 'on');
add_line(proto, 'one_h/1', 'hDen/1', 'autorouting', 'on');
add_line(proto, 'hK/1', 'hDen/2', 'autorouting', 'on');
add_line(proto, 'hDen/1', 'hInf/1', 'autorouting', 'on');
add_line(proto, 'hInf/1', 'dh/1', 'autorouting', 'on');
add_line(proto, 'hint/1', 'dh/2', 'autorouting', 'on');
add_line(proto, 'dh/1', 'invTauH/1', 'autorouting', 'on');
add_line(proto, 'invTauH/1', 'hint/1', 'autorouting', 'on');
% INa = GNa*mInf*h*(ENa-V)
add_block('simulink/Sources/Constant', [proto '/ENa_c'], 'Value', 'ENa', 'Position', [500 240 530 270]);
add_block('simulink/Math Operations/Sum', [proto '/naDrive'], 'Inputs', '+-', 'Position', [560 245 590 275]);
add_block('simulink/Math Operations/Product', [proto '/naProd'], 'Inputs', '3', 'Position', [680 230 710 300]);
add_block('simulink/Math Operations/Gain', [proto '/GNa_g'], 'Gain', 'GNa', 'Position', [740 245 770 275]);
add_line(proto, 'ENa_c/1', 'naDrive/1', 'autorouting', 'on');
add_line(proto, 'mInf/1', 'naProd/1', 'autorouting', 'on');
add_line(proto, 'hint/1', 'naProd/2', 'autorouting', 'on');
add_line(proto, 'naDrive/1', 'naProd/3', 'autorouting', 'on');
add_line(proto, 'naProd/1', 'GNa_g/1', 'autorouting', 'on');
% membrane
add_block('simulink/Math Operations/Sum', [proto '/dV'], 'Inputs', '++++', 'Position', [820 150 850 210]);
add_block('simulink/Math Operations/Gain', [proto '/membrane'], 'Gain', '1000/Cm', 'Position', [890 155 920 185]);
add_block('simulink/Continuous/Integrator', [proto '/Vint'], 'InitialCondition', 'Vrest', 'Position', [950 151 980 189]);
add_block('simulink/Sinks/Out1', [proto '/V_mV'], 'Port', '1', 'Position', [1010 155 1040 169]);
add_line(proto, 'Iapp/1', 'IappSum/1', 'autorouting', 'on');
add_line(proto, 'IappSum/1', 'dV/1', 'autorouting', 'on');
add_line(proto, 'Isyn/1', 'dV/2', 'autorouting', 'on');
add_line(proto, 'Gm/1', 'dV/3', 'autorouting', 'on');
add_line(proto, 'GNa_g/1', 'dV/4', 'autorouting', 'on');
add_line(proto, 'dV/1', 'membrane/1', 'autorouting', 'on');
add_line(proto, 'membrane/1', 'Vint/1', 'autorouting', 'on');
add_line(proto, 'Vint/1', 'V_mV/1', 'autorouting', 'on');
add_line(proto, 'Vint/1', 'gTimesV/2', 'autorouting', 'on');
add_line(proto, 'Vint/1', 'dVm/1', 'autorouting', 'on');
add_line(proto, 'Vint/1', 'mArg/2', 'autorouting', 'on');
add_line(proto, 'Vint/1', 'hArg/2', 'autorouting', 'on');
add_line(proto, 'Vint/1', 'naDrive/2', 'autorouting', 'on');
m = Simulink.Mask.create(proto);
m.Description = ['W2L persistent-Na half-center (sns_toolbox ' ...
    'NonSpikingNeuronWithPersistentSodiumChannel semantics, FIXED tau_h). ' ...
    'Cm dV/dt = Gm(Vrest-V) + Iapp + sum(g*Esyn) - V*sum(g) + GNa*m*h*(ENa-V); ' ...
    'm = 1/(1+kM*exp(sM*(eM-V))); dh/dt = (hInf-h)/tauH with tauH FIXED ' ...
    '(SNS_NumpyFixedTau / AnimatLab semantics). Params = spinal/params.py ' ...
    'NAP with tau_max_h 0.25 s (build_w2l_net.py), membrane tau 50 ms.'];
NAPP = {'Vrest','Resting potential Vrest (mV)','0'; ...
        'Gm','Membrane conductance Gm (uS)','1'; ...
        'Cm','Membrane capacitance Cm (nF)','50'; ...
        'GNa','Persistent-Na conductance GNa (uS)','12'; ...
        'ENa','Na reversal ENa (mV)','8'; ...
        'kM','m-gate k','1'; 'sM','m-gate slope','0.8'; 'eM','m-gate midpoint (mV)','2'; ...
        'kH','h-gate k','1'; 'sH','h-gate slope','-2'; 'eH','h-gate midpoint (mV)','3.5'; ...
        'tauH','h-gate time constant tauH (ms), FIXED','250'; ...
        'h0','initial h (h_inf(Vrest) for these params)','0.9990889'};
for k = 1:size(NAPP, 1)
    m.addParameter('Name', NAPP{k, 1}, 'Type', 'edit', ...
        'Prompt', NAPP{k, 2}, 'Value', NAPP{k, 3});
end

%% ---------------- 2) circuit cells ------------------------------------------
P = struct('pos', containers.Map('KeyType','char','ValueType','any'), ...
           'nxt', containers.Map('KeyType','char','ValueType','double'), ...
           'cnt', containers.Map('KeyType','char','ValueType','double'));
XC = 80; XIN = 340; XRG = 660; XPF = 980; XPFI = 1250; XMN = 1540;

sides = {'L', 60; 'R', 1060};
for s = 1:2
    sd = sides{s, 1}; y0 = sides{s, 2};
    % NaP half-centers
    P = nap(mdl, P, proto, sprintf('%s_RG_ext', sd), XRG, y0+80);
    P = nap(mdl, P, proto, sprintf('%s_RG_flx', sd), XRG, y0+320);
    % RG-layer INs (tau 50 ms): InE, InF, C1, V3
    P = neu(mdl, P, sprintf('%s_RG_ext_IN', sd), 50, XIN, y0+60);
    P = neu(mdl, P, sprintf('%s_RG_flx_IN', sd), 50, XIN, y0+260);
    P = neu(mdl, P, sprintf('%s_C1', sd), 50, XIN, y0+460);
    P = neu(mdl, P, sprintf('%s_V3', sd), 50, XIN, y0+660);
    % PF half-centers + cross-INs (tau 80 ms)
    P = neu(mdl, P, sprintf('%s_PF_ext', sd), 80, XPF, y0+80);
    P = neu(mdl, P, sprintf('%s_PF_flx', sd), 80, XPF, y0+320);
    P = neu(mdl, P, sprintf('%s_PF_ext_IN', sd), 80, XPFI, y0+60);
    P = neu(mdl, P, sprintf('%s_PF_flx_IN', sd), 80, XPFI, y0+260);
    % representative MNs (tau 30 ms)
    P = neu(mdl, P, sprintf('%s_MN_ext', sd), 30, XMN, y0+100);
    P = neu(mdl, P, sprintf('%s_MN_flx', sd), 30, XMN, y0+340);
    % heel contact neuron (tau 40 ms)
    P = neu(mdl, P, sprintf('%s_heel', sd), 40, XC, y0+600);
end
delete_block(proto);   % prototype served its purpose (4 copies made)

%% ---------------- 3) synapses ----------------------------------------------
% RG lamination (per side)
for s = 1:2
    sd = sides{s, 1};
    P = wire(mdl, P, sprintf('%s_RG_ext', sd), sprintf('%s_RG_ext_IN', sd), g_lam, EXC);
    P = wire(mdl, P, sprintf('%s_RG_ext_IN', sd), sprintf('%s_RG_flx', sd), g_lam, INH);
    P = wire(mdl, P, sprintf('%s_RG_flx', sd), sprintf('%s_RG_flx_IN', sd), g_lam, EXC);
    P = wire(mdl, P, sprintf('%s_RG_flx_IN', sd), sprintf('%s_RG_ext', sd), g_lam, INH);
    P = wire(mdl, P, sprintf('%s_RG_ext', sd), sprintf('%s_RG_flx', sd), g_dir, EXC);
    P = wire(mdl, P, sprintf('%s_RG_flx', sd), sprintf('%s_RG_ext', sd), g_dir, EXC);
    % commissural c1 + V3 (crossed)
    o = 'R'; if strcmp(sd, 'R'), o = 'L'; end
    P = wire(mdl, P, sprintf('%s_RG_flx', sd), sprintf('%s_C1', sd), g_c1pre, EXC);
    P = wire(mdl, P, sprintf('%s_C1', sd), sprintf('%s_RG_flx', o), g_c1, INH);
    P = wire(mdl, P, sprintf('%s_RG_ext', sd), sprintf('%s_V3', sd), g_v3pre, EXC);
    P = wire(mdl, P, sprintf('%s_V3', sd), sprintf('%s_RG_ext_IN', o), g_v3, EXC);
    % PF drive + PF cross lamination
    P = wire(mdl, P, sprintf('%s_RG_ext', sd), sprintf('%s_PF_ext', sd), g_pfdrive, EXC);
    P = wire(mdl, P, sprintf('%s_RG_flx', sd), sprintf('%s_PF_flx', sd), g_pfdrive, EXC);
    P = wire(mdl, P, sprintf('%s_PF_ext', sd), sprintf('%s_PF_ext_IN', sd), g_pfcross, EXC);
    P = wire(mdl, P, sprintf('%s_PF_ext_IN', sd), sprintf('%s_PF_flx', sd), g_pfcross, INH);
    P = wire(mdl, P, sprintf('%s_PF_flx', sd), sprintf('%s_PF_flx_IN', sd), g_pfcross, EXC);
    P = wire(mdl, P, sprintf('%s_PF_flx_IN', sd), sprintf('%s_PF_ext', sd), g_pfcross, INH);
    % MNs
    P = wire(mdl, P, sprintf('%s_PF_ext', sd), sprintf('%s_MN_ext', sd), g_pfmn, EXC);
    P = wire(mdl, P, sprintf('%s_PF_flx', sd), sprintf('%s_MN_flx', sd), g_pfmn, EXC);
    % heel contact -> ipsilateral extensor layers (RG/PF/MN ext)
    P = wire(mdl, P, sprintf('%s_heel', sd), sprintf('%s_RG_ext', sd), g_contact, EXC);
    P = wire(mdl, P, sprintf('%s_heel', sd), sprintf('%s_PF_ext', sd), g_contact, EXC);
    P = wire(mdl, P, sprintf('%s_heel', sd), sprintf('%s_MN_ext', sd), g_contact, EXC);
end
n_syn = sum(cell2mat(P.cnt.values));
fprintf(['SNS_W2L_CPG core: 4 NaP half-centers + 22 graded cells (incl. 2 ' ...
         'heel contact encoders) + %d synapses\n'], n_syn);

%% ---------------- 4) DRIVE inputs ------------------------------------------
for s = 1:2
    sd = sides{s, 1}; y0 = sides{s, 2};
    % tonic + kickoff -> RG_ext ; tonic -> RG_flx (+kick on L ext / R flx only)
    add_block('simulink/Sources/Constant', [mdl '/tonicE_' sd], 'Value', num2str(tonicE), ...
        'Position', [XRG-320 y0-40 XRG-280 y0-20]);
    add_block('simulink/Sources/Constant', [mdl '/tonicF_' sd], 'Value', num2str(tonicF), ...
        'Position', [XRG-320 y0+220 XRG-280 y0+240]);
    add_block('simulink/Math Operations/Sum', [mdl '/sumE_' sd], 'Inputs', '++', ...
        'Position', [XRG-240 y0-30 XRG-210 y0-10]);
    add_block('simulink/Math Operations/Sum', [mdl '/sumF_' sd], 'Inputs', '++', ...
        'Position', [XRG-240 y0+230 XRG-210 y0+250]);
    add_line(mdl, ['tonicE_' sd '/1'], ['sumE_' sd '/1'], 'autorouting', 'on');
    add_line(mdl, ['tonicF_' sd '/1'], ['sumF_' sd '/1'], 'autorouting', 'on');
    add_line(mdl, ['sumE_' sd '/1'], [sprintf('%s_RG_ext', sd) '/1'], 'autorouting', 'on');
    add_line(mdl, ['sumF_' sd '/1'], [sprintf('%s_RG_flx', sd) '/1'], 'autorouting', 'on');
    % heel pulse train (20 nA, 1 Hz, 30% width; R phase-delayed 0.5 s)
    pd = '0'; if strcmp(sd, 'R'), pd = '0.5'; end
    add_block('simulink/Sources/Pulse Generator', [mdl '/heel_' sd], ...
        'PulseType', 'Time based', 'Amplitude', num2str(heelI), ...
        'Period', num2str(heelPeriod), 'PulseWidth', num2str(heelWidth), ...
        'PhaseDelay', pd, 'Position', [XC-10 y0+520 XC+70 y0+550]);
    add_line(mdl, ['heel_' sd '/1'], [sprintf('%s_heel', sd) '/1'], 'autorouting', 'on');
end
% kickoff pair: 10 nA for 10 ms at t=0 (L RG ext + R RG flx)
kickPW = num2str(100*kickT/24, '%.10g');      % % of 24 s period
add_block('simulink/Sources/Pulse Generator', [mdl '/kick_L'], ...
    'PulseType', 'Time based', 'Amplitude', num2str(kickAmp), 'Period', '24', ...
    'PulseWidth', kickPW, 'PhaseDelay', '0', 'Position', [XRG-320 -30 XRG-240 0]);
add_line(mdl, 'kick_L/1', 'sumE_L/2', 'autorouting', 'on');
add_block('simulink/Sources/Pulse Generator', [mdl '/kick_R'], ...
    'PulseType', 'Time based', 'Amplitude', num2str(kickAmp), 'Period', '24', ...
    'PulseWidth', kickPW, 'PhaseDelay', '0', 'Position', [XRG-320 970 XRG-240 1000]);
add_line(mdl, 'kick_R/1', 'sumF_R/2', 'autorouting', 'on');

%% ---------------- 5) logging + annotation + save ---------------------------
logit = {'L_RG_ext', 'vL_RG_ext'; 'L_RG_flx', 'vL_RG_flx'; ...
         'R_RG_ext', 'vR_RG_ext'; 'R_RG_flx', 'vR_RG_flx'; ...
         'L_MN_ext', 'vL_MN_ext'; 'L_MN_flx', 'vL_MN_flx'; ...
         'R_MN_ext', 'vR_MN_ext'; 'R_MN_flx', 'vR_MN_flx'; ...
         'L_heel', 'vL_heel'; 'R_heel', 'vR_heel'};
for k = 1:size(logit, 1)
    ph = get_param([mdl '/' logit{k, 1}], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', 'DataLoggingName', logit{k, 2});
end
try
    a = Simulink.Annotation(mdl, ['REPRESENTATIVE CORE of the W2L contact-driven ' ...
        'bilateral CPG (NOT the full 95-neuron net). 4 NaP half-centers (fixed ' ...
        'tau_h 250 ms) + laminated mutual inhibition + c1/V3 commissurals + ' ...
        'PF E/F + 4 representative MNs + heel-contact drive. Gains = ' ...
        'build_w2l_net.py W2L_GAINS. Validation vs w2l_numpy_ref.json.']);
    a.Position = [80 1900; 80 1900];
catch
end
set_param(mdl, 'Description', ['SNS_W2L_CPG: representative core of the W2L ' ...
    'contact-driven bilateral CPG (build_w2l_net.py), SNS_Library blocks + ' ...
    'custom NapHC (persistent Na, fixed tau_h). Validated against the numpy ' ...
    'reference campaigns/20260930/w2l_numpy_ref.json.']);
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', 'FixedStep', '0.002');
save_system(mdl, slxDst);
fprintf('saved %s\n', slxDst);

%% ---------------- 6) validation run ----------------------------------------
set_param(mdl, 'StopTime', num2str(T_END), 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
g = out.sigs;
t   = squeeze(g.get('vL_RG_ext').Values.Time);
vLE = squeeze(g.get('vL_RG_ext').Values.Data);
vLF = squeeze(g.get('vL_RG_flx').Values.Data);
vRE = squeeze(g.get('vR_RG_ext').Values.Data);
vRF = squeeze(g.get('vR_RG_flx').Values.Data);
close_system(mdl, 0);

win = t >= SKIP;
lE = vLE(win); rE = vRE(win); tw = t(win);
sL = burst_starts(lE); sR = burst_starts(rE);
nL = numel(sL); nR = numel(sR);
if nL >= 2, perL = mean(diff(tw(sL)))*1; else, perL = NaN; end %#ok<NASGU>
if nR >= 2, perR = mean(diff(tw(sR)))*1; else, perR = NaN; end %#ok<NASGU>
r = corrcoef(lE(:), rE(:)); r = r(1, 2);
fprintf(['SIMSCAPE W2L CORE: window %.0f-%.0f s | L RG E max %.3f mV, ' ...
         'R RG E max %.3f mV\n'], SKIP, T_END, max(lE), max(rE));
fprintf('bursts: L=%d (period %.3f s)  R=%d (period %.3f s); antiphase r=%.3f\n', ...
        nL, perL, nR, perR, r);

% numpy reference gate (same protocol, full 95-neuron net)
refJ = jsondecode(fileread(fullfile(root, '..', '..', 'MuJoCo_SNS', 'spinal', ...
    'campaigns', '20260930', 'w2l_numpy_ref.json')));
fprintf('numpy reference: period_L %.3f s, period_R %.3f s, r %.3f (full net)\n', ...
        refJ.period_L_s, refJ.period_R_s, refJ.antiphase_r);
okL = ~isnan(perL) && abs(perL - refJ.period_L_s) / refJ.period_L_s < 0.20;
okR = ~isnan(perR) && abs(perR - refJ.period_R_s) / refJ.period_R_s < 0.20;
okA = r < -0.5;
okN = nL >= 4 && nR >= 4;
pass = okL && okR && okA && okN;
fprintf(['W2L_SIMSCAPE %s | period L %.3f vs %.3f s (%.1f%%) R %.3f vs %.3f s ' ...
         '(%.1f%%) | antiphase r %.3f (ref %.3f) | bursts %d/%d\n'], ...
        tf2(pass), perL, refJ.period_L_s, 100*abs(perL-refJ.period_L_s)/refJ.period_L_s, ...
        perR, refJ.period_R_s, 100*abs(perR-refJ.period_R_s)/refJ.period_R_s, ...
        r, refJ.antiphase_r, nL, nR);

save(fullfile(root, 'results', 'SNS_W2L_CPG_run_20260930.mat'), ...
     't', 'vLE', 'vLF', 'vRE', 'vRF', 'perL', 'perR', 'r', 'nL', 'nR', 'pass');
fprintf('saved results/SNS_W2L_CPG_run_20260930.mat\n');
if ~pass, error('W2L_SIMSCAPE validation FAILED'); end
end

% ---------------- helpers ----------------------------------------------------
function P = nap(mdl, P, proto, name, x, y)
add_block(proto, [mdl '/' name], 'Position', [x y x+110 y+90]);
P.pos(name) = [x y]; P.nxt(name) = 2; P.cnt(name) = 0;
end

function P = neu(mdl, P, name, Cm, x, y)
add_block('SNS_Library/NonSpikingNeuron', [mdl '/' name], ...
    'Vrest', '0', 'Gm', '1', 'Cm', num2str(Cm), 'Thr', '0', 'Slope', '5', ...
    'Position', [x y x+90 y+80]);
P.pos(name) = [x y]; P.nxt(name) = 2; P.cnt(name) = 0;
end

function P = wire(mdl, P, src, dst, g, esyn)
k = P.cnt(dst) + 1; P.cnt(dst) = k;
pd = P.pos(dst);
nm = sprintf('syn%d_%s_to_%s', k, src, dst);
add_block('SNS_Library/NonSpikingSynapse', [mdl '/' nm], ...
    'gmax', num2str(g, 12), 'Esyn', num2str(esyn, 12), ...
    'ThrPre', '0', 'SlopePre', '5', ...
    'Position', [pd(1)-210 pd(2)+40+42*(k-1) pd(1)-150 pd(2)+72+42*(k-1)]);
set_param([mdl '/' nm], 'ShowName', 'off');
add_line(mdl, [src '/1'], [nm '/1'], 'autorouting', 'on');
port = P.nxt(dst); P.nxt(dst) = port + 1;
add_line(mdl, [nm '/1'], sprintf('%s/%d', dst, port), 'autorouting', 'on');
end

function s = burst_starts(sig)
on = sig > 0.5 * max(sig);
s = find(on(2:end) & ~on(1:end-1)) + 1;
end

function s = tf2(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
