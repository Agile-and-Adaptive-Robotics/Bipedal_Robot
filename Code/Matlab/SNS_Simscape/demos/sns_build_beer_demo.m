%% sns_build_beer_demo.m — build BeerCupReflexDemo.slx
%
% Ben's request (2026-09-10): "do a demo with a simple arm holding a cup.
% Show that the reflex pathway can hold the cup level as it is filled with
% beer." 2026-09-22 rework (Ben): biceps origin corrected (shoulder ->
% forearm, animation side), TRICEPS added with reciprocal (antagonist) Ia
% inhibition, triceps activation plotted, and the model now STARTS IN
% EQUILIBRIUM (activation initialized at the hold value, descending drive
% sized for the EMPTY cup) so there is no startup settling spike — the old
% demo's 5.13/9.66 deg "max sag" was a t~0.15 s startup transient, not the
% pour (measured from the committed results mat).
%
% Setup: 1-DOF elbow, forearm horizontal, cup in hand.
%   theta = sag angle from level (deg->rad internally; 0 = cup level, + = sag)
%   I*thdd = T_biceps - T_triceps - b*thd + T_gravity(theta, t)
%   T_gravity(t) = (m_arm*Lac + m_cup(t)*Lcup)*g*cos(theta)   [beer pours in]
%   m_cup(t) = pourRate*t
%   biceps  length L_bi  = L_hold  + c_len*theta   (sag lengthens the flexor)
%   triceps length L_tri = L_hold_t - c_len_t*theta (sag shortens the extensor)
%   (c_len < r_arm models series/tendon compliance: a real musculotendon unit
%   is softer than its moment arm would suggest, which is what leaves real
%   work for the reflex.)
%
% Reflex pathways (SNS_Library blocks, real BPA_20mm actuators):
%   descending drive -> MN_biceps / MN_triceps (baseline, empty-cup equilibrium)
%   Ia(flex)  --Exc--> MN_biceps    (stretch reflex: sag stretches biceps)
%   Ib(flex)  --Inh--> MN_biceps    (autogenic inhibition)
%   Ia(ext)   --Inh--> MN_biceps    (reciprocal inhibition)
%   Ia(ext)   --Exc--> MN_triceps   (triceps stretch reflex)
%   Ib(ext)   --Inh--> MN_triceps   (autogenic inhibition)
%   Ia(flex)  --Inh--> MN_triceps   (ANTAGONIST inhibition — as beer pours and
%                the biceps stretches, its Ia afferent INHIBITS the triceps MN,
%                so triceps activation DECREASES while biceps activation rises)
%
% kReflex scales ALL six reflex synapses. The runner runs kReflex = 1 / 0
% to show sag with the reflex pathway ON vs OFF.

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder
mdl = 'BeerCupReflexDemo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '10', 'ScreenColor', 'white', ...
    'UnconnectedInputMsg', 'none');

%% ---- parameters + startup equilibrium (computed, then frozen into the model) ----
p = struct();
% elbow model
p.I_el = 0.06; p.b_el = 2.0; p.g = 9.81;
p.m_arm = 1.5; p.Lac = 0.15; p.Lcup = 0.30;
p.pourRate = 0.0625;              % kg/s (0.5 kg beer over 8 s)
% BPAs (Festo 20 mm surfaces, same as the BPA_20mm block)
p.Pmax = 620; p.tauAct = 150;     % kPa supply, ms activation
p.Rest = 0.20; p.Kmax = 0.165;
p.L_hold = 0.185; p.c_len = 0.02; p.r_bi = 0.04;    % biceps (flexor)
p.L_hold_t = 0.185; p.c_len_t = 0.005; p.r_tri = 0.010; % triceps (extensor)
p.a0 = 0.257852586017; p.a1 = 6.4766142989; p.a3 = 1.32087718059;
p.Fmax = 620*(1.4877*atan(0.0248*(p.Rest-0.0075)*620));   % maxBPAforce, 20 mm
% reflex
p.sagN = 1/0.3;                   % rad -> normalized stretch (~0.3 rad full span)
p.kReflex = 1;
p.A_tri0 = 0.40;                  % baseline triceps co-contraction

% BPA force helpers (identical to the block internals / festo4.m)
relOf = @(L) (p.Rest - L)/(p.Rest - p.Kmax);           % relative contraction
fOf = @(Pn, L) p.Fmax * max(0, p.a0*(exp(-p.a1*relOf(L)) - 1) + Pn*exp(-p.a3*relOf(L)^2));
pnFor = @(F, L) (F/p.Fmax - p.a0*(exp(-p.a1*relOf(L)) - 1)) / exp(-p.a3*relOf(L)^2);

% --- startup equilibrium (empty cup, theta = 0, thetad = 0) ---
Tg0 = p.m_arm*p.Lac*p.g;                       % gravity torque, empty cup
Ftri0 = fOf(p.A_tri0, p.L_hold_t);             % triceps baseline force
Fbi0 = (Tg0 + p.r_tri*Ftri0)/p.r_bi;           % biceps must hold both
p.A_bi0 = min(1, max(0, pnFor(Fbi0, p.L_hold)));       % -> activation (P = Pmax*A)
% descending drives: MN steady state V0 = Thr + A0*Slope with Gm = 0.5,
% Vrest = -52, Thr = -45, Slope = 1  ->  desc = Gm*(V0 - Vrest)
p.desc_bi = 0.5*(7 + p.A_bi0);
p.desc_tri = 0.5*(7 + p.A_tri0);
fprintf('beer equilibrium: Fmax %.1f N, Ftri0 %.1f N, Fbi0 %.1f N -> A_bi0 %.3f (P %.0f kPa), A_tri0 %.3f (P %.0f kPa)\n', ...
    p.Fmax, Ftri0, Fbi0, p.A_bi0, p.Pmax*p.A_bi0, p.A_tri0, p.Pmax*p.A_tri0);

gv = @(x) num2str(x, '%.10g');
paramCmd = strjoin({ ...
    '% --- elbow model ---' ...
    ['I_el  = ' gv(p.I_el) ';      % kg*m^2 forearm+cup inertia about elbow'] ...
    ['b_el  = ' gv(p.b_el) ';         % N*m*s/rad damping'] ...
    ['g     = ' gv(p.g) ';        % m/s^2'] ...
    ['m_arm = ' gv(p.m_arm) ';         % kg forearm mass'] ...
    ['Lac   = ' gv(p.Lac) ';        % m forearm COM distance from elbow'] ...
    ['Lcup  = ' gv(p.Lcup) ';        % m cup distance from elbow'] ...
    ['pourRate = ' gv(p.pourRate) ';   % kg/s beer pour (0.5 kg over 8 s)'] ...
    '% --- BPAs (real Festo 20 mm equations) ---' ...
    ['Pmax  = ' gv(p.Pmax) ';         % kPa supply pressure'] ...
    ['tauAct = ' gv(p.tauAct) ';      % ms activation'] ...
    ['L_hold = ' gv(p.L_hold) ';      % m biceps length at level hold'] ...
    ['c_len = ' gv(p.c_len) ';      % m biceps length change per rad sag (series compliance)'] ...
    ['r_bi = ' gv(p.r_bi) ';        % m biceps moment arm'] ...
    ['L_hold_t = ' gv(p.L_hold_t) ';   % m triceps length at level hold'] ...
    ['c_len_t = ' gv(p.c_len_t) ';    % m triceps length change per rad sag'] ...
    ['r_tri = ' gv(p.r_tri) ';       % m triceps moment arm'] ...
    ['A_bi0 = ' gv(p.A_bi0) ';       % biceps activation holding the EMPTY cup (startup equilibrium)'] ...
    ['A_tri0 = ' gv(p.A_tri0) ';      % triceps baseline co-contraction'] ...
    '% --- reflex ---' ...
    ['desc_bi  = ' gv(p.desc_bi) ';   % nA descending drive (empty-cup equilibrium)'] ...
    ['desc_tri = ' gv(p.desc_tri) ';  % nA descending drive'] ...
    ['sagN   = ' gv(p.sagN) ';      % normalize sag rad -> spindle input (~0.3 rad full span)'] ...
    ['kReflex = ' gv(p.kReflex) ';         % set 0 in runner to disable reflex synapses'] ...
    }, newline);
set_param(mdl, 'PreLoadFcn', paramCmd);
eval(paramCmd);

%% ---- elbow model subsystem: [T_bi, T_tri] -> [theta, thd, L_bi, L_tri, m_cup] ----
emdl = [mdl '/ElbowModel'];
add_block('simulink/Ports & Subsystems/Subsystem', emdl, 'Position', [940 300 1100 440]);
delete_line(emdl, 'In1/1', 'Out1/1');
delete_block([emdl '/In1']);
delete_block([emdl '/Out1']);
B = @(n, blk, v) add_block(blk, [emdl '/' n], 'Position', v);
B('T_bi',    'simulink/Sources/In1',                  [25 63 55 77]);
B('T_tri',   'simulink/Sources/In1',                  [25 108 55 122]);
B('negBi',   'simulink/Math Operations/Gain',         [95 58 125 88]);    % -T_biceps (flexor holds up)
B('netT',    'simulink/Math Operations/Sum',          [160 70 190 130]);  % Tg - T_bi + T_tri - b*thd
B('invI',    'simulink/Math Operations/Gain',         [230 85 260 115]);
B('thd_i',   'simulink/Continuous/Integrator',        [300 85 330 115]);
B('th_i',    'simulink/Continuous/Integrator',        [360 85 390 115]);
B('damp',    'simulink/Math Operations/Gain',         [300 150 330 180]);
B('mCupR',   'simulink/Sources/Ramp',                 [95 240 125 270]);  % pourRate*t
B('mgArm',   'simulink/Sources/Constant',             [95 285 165 315]);  % m_arm*Lac
B('cosT',    'simulink/Math Operations/Trigonometric Function', [300 240 330 270]);
B('LcupM',   'simulink/Math Operations/Gain',         [190 240 220 270]);
B('mSum',    'simulink/Math Operations/Sum',          [250 250 280 290]);
B('gMul',    'simulink/Math Operations/Product',      [350 250 380 280]);
B('g_c',     'simulink/Sources/Constant',             [300 300 330 330]);
B('Tg',      'simulink/Math Operations/Product',      [420 250 450 280]); % *cos(theta)
B('LbiSum',  'simulink/Math Operations/Sum',          [420 330 450 370]); % L_hold + c_len*theta
B('LtriSum', 'simulink/Math Operations/Sum',          [420 390 450 430]); % L_hold_t - c_len_t*theta
B('Lhold',   'simulink/Sources/Constant',             [340 345 380 375]);
B('LholdT',  'simulink/Sources/Constant',             [340 410 380 440]);
B('cLenG',   'simulink/Math Operations/Gain',         [340 300 370 330]); % will hold both c_len terms via sign
set_param([emdl '/netT'],   'Inputs', '+++-');
set_param([emdl '/mSum'],   'Inputs', '++');
set_param([emdl '/LbiSum'], 'Inputs', '++');
set_param([emdl '/LtriSum'], 'Inputs', '+-');
set_param([emdl '/mgArm'],  'Value', 'm_arm*Lac');
set_param([emdl '/g_c'],    'Value', 'g');
set_param([emdl '/mCupR'],  'Slope', 'pourRate');
set_param([emdl '/LcupM'],  'Gain', 'Lcup');
set_param([emdl '/cosT'],   'Operator', 'cos');
set_param([emdl '/invI'],   'Gain', '1/I_el');
set_param([emdl '/damp'],   'Gain', 'b_el');
set_param([emdl '/negBi'],  'Gain', '-1');
set_param([emdl '/Lhold'],  'Value', 'L_hold');
set_param([emdl '/LholdT'], 'Value', 'L_hold_t');
set_param([emdl '/th_i'],   'InitialCondition', '0');
set_param([emdl '/thd_i'],  'InitialCondition', '0');
set_param([emdl '/cLenG'],  'Gain', 'c_len');
set_param([emdl '/LtriSum'], 'Inputs', '+-');
% second gain for the triceps c_len_t (separate block; negated inside LtriSum)
add_block('simulink/Math Operations/Gain', [emdl '/cLenTG'], 'Gain', 'c_len_t', 'Position', [340 455 370 485]);
outs = {'O_th',1; 'O_thd',2; 'O_Lbi',3; 'O_Ltri',4; 'O_mCup',5};
for k = 1:size(outs, 1)
    add_block('simulink/Sinks/Out1', [emdl '/' outs{k,1}], 'Port', num2str(outs{k,2}), ...
        'Position', [540 33+46*(k-1) 570 47+46*(k-1)]);
end
el = @(a, b) add_line(emdl, a, b, 'autorouting', 'on');
el('T_bi/1', 'negBi/1');
el('T_tri/1', 'netT/2');
el('negBi/1', 'netT/1');
el('netT/1', 'invI/1');
el('invI/1', 'thd_i/1');
el('thd_i/1', 'th_i/1');
el('thd_i/1', 'damp/1');
el('damp/1', 'netT/4');       % -b*thd (netT port 4 carries the minus sign)
el('th_i/1', 'cosT/1');
el('mCupR/1', 'LcupM/1');
el('LcupM/1', 'mSum/1');
el('mgArm/1', 'mSum/2');
el('mSum/1', 'gMul/1');
el('g_c/1', 'gMul/2');
el('gMul/1', 'Tg/1');
el('cosT/1', 'Tg/2');
el('Tg/1', 'netT/3');         % +T_gravity (sag direction)
el('th_i/1', 'O_th/1');
el('thd_i/1', 'O_thd/1');
el('th_i/1', 'cLenG/1');
el('Lhold/1', 'LbiSum/1');
el('cLenG/1', 'LbiSum/2');
el('LbiSum/1', 'O_Lbi/1');
el('th_i/1', 'cLenTG/1');
el('LholdT/1', 'LtriSum/1');
el('cLenTG/1', 'LtriSum/2');  % L_hold_t - c_len_t*theta
el('LtriSum/1', 'O_Ltri/1');
el('mCupR/1', 'O_mCup/1');
m = Simulink.Mask.create(emdl);
m.Type = 'SNS Elbow Model (1-DOF)';
m.Description = ['1-DOF elbow with cup: I*thdd = T_biceps - T_triceps - b*thd + T_grav, ' ...
    'T_grav = (m_cup(t)*Lcup + m_arm*Lac)*g*cos(theta); beer pours at pourRate. ' ...
    'Outputs: theta (sag, rad), thd, biceps length, triceps length, m_cup.'];
m.Display = elbowIconCode();
set_param(emdl, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
    'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');

%% ---- afferents ----
snsInst = @(pth, n, pos, varargin) add_block(pth, [mdl '/' n], 'Position', pos, varargin{:});
snsInst('SNS_Library/IaMuscleSpindle', 'Ia_biceps', [60 40 130 100], ...
    'Imax', '10', 'Wl', '24', 'Wv', '8');
snsInst('SNS_Library/IbGolgiTendon', 'Ib_biceps', [60 160 130 220], ...
    'Imax', '10', 'Kf', '0.008');
snsInst('SNS_Library/IaMuscleSpindle', 'Ia_triceps', [60 300 130 360], ...
    'Imax', '10', 'Wl', '24', 'Wv', '8');
snsInst('SNS_Library/IbGolgiTendon', 'Ib_triceps', [60 420 130 480], ...
    'Imax', '10', 'Kf', '0.008');
% spindle inputs: biceps stretches with sag; triceps stretches when the arm rises
add_block('simulink/Math Operations/Gain', [mdl '/thScl'],  'Gain', 'sagN',  'Position', [760 40 790 70]);
add_block('simulink/Math Operations/Gain', [mdl '/thdScl'], 'Gain', 'sagN',  'Position', [760 110 790 140]);
add_block('simulink/Math Operations/Gain', [mdl '/thSclT'], 'Gain', '-sagN', 'Position', [760 300 790 330]);
add_block('simulink/Math Operations/Gain', [mdl '/thdSclT'], 'Gain', '-sagN', 'Position', [760 370 790 400]);
add_line(mdl, 'ElbowModel/1', 'thScl/1', 'autorouting', 'on');
add_line(mdl, 'ElbowModel/2', 'thdScl/1', 'autorouting', 'on');
add_line(mdl, 'thScl/1', 'Ia_biceps/1', 'autorouting', 'on');
add_line(mdl, 'thdScl/1', 'Ia_biceps/2', 'autorouting', 'on');
add_line(mdl, 'ElbowModel/1', 'thSclT/1', 'autorouting', 'on');
add_line(mdl, 'ElbowModel/2', 'thdSclT/1', 'autorouting', 'on');
add_line(mdl, 'thSclT/1', 'Ia_triceps/1', 'autorouting', 'on');
add_line(mdl, 'thdSclT/1', 'Ia_triceps/2', 'autorouting', 'on');

%% ---- sensory neurons ----
snsInst('SNS_Library/NonSpikingNeuron', 'SN_Ia_bi',  [230 40 310 120], ...
    'Vrest', '-52', 'Gm', '0.12', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst('SNS_Library/NonSpikingNeuron', 'SN_Ib_bi',  [230 160 310 240], ...
    'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst('SNS_Library/NonSpikingNeuron', 'SN_Ia_tri', [230 300 310 380], ...
    'Vrest', '-52', 'Gm', '0.12', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst('SNS_Library/NonSpikingNeuron', 'SN_Ib_tri', [230 420 310 500], ...
    'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
add_line(mdl, 'Ia_biceps/1', 'SN_Ia_bi/1', 'autorouting', 'on');
add_line(mdl, 'Ib_biceps/1', 'SN_Ib_bi/1', 'autorouting', 'on');
add_line(mdl, 'Ia_triceps/1', 'SN_Ia_tri/1', 'autorouting', 'on');
add_line(mdl, 'Ib_triceps/1', 'SN_Ib_tri/1', 'autorouting', 'on');

%% ---- motoneurons ----
snsInst('SNS_Library/NonSpikingNeuron', 'MN_biceps', [640 40 740 200], ...
    'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
snsInst('SNS_Library/NonSpikingNeuron', 'MN_triceps', [640 300 740 460], ...
    'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
add_block('simulink/Sources/Constant', [mdl '/desc_bi_c'], 'Value', 'desc_bi', 'Position', [560 -10 600 14]);
add_block('simulink/Sources/Constant', [mdl '/desc_tri_c'], 'Value', 'desc_tri', 'Position', [560 250 600 274]);
add_line(mdl, 'desc_bi_c/1', 'MN_biceps/1', 'autorouting', 'on');
add_line(mdl, 'desc_tri_c/1', 'MN_triceps/1', 'autorouting', 'on');

%% ---- synapses: small, against the MN they synapse onto ----
% MN_biceps syn1..syn3
snsInst('SNS_Library/NonSpikingSynapse', 'syn_IaBi_exc',     [560 58 600 90], ...
    'gmax', '0.0005*kReflex',  'Esyn', '0',   'ThrPre', '-45', 'SlopePre', '0.5');
snsInst('SNS_Library/NonSpikingSynapse', 'syn_IbBi_inh',     [560 100 600 132], ...
    'gmax', '0.0006*kReflex', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst('SNS_Library/NonSpikingSynapse', 'syn_IaTri_inh_onBi', [560 142 600 174], ...
    'gmax', '0.0008*kReflex', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% MN_triceps syn1..syn3
snsInst('SNS_Library/NonSpikingSynapse', 'syn_IaTri_exc',    [560 318 600 350], ...
    'gmax', '0.001*kReflex',  'Esyn', '0',   'ThrPre', '-45', 'SlopePre', '0.5');
snsInst('SNS_Library/NonSpikingSynapse', 'syn_IbTri_inh',    [560 360 600 392], ...
    'gmax', '0.0004*kReflex', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst('SNS_Library/NonSpikingSynapse', 'syn_IaBi_inh_onTri', [560 402 600 434], ...
    'gmax', '0.0025*kReflex',  'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');
% presynaptic wiring
add_line(mdl, 'SN_Ia_bi/1', 'syn_IaBi_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib_bi/1', 'syn_IbBi_inh/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_tri/1', 'syn_IaTri_inh_onBi/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_tri/1', 'syn_IaTri_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib_tri/1', 'syn_IbTri_inh/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia_bi/1', 'syn_IaBi_inh_onTri/1', 'autorouting', 'on');
% synapse outputs -> postsynaptic syn ports
add_line(mdl, 'syn_IaBi_exc/1', 'MN_biceps/2', 'autorouting', 'on');
add_line(mdl, 'syn_IbBi_inh/1', 'MN_biceps/3', 'autorouting', 'on');
add_line(mdl, 'syn_IaTri_inh_onBi/1', 'MN_biceps/4', 'autorouting', 'on');
add_line(mdl, 'syn_IaTri_exc/1', 'MN_triceps/2', 'autorouting', 'on');
add_line(mdl, 'syn_IbTri_inh/1', 'MN_triceps/3', 'autorouting', 'on');
add_line(mdl, 'syn_IaBi_inh_onTri/1', 'MN_triceps/4', 'autorouting', 'on');
for nm = {'syn_IaBi_exc', 'syn_IbBi_inh', 'syn_IaTri_inh_onBi', 'syn_IaTri_exc', 'syn_IbTri_inh', 'syn_IaBi_inh_onTri'}
    try, set_param([mdl '/' nm{1}], 'ShowName', 'off'); catch, end
end

%% ---- muscles: MN S -> pentagon activation (A0 = equilibrium) -> BPA_20mm -> torque ----
snsInst('SNS_Library/MuscleActivation', 'Act_biceps', [800 60 860 120], ...
    'tauAct', 'tauAct', 'A0', 'A_bi0');
snsInst('SNS_Library/MuscleActivation', 'Act_triceps', [800 320 860 380], ...
    'tauAct', 'tauAct', 'A0', 'A_tri0');
snsInst('SNS_Library/BPA_20mm', 'BPA_biceps', [900 60 980 120], ...
    'Rest', '0.20', 'Kmax', '0.165', 'TendonL', '0', 'FittingL', '0');
snsInst('SNS_Library/BPA_20mm', 'BPA_triceps', [900 320 980 380], ...
    'Rest', '0.20', 'Kmax', '0.165', 'TendonL', '0', 'FittingL', '0');
add_block('simulink/Math Operations/Gain', [mdl '/P2kPa_bi'],  'Gain', 'Pmax', 'Position', [870 65 900 95]);
add_block('simulink/Math Operations/Gain', [mdl '/P2kPa_tri'], 'Gain', 'Pmax', 'Position', [870 325 900 355]);
add_block('simulink/Math Operations/Gain', [mdl '/rBi'],  'Gain', 'r_bi',  'Position', [1010 70 1040 100]);
add_block('simulink/Math Operations/Gain', [mdl '/rTri'], 'Gain', 'r_tri', 'Position', [1010 330 1040 360]);
add_line(mdl, 'MN_biceps/2', 'Act_biceps/1', 'autorouting', 'on');
add_line(mdl, 'MN_triceps/2', 'Act_triceps/1', 'autorouting', 'on');
add_line(mdl, 'Act_biceps/1', 'P2kPa_bi/1', 'autorouting', 'on');
add_line(mdl, 'Act_triceps/1', 'P2kPa_tri/1', 'autorouting', 'on');
add_line(mdl, 'P2kPa_bi/1', 'BPA_biceps/1', 'autorouting', 'on');
add_line(mdl, 'P2kPa_tri/1', 'BPA_triceps/1', 'autorouting', 'on');
add_line(mdl, 'ElbowModel/3', 'BPA_biceps/2', 'autorouting', 'on');
add_line(mdl, 'ElbowModel/4', 'BPA_triceps/2', 'autorouting', 'on');
add_line(mdl, 'BPA_biceps/1', 'rBi/1', 'autorouting', 'on');
add_line(mdl, 'BPA_triceps/1', 'rTri/1', 'autorouting', 'on');
add_line(mdl, 'rBi/1', 'ElbowModel/1', 'autorouting', 'on');
add_line(mdl, 'rTri/1', 'ElbowModel/2', 'autorouting', 'on');
% Ib sensing (force)
add_line(mdl, 'BPA_biceps/1', 'Ib_biceps/1', 'autorouting', 'on');
add_line(mdl, 'BPA_triceps/1', 'Ib_triceps/1', 'autorouting', 'on');

%% ---- logging ----
logNames = {'th', 'mCup', 'F_bi', 'A_bi', 'V_MN_bi', 'P_bi', 'F_tri', 'A_tri', 'V_MN_tri'};
for k = 1:numel(logNames)
    ypos = 560 + 45*(k-1);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logNames{k}], ...
        'VariableName', ['log_' logNames{k}], 'SaveFormat', 'Timeseries', ...
        'Position', [1180 ypos 1250 ypos+30]);
end
add_line(mdl, 'ElbowModel/1', 'log_th/1', 'autorouting', 'on');
add_line(mdl, 'ElbowModel/5', 'log_mCup/1', 'autorouting', 'on');
add_line(mdl, 'BPA_biceps/1', 'log_F_bi/1', 'autorouting', 'on');
add_line(mdl, 'Act_biceps/1', 'log_A_bi/1', 'autorouting', 'on');
add_line(mdl, 'MN_biceps/1', 'log_V_MN_bi/1', 'autorouting', 'on');
add_line(mdl, 'P2kPa_bi/1', 'log_P_bi/1', 'autorouting', 'on');
add_line(mdl, 'BPA_triceps/1', 'log_F_tri/1', 'autorouting', 'on');
add_line(mdl, 'Act_triceps/1', 'log_A_tri/1', 'autorouting', 'on');
add_line(mdl, 'MN_triceps/1', 'log_V_MN_tri/1', 'autorouting', 'on');

try
    anno = Simulink.Annotation(mdl, ...
        'Beer-cup reflex demo: biceps + triceps with reciprocal Ia inhibition hold the cup level while beer pours in (BPA_20mm actuators)');
    anno.Position = [40 -90 1100 -50];
catch
end

save_system(mdl);
fprintf('BeerCupReflexDemo.slx built (biceps + triceps, equilibrium start).\n');

%% ---------------- local functions ----------------
function s = elbowIconCode()
    % Heavy-outline box with elbow + cup glyph.
    s = strjoin({ ...
        'patch([-1 1 1 -1], [-1 -1 1 1], [0 0 0]);' ...
        'patch([-0.92 0.92 0.92 -0.92], [-0.84 -0.84 0.84 0.84], [0.93 0.93 0.93]);' ...
        'color(''black'');' ...
        'plot([-0.55 -0.45], [0.45 0.45]);' ...
        'plot([-0.45 -0.45], [0.45 -0.10]);' ...
        'plot([-0.45 -0.05], [-0.10 -0.05]);' ...
        'plot([-0.05 -0.05], [-0.05 -0.35]);' ...
        'patch(-0.05 + 0.08*cos(0.3927 + (0:0.0628:1.5708)), -0.05 + 0.08*sin(0.3927 + (0:0.0628:1.5708)), [1 1 1]);' ...
        'patch(-0.05 + 0.05*cos(0.3927 + (0:0.0628:1.5708)), -0.05 + 0.05*sin(0.3927 + (0:0.0628:1.5708)), [0 0 0]);' ...
        'plot([-0.15 0.05 0.05 -0.15 -0.15], [-0.35 -0.35 -0.55 -0.55 -0.35]);' ...
        'disp(''elbow'');' ...
        }, newline);
end
