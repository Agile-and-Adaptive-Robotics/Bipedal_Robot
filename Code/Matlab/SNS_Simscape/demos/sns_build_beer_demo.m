%% sns_build_beer_demo.m — build BeerCupReflexDemo.slx
%
% Ben's request (2026-09-10): "do a demo with a simple arm holding a cup.
% Show that the reflex pathway can hold the cup level as it is filled with
% beer."
%
% Setup: 1-DOF elbow, forearm horizontal, cup in hand.
%   theta = sag angle from level (deg->rad internally; 0 = cup level, + = sag)
%   I*thdd = T_biceps - T_gravity - b*thd
%   T_gravity(t) = (m_arm*Lac + m_cup(t)*Lcup)*g*cos(theta)   [beer pours in]
%   m_cup(t) = pourRate*t  (0 -> m_full over the pour)
%   biceps length L = L_hold + c_len*theta  (sag lengthens the biceps)
%
% Reflex pathway (SNS_Library blocks, real BPA_20mm actuator):
%   descending drive -> MN_biceps (baseline activation)
%   Ia spindle (stretch+vel of theta) --Exc--> MN   (stretch reflex: sag ->
%       biceps stretch -> more drive -> more pressure -> more force)
%   Ib GTO (biceps force)             --Inh--> MN   (autogenic inhibition)
%   MN S -> MuscleActivation -> pressure P = Pmax*A -> BPA_20mm(P, L) -> torque
%
% kIa/kIb scale the reflex synapses. Runner runs twice (kReflex = 1 / 0)
% to show sag with the reflex pathway ON vs OFF.

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder
mdl = 'BeerCupReflexDemo';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '10', 'ScreenColor', 'white');

%% ---- parameters ----
paramCmd = strjoin({ ...
    '% --- elbow plant ---' ...
    'I_el  = 0.06;        % kg*m^2 forearm+cup inertia about elbow' ...
    'b_el  = 2.0;         % N*m*s/rad damping (heavy: damps the hold transient)' ...
    'g     = 9.81;        % m/s^2' ...
    'm_arm = 1.5;         % kg forearm mass' ...
    'Lac   = 0.15;        % m forearm COM distance from elbow' ...
    'Lcup  = 0.30;        % m cup distance from elbow' ...
    'pourRate = 0.0625;   % kg/s beer pour (0.5 kg over 8 s)' ...
    '% --- biceps BPA (real Festo equations) ---' ...
    'Pmax  = 620;         % kPa supply pressure' ...
    'tauAct = 40;         % ms activation' ...
    'L_hold = 0.185;      % m biceps length at level hold' ...
    'c_len = 0.06;        % m biceps length change per rad of sag' ...
    'r_arm = 0.04;        % m biceps moment arm about elbow' ...
    '% --- reflex ---' ...
          'desc  = 3.685;        % nA descending drive: S=0.42 -> P=261 kPa -> holds full cup at theta=0' ...
    'sagN   = 1/0.3;      % normalize sag rad -> spindle input (~0.3 rad full span)' ...
    'kReflex = 1;         % set 0 in runner to disable reflex synapses' ...
    }, newline);
set_param(mdl, 'PreLoadFcn', paramCmd);
eval(paramCmd);

%% ---- elbow plant (masked subsystem): [T_biceps] -> [theta, thd, L_mus, m_cup, T_grav] ----
plant = [mdl '/ElbowPlant'];
add_block('simulink/Ports & Subsystems/Subsystem', plant, 'Position', [900 300 1060 440]);
delete_line(plant, 'In1/1', 'Out1/1');
delete_block([plant '/In1']);
delete_block([plant '/Out1']);
B = @(n, p, v) add_block(p, [plant '/' n], 'Position', v);
B('T_bi',    'simulink/Sources/In1',                  [25 63 55 77]);
B('netT',    'simulink/Math Operations/Sum',          [120 55 150 95]);   % T_bi - T_g - b*thd
B('invI',    'simulink/Math Operations/Gain',         [190 60 220 90]);
B('thd_i',   'simulink/Continuous/Integrator',        [250 60 280 90]);
B('th_i',    'simulink/Continuous/Integrator',        [310 60 340 90]);
B('damp',    'simulink/Math Operations/Gain',         [250 130 280 160]);
B('negD',    'simulink/Math Operations/Gain',         [180 170 210 200]); % sign flip into netT
B('mCupR',   'simulink/Sources/Ramp',                 [120 230 150 260]); % pourRate*t
B('mgArm',   'simulink/Sources/Constant',             [120 280 190 310]); % m_arm*Lac*g
B('cosT',    'simulink/Math Operations/Trigonometric Function', [310 130 340 160]);
B('LcupM',   'simulink/Math Operations/Gain',         [190 230 220 260]); % *Lcup
B('LarmM',   'simulink/Math Operations/Gain',         [190 280 220 310]); % *Lac? precomputed
B('mSum',    'simulink/Math Operations/Sum',          [250 245 280 285]); % mCup*Lcup + mArm*Lac
B('gMul',    'simulink/Math Operations/Product',      [310 245 340 275]);
B('g_c',     'simulink/Sources/Constant',             [280 290 310 320]); % g
B('negBi',   'simulink/Math Operations/Gain',         [430 190 460 220]); % -T_biceps
B('Tg',      'simulink/Math Operations/Product',      [370 150 400 180]); % * cos(theta)
B('Lsum',    'simulink/Math Operations/Sum',          [400 250 430 290]); % L_hold + c_len*theta
B('Lhold',   'simulink/Sources/Constant',             [310 300 370 330]);
B('cLenG',   'simulink/Math Operations/Gain',         [340 250 370 280]);
add_block('simulink/Sinks/Out1', [plant '/O_th'],   'Port', '1', 'Position', [520 48 550 62]);
add_block('simulink/Sinks/Out1', [plant '/O_thd'],  'Port', '2', 'Position', [520 118 550 132]);
add_block('simulink/Sinks/Out1', [plant '/O_L'],    'Port', '3', 'Position', [520 258 550 272]);
add_block('simulink/Sinks/Out1', [plant '/O_mCup'], 'Port', '4', 'Position', [520 328 550 342]);
set_param([plant '/netT'],   'Inputs', '+++');
set_param([plant '/mSum'],   'Inputs', '++');
set_param([plant '/Lsum'],   'Inputs', '++');
set_param([plant '/mgArm'],  'Value', 'm_arm*Lac');
set_param([plant '/g_c'],    'Value', 'g');
set_param([plant '/mCupR'],  'Slope', 'pourRate');
set_param([plant '/LcupM'],  'Gain', 'Lcup');
set_param([plant '/cosT'],   'Operator', 'cos');
set_param([plant '/invI'],   'Gain', '1/I_el');
set_param([plant '/damp'],   'Gain', 'b_el');
set_param([plant '/cLenG'],  'Gain', 'c_len');
set_param([plant '/Lhold'],  'Value', 'L_hold');
set_param([plant '/negD'],   'Gain', '-1');
set_param([plant '/negBi'],  'Gain', '-1');
set_param([plant '/th_i'],   'InitialCondition', '0');
pl = @(a, b) add_line(plant, a, b, 'autorouting', 'on');
% netT = T_grav(+sag) - b*thd - T_biceps: theta positive = cup sags down
pl('Tg/1', 'netT/1');
pl('negD/1', 'netT/2');
pl('negBi/1', 'netT/3');
pl('T_bi/1', 'negBi/1');
pl('netT/1', 'invI/1');
pl('invI/1', 'thd_i/1');
pl('thd_i/1', 'th_i/1');
pl('thd_i/1', 'damp/1');
pl('damp/1', 'negD/1');
pl('th_i/1', 'cosT/1');
pl('mCupR/1', 'LcupM/1');
pl('LcupM/1', 'mSum/1');
pl('mgArm/1', 'mSum/2');
pl('mSum/1', 'gMul/1');
pl('g_c/1', 'gMul/2');
pl('th_i/1', 'O_th/1');
pl('thd_i/1', 'O_thd/1');
pl('gMul/1', 'Tg/1');
pl('cosT/1', 'Tg/2');
pl('th_i/1', 'cLenG/1');
pl('Lhold/1', 'Lsum/1');
pl('cLenG/1', 'Lsum/2');
pl('Lsum/1', 'O_L/1');
pl('mCupR/1', 'O_mCup/1');
% netT ports: 1=T_biceps(+), 2=-(b*thd), 3=-(T_grav)

m = Simulink.Mask.create(plant);
m.Type = 'SNS Elbow Plant (1-DOF)';
m.Description = ['1-DOF elbow with cup: I*thdd = T_biceps - b*thd - T_grav, ' ...
    'T_grav = (m_cup(t)*Lcup + m_arm*Lac)*g*cos(theta); beer pours at pourRate. ' ...
    'Outputs: theta (sag, rad), thd, biceps length L, m_cup.'];
m.Display = plantIconCode();
set_param(plant, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
    'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');

%% ---- reflex circuit ----
snsInst = @(p, n, pos, varargin) add_block(p, [mdl '/' n], 'Position', pos, varargin{:});
snsInst('SNS_Library/NonSpikingNeuron', 'MN_biceps', [420 60 530 170], ...
    'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
snsInst('SNS_Library/IaMuscleSpindle', 'Ia_biceps', [120 60 220 140], ...
    'Imax', '10', 'Wl', '6', 'Wv', '8');
snsInst('SNS_Library/IbGolgiTendon', 'Ib_biceps', [120 200 220 280], ...
    'Imax', '10', 'Kf', '0.008');
snsInst('SNS_Library/NonSpikingNeuron', 'SN_Ia', [270 50 370 150], ...
    'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst('SNS_Library/NonSpikingNeuron', 'SN_Ib', [270 190 370 290], ...
    'Vrest', '-52', 'Gm', '0.4', 'Cm', '2', 'Thr', '-55', 'Slope', '1');
snsInst('SNS_Library/NonSpikingSynapse', 'syn_Ia_exc', [395 230 465 300], ...
    'gmax', '0.005', 'Esyn', '0', 'ThrPre', '-45', 'SlopePre', '0.5');
snsInst('SNS_Library/NonSpikingSynapse', 'syn_Ib_inh', [395 320 465 390], ...
    'gmax', '0.0015', 'Esyn', '-72', 'ThrPre', '-45', 'SlopePre', '0.5');

add_block('simulink/Math Operations/Sum', [mdl '/sumMN'], 'Inputs', '+++', 'Position', [380 90 410 140]);
add_block('simulink/Sources/Constant', [mdl '/desc_c'], 'Value', 'desc', 'Position', [300 90 340 120]);
add_line(mdl, 'desc_c/1', 'sumMN/1', 'autorouting', 'on');
add_line(mdl, 'sumMN/1', 'MN_biceps/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ia/1', 'syn_Ia_exc/1', 'autorouting', 'on');
add_line(mdl, 'SN_Ib/1', 'syn_Ib_inh/1', 'autorouting', 'on');
add_line(mdl, 'Ia_biceps/1', 'SN_Ia/1', 'autorouting', 'on');
add_line(mdl, 'Ib_biceps/1', 'SN_Ib/1', 'autorouting', 'on');
add_line(mdl, 'syn_Ia_exc/1', 'sumMN/2', 'autorouting', 'on');
add_line(mdl, 'syn_Ib_inh/1', 'sumMN/3', 'autorouting', 'on');
add_line(mdl, 'MN_biceps/1', 'syn_Ia_exc/2', 'autorouting', 'on');
add_line(mdl, 'MN_biceps/1', 'syn_Ib_inh/2', 'autorouting', 'on');

%% ---- muscle chain: MN S -> activation -> pressure -> BPA_20mm -> torque ----
snsInst('SNS_Library/MuscleActivation', 'Act_biceps', [570 70 650 120], 'tauAct', '150');
snsInst('SNS_Library/BPA_20mm', 'BPA_biceps', [730 60 830 160], ...
    'Rest', '0.20', 'Kmax', '0.165', 'TendonL', '0', 'FittingL', '0');
add_block('simulink/Math Operations/Gain', [mdl '/P2kPa'], 'Gain', 'Pmax', 'Position', [690 75 720 105]);
add_block('simulink/Math Operations/Gain', [mdl '/rBi'], 'Gain', 'r_arm', 'Position', [860 75 890 105]);
add_line(mdl, 'MN_biceps/2', 'Act_biceps/1', 'autorouting', 'on');
add_line(mdl, 'Act_biceps/1', 'P2kPa/1', 'autorouting', 'on');
add_line(mdl, 'P2kPa/1', 'BPA_biceps/1', 'autorouting', 'on');
add_line(mdl, 'ElbowPlant/3', 'BPA_biceps/2', 'autorouting', 'on');
add_line(mdl, 'BPA_biceps/1', 'rBi/1', 'autorouting', 'on');
add_line(mdl, 'rBi/1', 'ElbowPlant/1', 'autorouting', 'on');
% afferent sensing
% normalize theta/thd (rad, rad/s) before the spindle: sag span ~0.3 rad
add_block('simulink/Math Operations/Gain', [mdl '/thScl'], 'Gain', 'sagN', 'Position', [940 55 970 85]);
add_block('simulink/Math Operations/Gain', [mdl '/thdScl'], 'Gain', 'sagN', 'Position', [940 125 970 155]);
add_line(mdl, 'ElbowPlant/1', 'thScl/1', 'autorouting', 'on');
add_line(mdl, 'ElbowPlant/2', 'thdScl/1', 'autorouting', 'on');
add_line(mdl, 'thScl/1', 'Ia_biceps/1', 'autorouting', 'on');
add_line(mdl, 'thdScl/1', 'Ia_biceps/2', 'autorouting', 'on');
add_line(mdl, 'BPA_biceps/1', 'Ib_biceps/1', 'autorouting', 'on');

%% ---- logging ----
logNames = {'th', 'mCup', 'F_bi', 'A_bi', 'V_MN', 'P_bi'};
for k = 1:numel(logNames)
    ypos = 480 + 60*(k-1);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logNames{k}], ...
        'VariableName', ['log_' logNames{k}], 'SaveFormat', 'Timeseries', ...
        'Position', [1150 ypos 1220 ypos+30]);
end
add_line(mdl, 'ElbowPlant/1', 'log_th/1', 'autorouting', 'on');
add_line(mdl, 'ElbowPlant/4', 'log_mCup/1', 'autorouting', 'on');
add_line(mdl, 'BPA_biceps/1', 'log_F_bi/1', 'autorouting', 'on');
add_line(mdl, 'Act_biceps/1', 'log_A_bi/1', 'autorouting', 'on');
add_line(mdl, 'MN_biceps/1', 'log_V_MN/1', 'autorouting', 'on');
add_line(mdl, 'P2kPa/1', 'log_P_bi/1', 'autorouting', 'on');

try
    anno = Simulink.Annotation(mdl, ...
        'Beer-cup reflex demo: Ia stretch reflex + Ib autogenic inhibition hold the cup level while beer pours in (BPA_20mm biceps)');
    anno.Position = [40 -90 1100 -50];
catch
end

save_system(mdl);
fprintf('BeerCupReflexDemo.slx built.\n');

%% ---------------- local functions ----------------
function s = plantIconCode()
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
