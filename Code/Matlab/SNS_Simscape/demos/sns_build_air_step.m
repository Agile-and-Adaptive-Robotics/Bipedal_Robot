%% sns_build_air_step.m - BPA biped AIR-STEPPING model (curriculum stages 1-2)
%
% Ben's order (2026-10-24, curriculum definition): stages 1 and 2 are AIR
% STEPPING. This is the Simulink (base, no Simscape license needed)
% implementation: two planar legs (hip + knee each), pelvis FIXED (the
% harness analog of infant-support / air-stepping preparations - Thelen
% 1986), no ground contact, gravity ON (legs swing loaded only by their
% own weight).
%
% Architecture (all per-target conductances, no family collapse):
%   - ONE shared Deng RG (the VERIFIED SNS_Deng_RG.slx, self-running after
%     a 10 nA kick, period 1.938 s) drives TWO PF layers (SNS_Deng_PF.slx
%     Model blocks). The RIGHT leg's RG inputs are SWAPPED (ext<->flx) -
%     antiphase legs from one shared oscillator (single-shared-CPG v1;
%     two commissurally coupled RGs is the documented next step).
%   - Per-leg MNs: hip flx/ext, knee flx/ext (8 MNs), each driven from its
%     PF cell by a synapse carrying the PAPER'S OWN per-joint conductance
%     (Deng 2019 ModelDB 264844 mine, 2026-10-25):
%       hip flx 3.6325 uS, hip ext 2.6563 uS, knee flx 1.5165 uS,
%       knee ext 4.9296 uS (all Esyn -10 mV, MN Vrest -100, Deng Fig 6B
%       activation map - same as the verified CPG demo).
%   - MN S -> MuscleActivation (tau 150 ms) -> P = 620*A kPa -> BPA_20mm
%     (real Festo surfaces) with L = L0 -/+ r*theta (flx shortens flexor)
%     -> torque r*F about each joint.
%   - Mechanics: planar 2-link pendulum per leg (thigh+shank, lumped), base
%     Simulink integrators, gravity torques, viscous damping:
%       I_h*hipdd = T_hip - b*hipd - m_t*g*l_c*sin(hip) - coupling from shank (approx)
%       I_k*kneedd = T_knee - b*kneed - m_s*g*l_cs*sin(hip+knee)
%     (decoupled-joint approximation - DOCUMENTED; a full 2-link Lagrangian
%     is the upgrade. Air stepping tolerates this: no stance loads.)
%
% Output: demos\SNS_AirStep.slx. Runner: dev\run_air_step_20261025.m
% (checks: both hips oscillate, L-vs-R hip antiphase, knee cycling).

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
load_system('SNS_Library');
load_system(fullfile(here, 'SNS_Deng_RG.slx'));
load_system(fullfile(here, 'SNS_Deng_PF.slx'));

mdl = 'SNS_AirStep';
if bdIsLoaded(mdl), close_system(mdl, 0); end
slxPath = fullfile(here, [mdl '.slx']);
if exist(slxPath, 'file'), delete(slxPath); end
new_system(mdl);
set_param(mdl, 'UnconnectedInputMsg', 'none', 'ScreenColor', 'white');

% ---------- parameters (frozen into PreLoadFcn) ----------
P = struct( ...
  'I_h', 0.25, 'b_h', 0.8, ...          % thigh inertia/damping about hip (kg m^2, N m s)
  'I_k', 0.06, 'b_k', 0.4, ...          % shank about knee
  'm_t', 5.0, 'l_c', 0.20, ...          % thigh mass (kg), COM dist (m)
  'm_s', 2.0, 'l_cs', 0.22, ...         % shank mass, COM dist from knee
  'g', 9.81, ...
  'L0', 0.185, 'r_h', 0.045, 'r_k', 0.035, ...  % BPA rest lengths, moment arms
  'tauAct', 150, 'Pmax', 620, 'Fmax', 1148.3, 'Rest', 0.20, 'Kmax', 0.165);
gv = @(x) num2str(x, '%.10g');
pl = '';
fns = fieldnames(P);
for k = 1:numel(fns)
    pl = [pl newline fns{k} ' = ' gv(P.(fns{k})) ';']; %#ok<AGROW>
end
set_param(mdl, 'PreLoadFcn', pl);
eval(pl);

% ---------- shared Deng RG + two PF layers ----------
add_block('simulink/Ports & Subsystems/Model', [mdl '/RG'], ...
    'ModelFile', 'SNS_Deng_RG.slx', 'Position', [60 200 160 300]);
add_block('simulink/Sources/Pulse Generator', [mdl '/kick'], ...
    'PulseType', 'Time based', 'Amplitude', '10', 'Period', '500', ...
    'PulseWidth', '0.002', 'PhaseDelay', '0.1', 'Position', [20 240 60 280]);
add_line(mdl, 'kick/1', 'RG/1', 'autorouting', 'on');
add_block('simulink/Ports & Subsystems/Model', [mdl '/PF_L'], ...
    'ModelFile', 'SNS_Deng_PF.slx', 'Position', [260 120 360 260]);
add_block('simulink/Ports & Subsystems/Model', [mdl '/PF_R'], ...
    'ModelFile', 'SNS_Deng_PF.slx', 'Position', [260 360 360 500]);
add_line(mdl, 'RG/1', 'PF_L/1', 'autorouting', 'on');   % V_ext -> L PF ext
add_line(mdl, 'RG/2', 'PF_L/2', 'autorouting', 'on');
add_line(mdl, 'RG/2', 'PF_R/1', 'autorouting', 'on');   % SWAPPED: flx drives R ext
add_line(mdl, 'RG/1', 'PF_R/2', 'autorouting', 'on');

% ---------- per-joint channel builder (MN + activation + BPA pair + joint) ----------
% PF outs: 1 hip_e, 2 hip_f, 3 ka_e, 4 ka_f. Paper per-joint PF->MN gains.
GAIN = struct('hip_f', 3.6325, 'hip_e', 2.6563, 'knee_f', 1.5165, 'knee_e', 4.9296);
build_joint(mdl, 'L', 'hip', 1, 2, GAIN, P, [520 60]);
build_joint(mdl, 'L', 'knee', 3, 4, GAIN, P, [520 300]);
build_joint(mdl, 'R', 'hip', 2, 1, GAIN, P, [1000 60]);   % swapped ext/flx outs
build_joint(mdl, 'R', 'knee', 4, 3, GAIN, P, [1000 300]);
% PF->MN wiring: ext cell -> ext MN, flx cell -> flx MN. The RIGHT leg's
% antiphase comes from the swapped RG INPUTS above; outputs stay in
% natural order (double-swapping cancels the phase inversion - caught by
% the +1.000 hip correlation in the first smoke run).
wire_pf_mn(mdl, 'L', 'hip',  'PF_L', 1, 2);
wire_pf_mn(mdl, 'L', 'knee', 'PF_L', 3, 4);
wire_pf_mn(mdl, 'R', 'hip',  'PF_R', 1, 2);
wire_pf_mn(mdl, 'R', 'knee', 'PF_R', 3, 4);

set_param(mdl, 'SolverType', 'Variable-step', 'Solver', 'ode45', ...
    'MaxStep', '5e-4', 'StopTime', '40', 'UnconnectedInputMsg', 'none');
save_system(mdl, slxPath);
fprintf('built demos/SNS_AirStep.slx (shared Deng RG, 2 PF layers, 8 MNs, 8 BPA_20mm)\n');

%% ---------------- local functions ----------------
function build_joint(mdl, side, jname, pe, pf, GAIN, P, pos) %#ok<INUSD>
% One joint: flx/ext MNs + activations + BPA pair + planar pendulum.
x = pos(1); y = pos(2);
nm = [side '_' jname];
b = @(n, blk, v, varargin) add_block(blk, [mdl '/' n], ...
    'Position', v, varargin{:});
% MNs (Deng Fig 6B convention: Vrest -100, Thr -60, Slope 10)
b([nm '_MN_f'], 'SNS_Library/NonSpikingNeuron', [x y x+90 y+70], ...
  'Vrest', '-100', 'Gm', '1', 'Cm', '5', 'Thr', '-60', 'Slope', '10');
b([nm '_MN_e'], 'SNS_Library/NonSpikingNeuron', [x y+130 x+90 y+200], ...
  'Vrest', '-100', 'Gm', '1', 'Cm', '5', 'Thr', '-60', 'Slope', '10');
% activations
b([nm '_Act_f'], 'SNS_Library/MuscleActivation', [x+130 y+10 x+190 y+60], ...
  'tauAct', 'tauAct', 'A0', '0');
b([nm '_Act_e'], 'SNS_Library/MuscleActivation', [x+130 y+140 x+190 y+190], ...
  'tauAct', 'tauAct', 'A0', '0');
% pressures
b([nm '_P_f'], 'simulink/Math Operations/Gain', [x+220 y+15 x+250 y+55], 'Gain', 'Pmax');
b([nm '_P_e'], 'simulink/Math Operations/Gain', [x+220 y+145 x+250 y+185], 'Gain', 'Pmax');
% BPA blocks (fixed rest length here; length feedback below)
b([nm '_BPA_f'], 'SNS_Library/BPA_20mm', [x+290 y+5 x+370 y+65], ...
  'Rest', 'Rest', 'Kmax', 'Kmax', 'TendonL', '0', 'FittingL', '0');
b([nm '_BPA_e'], 'SNS_Library/BPA_20mm', [x+290 y+135 x+370 y+195], ...
  'Rest', 'Rest', 'Kmax', 'Kmax', 'TendonL', '0', 'FittingL', '0');
% net torque: flexion positive
b([nm '_rF'], 'simulink/Math Operations/Gain', [x+410 y+60 x+440 y+100], 'Gain', 'r_h');
b([nm '_rE'], 'simulink/Math Operations/Gain', [x+410 y+120 x+440 y+160], 'Gain', '-r_h');
if strcmp(jname, 'knee')
    set_param([mdl '/' nm '_rF'], 'Gain', 'r_k');
    set_param([mdl '/' nm '_rE'], 'Gain', '-r_k');
end
b([nm '_Tsum'], 'simulink/Math Operations/Sum', [x+470 y+85 x+500 y+135], 'Inputs', '++');
% pendulum: I*thdd = T - b*thd - grav(  + knee uses hip+theta abs angle)
b([nm '_grav'], 'simulink/User-Defined Functions/Fcn', ...
    [x+470 y+170 x+560 y+210], 'Expression', '-m_t*g*l_c*sin(u)');
if strcmp(jname, 'knee')
    % decoupled-joint v1 (documented): shank gravity on knee angle alone
    set_param([mdl '/' nm '_grav'], 'Expression', '-m_s*g*l_cs*sin(u)');
end
b([nm '_Iinv'], 'simulink/Math Operations/Gain', [x+590 y+90 x+620 y+130], 'Gain', '1/I_h');
if strcmp(jname, 'knee'), set_param([mdl '/' nm '_Iinv'], 'Gain', '1/I_k'); end
b([nm '_damp'], 'simulink/Math Operations/Gain', [x+590 y+160 x+620 y+200], 'Gain', '-b_h');
if strcmp(jname, 'knee'), set_param([mdl '/' nm '_damp'], 'Gain', '-b_k'); end
b([nm '_sumT'], 'simulink/Math Operations/Sum', [x+650 y+100 x+680 y+170], 'Inputs', '+++');
b([nm '_thdd_i'], 'simulink/Continuous/Integrator', [x+710 y+105 x+740 y+145]);
b([nm '_thd_i'], 'simulink/Continuous/Integrator', [x+770 y+105 x+800 y+145]);
% length feedback: theta -> BPA lengths (flexor shortens in flexion)
b([nm '_Lf'], 'simulink/User-Defined Functions/Fcn', ...
    [x+290 y-30 x+370 y], 'Expression', 'L0 - r_h*u');
if strcmp(jname, 'knee')
    set_param([mdl '/' nm '_Lf'], 'Expression', 'L0 - r_k*u');
end
b([nm '_Le'], 'simulink/User-Defined Functions/Fcn', ...
    [x+290 y+205 x+370 y+235], 'Expression', 'L0 + r_h*u');
if strcmp(jname, 'knee')
    set_param([mdl '/' nm '_Le'], 'Expression', 'L0 + r_k*u');
end
l = @(a, b_) add_line(mdl, a, b_, 'autorouting', 'on');
l([nm '_MN_f/2'], [nm '_Act_f/1']);
l([nm '_MN_e/2'], [nm '_Act_e/1']);
l([nm '_Act_f/1'], [nm '_P_f/1']);
l([nm '_Act_e/1'], [nm '_P_e/1']);
l([nm '_P_f/1'], [nm '_BPA_f/1']);
l([nm '_P_e/1'], [nm '_BPA_e/1']);
l([nm '_BPA_f/1'], [nm '_rF/1']);
l([nm '_BPA_e/1'], [nm '_rE/1']);
l([nm '_rF/1'], [nm '_Tsum/1']);
l([nm '_rE/1'], [nm '_Tsum/2']);
l([nm '_Tsum/1'], [nm '_sumT/1']);
l([nm '_thd_i/1'], [nm '_damp/1']);
l([nm '_damp/1'], [nm '_sumT/2']);
l([nm '_sumT/1'], [nm '_Iinv/1']);
l([nm '_Iinv/1'], [nm '_thdd_i/1']);
l([nm '_thdd_i/1'], [nm '_thd_i/1']);
l([nm '_thd_i/1'], [nm '_Lf/1']);
l([nm '_thd_i/1'], [nm '_Le/1']);
l([nm '_Lf/1'], [nm '_BPA_f/2']);
l([nm '_Le/1'], [nm '_BPA_e/2']);
% logging
ph = get_param([mdl '/' nm '_thd_i'], 'PortHandles');
set_param(ph.Outport(1), 'DataLogging', 'on', ...
    'DataLoggingNameMode', 'Custom', 'DataLoggingName', ['theta_' nm]);
end

function wire_pf_mn(mdl, side, jname, pfm, oe, of)
% PF cell -> per-target synapse (PAPER per-joint conductance) -> MN syn1.
G = struct('hip_f', 3.6325, 'hip_e', 2.6563, 'knee_f', 1.5165, 'knee_e', 4.9296);
nm = [side '_' jname];
for half = {'e', 'f'}
    h = half{1};
    sb = [nm '_synpf_' h];
    add_block('SNS_Library/NonSpikingSynapse', [mdl '/' sb], ...
        'gmax', num2str(G.([jname '_' h])), 'Esyn', '-10', ...
        'ThrPre', '-60', 'SlopePre', '10', 'ShowName', 'off', ...
        'Position', [400 + 480 * strcmp(side, 'R') 200 + 60 * strcmp(h, 'f') ...
                     445 + 480 * strcmp(side, 'R') 234 + 60 * strcmp(h, 'f')]);
    if strcmp(h, 'e'), p = oe; else, p = of; end
    add_line(mdl, sprintf('%s/%d', pfm, p), [sb '/1'], 'autorouting', 'on');
    add_line(mdl, [sb '/1'], [nm '_MN_' h '/2'], 'autorouting', 'on');
end
end
