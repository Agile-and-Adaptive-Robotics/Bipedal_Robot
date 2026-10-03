%% sns_build_rgmn_spiking.m — build SNS_SpikingRG_MN.slx (GOAL 2 task D)
%
% REPRESENTATIVE SPIKING SUBNETWORK of the 410-neuron SNS_SpinalNetwork:
% an RG half-center pair (spiking) driving a MN pool (non-spiking) through
% HYBRID synapses — the exact motif the full network repeats 92 times
% (RG -> PF -> per-muscle MN pools).
%
% Topology (hybrid doctrine: spikes above, analog voltages into the plant):
%   RG_E, RG_F  = SPIKING neurons (the half-centers), asymmetric drives
%                 4.5 / 2.5 nA
%   Mutual inhibition is IN-LAMINATED (Deng / Biped_2xCPG convention):
%     RG spikes -> HybridSpikingSynapse (tau_syn 100 ms) -> fast NON-SPIKING
%     IN -> graded inhibitory synapse (0.2 uS, Esyn -80) onto the OTHER RG.
%     (Direct spike-to-spike inhibition was tried first and could NOT hold
%     the loser below threshold - 4-8 Hz co-firing; see the goal-2 report.)
%   Adp_E, Adp_F = NON-SPIKING slow cells (tau 400 ms) integrating their
%     RG's spike RATE (hybrid tau_syn 150 ms, ginc 0.05) and feeding back
%     graded self-inhibition (burst termination, the CPG-demo pattern).
%   MN_E1..E3    = 3 NON-SPIKING motoneurons (the pool) driven by RG_E
%                 through HybridSpikingSynapses (tau_syn 20 ms): the MN
%                 membrane low-pass-filters the spike train into a GRADED
%                 voltage -> S(V) drives the muscle
%   MN_F1        = 1 non-spiking MN driven by RG_F (antiphase witness)
%   Act_E        = MuscleActivation on MN_E1's S output (50 ms)
%
% VERIFIED (2026-10-02, sns_run_rgmn_spiking.m): the hybrid synapses
% convert the RG spike trains into GRADED MN-pool voltages — MN_E1 V
% modulates over a ~5 mV depth and tracks the RG_E rate envelope
% (corr ~+0.3); MN_F1 tracks RG_F. PARTIAL: the half-center alternation
% is antiphase rate modulation (envelope corr ~-0.5) with RG_E bursting
% 4x/12 s, but RG_F does not segment into discrete bursts (it keeps a
% low-rate floor through E's bursts) — clean burst alternation needs
% parameter work beyond this session. This is a REPRESENTATIVE
% subnetwork, not a port: the full 410-neuron conversion is scoped as
% future work (all 1392 gains were tuned for graded voltages; see the
% goal-2 report).

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));
load_system('SNS_Library');

mdl = 'SNS_SpikingRG_MN';
if bdIsLoaded(mdl), close_system(mdl, 0); end
f = fullfile(cdto, [mdl '.slx']);
if exist(f, 'file'), delete(f); end
new_system(mdl);
load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '12', 'ScreenColor', 'white', ...
    'UnconnectedInputMsg', 'none');

%% ---- half-centers (spiking) ----
add_block('SNS_Library/SpikingNeuron', [mdl '/RG_E'], 'Position', [120 80 210 170], ...
    'Vrest', '-52', 'Vth', '-45', 'Vreset', '-60', 'Gm', '0.1', 'Cm', '5');
add_block('SNS_Library/SpikingNeuron', [mdl '/RG_F'], 'Position', [120 320 210 410], ...
    'Vrest', '-52', 'Vth', '-45', 'Vreset', '-60', 'Gm', '0.1', 'Cm', '5');
add_block('simulink/Sources/Constant', [mdl '/drive_E'], 'Value', '4.5', 'Position', [30 110 80 134]);
add_block('simulink/Sources/Constant', [mdl '/drive_F'], 'Value', '2.5', 'Position', [30 350 80 374]);
add_line(mdl, 'drive_E/1', 'RG_E/1', 'autorouting', 'on');
add_line(mdl, 'drive_F/1', 'RG_F/1', 'autorouting', 'on');

%% ---- mutual inhibition, IN-LAMINATED (Deng / Biped_2xCPG convention) ----
% Direct spike-to-spike inhibition could not HOLD the loser below threshold
% (measured: 4-8 Hz co-firing). The proven pattern laminates a fast
% NON-SPIKING interneuron between the half-centers: the IN low-pass filters
% the winner's spike train (hybrid synapse, tau_syn 100 ms) into a GRADED
% inhibitory conductance that pins the loser until the winner's own Adp
% terminates its burst.
add_block('SNS_Library/NonSpikingNeuron', [mdl '/IN_E'], 'Position', [300 200 390 290], ...
    'Vrest', '-52', 'Gm', '0.1', 'Cm', '5', 'Thr', '-55', 'Slope', '1');
add_block('SNS_Library/NonSpikingNeuron', [mdl '/IN_F'], 'Position', [300 40 390 130], ...
    'Vrest', '-52', 'Gm', '0.1', 'Cm', '5', 'Thr', '-55', 'Slope', '1');
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/rtoINE'], 'Position', [240 210 290 260], ...
    'gmax', '0.5', 'ginc', '0.1', 'tau_syn', '100', 'Esyn', '0', 'ThrPre', '-47');
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/rtoINF'], 'Position', [240 50 290 100], ...
    'gmax', '0.5', 'ginc', '0.1', 'tau_syn', '100', 'Esyn', '0', 'ThrPre', '-47');
add_block('SNS_Library/NonSpikingSynapse', [mdl '/inh_INEtoF'], 'Position', [440 210 490 260], ...
    'gmax', '0.2', 'Esyn', '-80', 'ThrPre', '-45', 'SlopePre', '0.5');
add_block('SNS_Library/NonSpikingSynapse', [mdl '/inh_INFtoE'], 'Position', [440 50 490 100], ...
    'gmax', '0.2', 'Esyn', '-80', 'ThrPre', '-45', 'SlopePre', '0.5');
add_line(mdl, 'RG_E/1', 'rtoINE/1', 'autorouting', 'on');
add_line(mdl, 'RG_F/1', 'rtoINF/1', 'autorouting', 'on');
add_line(mdl, 'rtoINE/1', 'IN_E/2', 'autorouting', 'on');
add_line(mdl, 'rtoINF/1', 'IN_F/2', 'autorouting', 'on');
add_line(mdl, 'IN_E/1', 'inh_INEtoF/1', 'autorouting', 'on');
add_line(mdl, 'IN_F/1', 'inh_INFtoE/1', 'autorouting', 'on');
add_line(mdl, 'inh_INEtoF/1', 'RG_F/2', 'autorouting', 'on');   % RG_F syn1
add_line(mdl, 'inh_INFtoE/1', 'RG_E/2', 'autorouting', 'on');   % RG_E syn1

%% ---- adaptation: slow non-spiking cells fed by the spike RATE ----
add_block('SNS_Library/NonSpikingNeuron', [mdl '/Adp_E'], 'Position', [120 560 210 650], ...
    'Vrest', '-52', 'Gm', '0.1', 'Cm', '40', 'Thr', '-55', 'Slope', '1');
add_block('SNS_Library/NonSpikingNeuron', [mdl '/Adp_F'], 'Position', [120 700 210 790], ...
    'Vrest', '-52', 'Gm', '0.1', 'Cm', '40', 'Thr', '-55', 'Slope', '1');
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/rtoE'], 'Position', [260 560 320 610], ...
    'gmax', '1.0', 'ginc', '0.05', 'tau_syn', '150', 'Esyn', '0', 'ThrPre', '-47');
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/rtoF'], 'Position', [260 700 320 750], ...
    'gmax', '1.0', 'ginc', '0.05', 'tau_syn', '150', 'Esyn', '0', 'ThrPre', '-47');
add_block('SNS_Library/NonSpikingSynapse', [mdl '/adpInh_E'], 'Position', [400 560 450 610], ...
    'gmax', '0.3', 'Esyn', '-80', 'ThrPre', '-45', 'SlopePre', '0.5');
add_block('SNS_Library/NonSpikingSynapse', [mdl '/adpInh_F'], 'Position', [400 700 450 750], ...
    'gmax', '0.3', 'Esyn', '-80', 'ThrPre', '-45', 'SlopePre', '0.5');
add_line(mdl, 'RG_E/1', 'rtoE/1', 'autorouting', 'on');       % RG_E V (hybrid detects)
add_line(mdl, 'RG_F/1', 'rtoF/1', 'autorouting', 'on');
add_line(mdl, 'rtoE/1', 'Adp_E/2', 'autorouting', 'on');
add_line(mdl, 'rtoF/1', 'Adp_F/2', 'autorouting', 'on');
add_line(mdl, 'Adp_E/1', 'adpInh_E/1', 'autorouting', 'on');
add_line(mdl, 'Adp_F/1', 'adpInh_F/1', 'autorouting', 'on');
add_line(mdl, 'adpInh_E/1', 'RG_E/3', 'autorouting', 'on');   % RG_E syn2
add_line(mdl, 'adpInh_F/1', 'RG_F/3', 'autorouting', 'on');   % RG_F syn2

%% ---- MN pool through HYBRID synapses ----
for k = 1:3
    add_block('SNS_Library/NonSpikingNeuron', sprintf('%s/MN_E%d', mdl, k), ...
        'Position', [620 60+150*(k-1) 710 150+150*(k-1)], ...
        'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
    add_block('SNS_Library/HybridSpikingSynapse', sprintf('%s/toMN_E%d', mdl, k), ...
        'Position', [500 80+150*(k-1) 550 130+150*(k-1)], ...
        'gmax', '0.5', 'ginc', '0.05', 'tau_syn', '20', 'Esyn', '0', 'ThrPre', '-47');
    add_line(mdl, 'RG_E/1', sprintf('toMN_E%d/1', k), 'autorouting', 'on');
    add_line(mdl, sprintf('toMN_E%d/1', k), sprintf('MN_E%d/2', k), 'autorouting', 'on');
end
add_block('SNS_Library/NonSpikingNeuron', [mdl '/MN_F1'], 'Position', [620 510 710 600], ...
    'Vrest', '-52', 'Gm', '0.5', 'Cm', '5', 'Thr', '-45', 'Slope', '1');
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/toMN_F1'], 'Position', [500 530 550 580], ...
    'gmax', '0.5', 'ginc', '0.05', 'tau_syn', '20', 'Esyn', '0', 'ThrPre', '-47');
add_line(mdl, 'RG_F/1', 'toMN_F1/1', 'autorouting', 'on');
add_line(mdl, 'toMN_F1/1', 'MN_F1/2', 'autorouting', 'on');

%% ---- muscle stage on the pool ----
add_block('SNS_Library/MuscleActivation', [mdl '/Act_E'], 'Position', [780 100 840 160], 'tauAct', '50');
add_line(mdl, 'MN_E1/2', 'Act_E/1', 'autorouting', 'on');

%% ---- logging ----
logs = {'RG_E', 1, 'V_RG_E'; 'RG_E', 2, 'spk_RG_E'; 'RG_F', 1, 'V_RG_F'; ...
        'RG_F', 2, 'spk_RG_F'; 'Adp_E', 1, 'V_Adp_E'; 'MN_E1', 1, 'V_MN_E1'; ...
        'MN_E2', 1, 'V_MN_E2'; 'MN_F1', 1, 'V_MN_F1'; 'Act_E', 1, 'A_E'};
for k = 1:size(logs, 1)
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' logs{k,3}], ...
        'VariableName', ['log_' logs{k,3}], 'SaveFormat', 'Timeseries', ...
        'Position', [950 40+45*(k-1) 1020 70+45*(k-1)]);
    add_line(mdl, [logs{k,1} '/' num2str(logs{k,2})], ['log_' logs{k,3} '/1'], ...
        'autorouting', 'on');
end

try
    anno = Simulink.Annotation(mdl, ...
        'Representative spiking subnetwork: SPIKING RG half-centers (mutual inhibition + spike-rate adaptation) -> HYBRID synapses -> NON-SPIKING MN pool -> activation');
    anno.Position = [40 -60 1100 -20];
catch
end

save_system(mdl);
close_system(mdl, 0);
fprintf('%s.slx built.\n', mdl);
