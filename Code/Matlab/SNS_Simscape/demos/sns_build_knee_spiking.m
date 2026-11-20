%% sns_build_knee_spiking.m — build KneeReflexDemo_Spiking.slx (GOAL 2, 2026-10-02)
%
% SPIKING twin of KneeReflexDemo: identical knee model, afferents (Ia/Ib),
% motoneurons, muscles, and logging — ONLY the neural encoding of the
% reflex pathways changes (the hybrid doctrine from sns_toolbox Tutorials
% 3/5/9: SPIKING interneurons above, ANALOG voltages into the plant):
%
%   baseline:  afferent current -> NON-SPIKING sensory neuron (graded V)
%              -> NonSpikingSynapse (graded g)            -> MN syn port
%   spiking:   afferent current -> SPIKING interneuron (rate-coded spikes)
%              -> HybridSpikingSynapse (g jumps per spike) -> MN syn port
%
% The MNs stay NON-SPIKING (graded S(V) drive to the muscles) — spikes live
% in the sensory/interneuron layer only. Build = copy of the committed
% KneeReflexDemo.slx with the SN_* + syn_* blocks surgically replaced; the
% demo on disk is never touched.
%
% Gain matching (first-order, then hand-tuned): at the baseline settle
% point (theta ~ 43.5 deg) the afferent currents are Ia_ext ~2.9 nA,
% Ia_flex ~3.1 nA, Ib_ext ~5.2 nA, Ib_flex ~5.3 nA -> IN firing rates
% ~83/92/200/204 Hz; the hybrid-synapse ginc values are set so the MEAN
% synaptic conductance lands the MNs at the same operating point
% (V_MN ~ -44.5 mV -> A ~ 0.53) with the descending drives unchanged
% (desc_ext 4.2 / desc_flex 2.5 nA, same as baseline).
%
% SPIKING PARAMETERS:
%   IN (SpikingNeuron): Vrest -52, Vth -46 (6 mV above rest -> Ia fires
%     ~80-90 Hz, Ib ~200 Hz at baseline), Vreset -58, Gm 0.4 uS, Cm 2 nF
%     (tau 5 ms, same membrane speed as the baseline sensory neurons).
%   HybridSpikingSynapse: ThrPre -47 (comfortably BELOW Vth -46: a spiking
%     cell resets at threshold, so its V never renders >= Vth), tau_syn
%     10 ms, gmax 0.3 (cap above all operating conductances).

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));
load_system('SNS_Library');

src = 'KneeReflexDemo';
mdl = 'KneeReflexDemo_Spiking';
if bdIsLoaded(mdl), close_system(mdl, 0); end
f = fullfile(cdto, [mdl '.slx']);
if exist(f, 'file'), delete(f); end
load_system(fullfile(cdto, [src '.slx']));
save_system(src, f);          % in-memory copy, renamed (demo file untouched)

%% ---- remove the graded sensory layer (lines + SN_* + syn_*) ----
oldLines = { ...
    'Ia_ext/1',  'SN_Ia_ext/1';  'Ia_flex/1', 'SN_Ia_flex/1'; ...
    'Ib_ext/1',  'SN_Ib_ext/1';  'Ib_flex/1', 'SN_Ib_flex/1'; ...
    'SN_Ia_ext/1',  'syn_Iaext_exc/1';         'SN_Ib_ext/1',  'syn_Ibext_inh/1'; ...
    'SN_Ia_flex/1', 'syn_Iaflex_inh_on_ext/1'; 'SN_Ia_flex/1', 'syn_Iaflex_exc/1'; ...
    'SN_Ib_flex/1', 'syn_Ibflex_inh/1';        'SN_Ia_ext/1',  'syn_Iaext_inh_on_flex/1'; ...
    'syn_Iaext_exc/1',         'MN_ext/2'; 'syn_Ibext_inh/1', 'MN_ext/3'; ...
    'syn_Iaflex_inh_on_ext/1', 'MN_ext/4'; 'syn_Iaflex_exc/1', 'MN_flex/2'; ...
    'syn_Ibflex_inh/1',        'MN_flex/3'; 'syn_Iaext_inh_on_flex/1', 'MN_flex/4'};
for k = 1:size(oldLines, 1)
    try, delete_line(mdl, oldLines{k,1}, oldLines{k,2}); catch, end
end
for b = {'SN_Ia_ext','SN_Ia_flex','SN_Ib_ext','SN_Ib_flex', ...
         'syn_Iaext_exc','syn_Ibext_inh','syn_Iaflex_inh_on_ext', ...
         'syn_Iaflex_exc','syn_Ibflex_inh','syn_Iaext_inh_on_flex'}
    delete_block([mdl '/' b{1}]);
end

%% ---- spiking interneurons (afferent current -> spike rate) ----
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ia_ext'],  'Position', [180 40 270 130], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.4', 'Cm', '2');
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ia_flex'], 'Position', [180 190 270 280], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.4', 'Cm', '2');
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ib_ext'],  'Position', [180 340 270 430], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.4', 'Cm', '2');
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ib_flex'], 'Position', [180 440 270 530], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.4', 'Cm', '2');
add_line(mdl, 'Ia_ext/1',  'IN_Ia_ext/1',  'autorouting', 'on');
add_line(mdl, 'Ia_flex/1', 'IN_Ia_flex/1', 'autorouting', 'on');
add_line(mdl, 'Ib_ext/1',  'IN_Ib_ext/1',  'autorouting', 'on');
add_line(mdl, 'Ib_flex/1', 'IN_Ib_flex/1', 'autorouting', 'on');

%% ---- hybrid spiking synapses (spike events -> conductance -> MN syn port)
% GAIN MATCHING (final, hand-tuned in 3 iterations): match each pathway's
% mean spiking conductance gbar = ginc*f*tau_syn (observed rates Ia_ext
% ~100 Hz, Ia_flex ~95 Hz, tau_syn 10 ms) to the baseline graded synapse's
% effective conductance gmax*Sat(Vpre_SN) at the baseline settle point.
% tau_syn 10 ms: the graded baseline synapse has NO dynamics, so the
% shortest practical tau keeps the loop lag close (tau_syn 30 ms produced a
% 9-deg 4.7 Hz limit cycle; 10 ms ships a 3-5-dec ~7 Hz alternation).
% SHIPPED ginc [uS/spike]: exc_ext 0.087, exc_flex 0.126, Ib autogenic
% 0.078, reciprocal on ext 0.108, reciprocal on flex 0.072; gmax cap 0.3.
% RESULT vs baseline (2026-10-02): rise 15.79 deg (identical), settle mean
% 44.4 vs 43.5 deg, A_ext/A_flex 0.58/0.56 vs 0.53/0.53.
gIb  = 0.078;
HS = {'gmax', '0.3', 'tau_syn', '10', 'ThrPre', '-47'};
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_Iaext_exc'],         'Position', [580 58 620 90],  ...
    'ginc', '0.087', 'Esyn', '0',   HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_Ibext_inh'],         'Position', [580 100 620 132], ...
    'ginc', num2str(gIb, '%.6g'), 'Esyn', '-72', HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_Iaflex_inh_on_ext'], 'Position', [580 142 620 174], ...
    'ginc', '0.108', 'Esyn', '-72', HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_Iaflex_exc'],        'Position', [580 338 620 370], ...
    'ginc', '0.126', 'Esyn', '0',   HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_Ibflex_inh'],        'Position', [580 380 620 412], ...
    'ginc', num2str(gIb, '%.6g'), 'Esyn', '-72', HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_Iaext_inh_on_flex'], 'Position', [580 422 620 454], ...
    'ginc', '0.072', 'Esyn', '-72', HS{:});
add_line(mdl, 'IN_Ia_ext/1',  'syn_Iaext_exc/1',         'autorouting', 'on');
add_line(mdl, 'IN_Ib_ext/1',  'syn_Ibext_inh/1',         'autorouting', 'on');
add_line(mdl, 'IN_Ia_flex/1', 'syn_Iaflex_inh_on_ext/1', 'autorouting', 'on');
add_line(mdl, 'IN_Ia_flex/1', 'syn_Iaflex_exc/1',        'autorouting', 'on');
add_line(mdl, 'IN_Ib_flex/1', 'syn_Ibflex_inh/1',        'autorouting', 'on');
add_line(mdl, 'IN_Ia_ext/1',  'syn_Iaext_inh_on_flex/1', 'autorouting', 'on');
add_line(mdl, 'syn_Iaext_exc/1',         'MN_ext/2',  'autorouting', 'on');
add_line(mdl, 'syn_Ibext_inh/1',         'MN_ext/3',  'autorouting', 'on');
add_line(mdl, 'syn_Iaflex_inh_on_ext/1', 'MN_ext/4',  'autorouting', 'on');
add_line(mdl, 'syn_Iaflex_exc/1',        'MN_flex/2', 'autorouting', 'on');
add_line(mdl, 'syn_Ibflex_inh/1',        'MN_flex/3', 'autorouting', 'on');
add_line(mdl, 'syn_Iaext_inh_on_flex/1', 'MN_flex/4', 'autorouting', 'on');
for nm = {'syn_Iaext_exc', 'syn_Ibext_inh', 'syn_Iaflex_inh_on_ext', ...
          'syn_Iaflex_exc', 'syn_Ibflex_inh', 'syn_Iaext_inh_on_flex'}
    try, set_param([mdl '/' nm{1}], 'ShowName', 'off', 'NamePlacement', 'alternate'); catch, end
end

%% ---- spike-rate logging (IN spike lines) ----
add_block('simulink/Sinks/To Workspace', [mdl '/log_IN_Ia_ext_spk'], ...
    'VariableName', 'log_IN_Ia_ext_spk', 'SaveFormat', 'Timeseries', 'Position', [1160 610 1230 640]);
add_block('simulink/Sinks/To Workspace', [mdl '/log_IN_Ib_ext_spk'], ...
    'VariableName', 'log_IN_Ib_ext_spk', 'SaveFormat', 'Timeseries', 'Position', [1160 650 1230 680]);
add_line(mdl, 'IN_Ia_ext/2', 'log_IN_Ia_ext_spk/1', 'autorouting', 'on');
add_line(mdl, 'IN_Ib_ext/2', 'log_IN_Ib_ext_spk/1', 'autorouting', 'on');

%% ---- annotation ----
try
    anno = Simulink.Annotation(mdl, ...
        'SPIKING knee reflex: afferents -> SPIKING interneurons (rate coding) -> hybrid spiking synapses -> NON-SPIKING MNs -> BPAs (analog voltages into the plant)');
    anno.Position = [40 -75 1100 -35];
catch
end

save_system(mdl);
close_system(mdl, 0);
close_system(src, 0);
fprintf('%s.slx built (spiking interneuron layer on the KneeReflexDemo model).\n', mdl);
