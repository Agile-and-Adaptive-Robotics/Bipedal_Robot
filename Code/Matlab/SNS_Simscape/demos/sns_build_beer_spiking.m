%% sns_build_beer_spiking.m — build BeerCupReflexDemo_Spiking.slx (GOAL 2, 2026-10-02)
%
% SPIKING twin of BeerCupReflexDemo: identical elbow model, pour schedule,
% afferents, motoneurons, BPA_20mm actuators, and logging — only the
% sensory-neuron + synapse layer becomes SPIKING (hybrid doctrine: spikes
% above, analog voltages into the plant; MNs stay non-spiking).
%
% EQUILIBRIUM START carries over exactly: at t=0 (empty cup, level arm)
% every afferent current is ~0 nA -> the baseline SNs sit at Vrest -52 so
% their synapses are OFF (Sat<0); the spiking INs are BELOW THRESHOLD
% (silent) so the hybrid synapses are off too -> the MNs are driven by the
% same desc_bi/desc_tri descending drives at the same A0 activations.
%
% Gain matching: at the baseline's operating sag (~2-8 deg) the active
% baseline synapses are the two Ia_biceps paths (graded, SATURATED:
% SN_Ia_bi Vpre > ThrPre once Ia_biceps current > 0.84 nA, i.e. any sag
% > 0.004 rad), so their effective g = gmax*kReflex. The spiking match:
% gbar = ginc*f*tau_syn with Ia_bi firing f ~ 114 Hz at 2 deg sag
% (Gm 0.12 -> Vss -28.7 mV, tau_m 16.7 ms, Vth -46), tau_syn 10 ms
% -> ginc = gmax*1.14 (kReflex-scaled exactly like the baseline gmax).
% The four other pathways stay quiescent during the pour (Ib currents
% < 1 nA < threshold; Ia_triceps input clipped at 0 for positive sag) —
% same as baseline, their ginc values are carried for completeness.
%
% IN params: Ia INs keep the baseline SN Gm 0.12 (Cm 2 -> tau 16.7 ms),
% Ib INs keep Gm 0.4 (Cm 2); Vrest -52, Vth -46 (6 mV above rest),
% Vreset -58. Hybrid synapses: ThrPre -47, tau_syn 10 ms, gmax cap 0.01.

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(fileparts(cdto));
load_system('SNS_Library');

src = 'BeerCupReflexDemo';
mdl = 'BeerCupReflexDemo_Spiking';
if bdIsLoaded(mdl), close_system(mdl, 0); end
f = fullfile(cdto, [mdl '.slx']);
if exist(f, 'file'), delete(f); end
load_system(fullfile(cdto, [src '.slx']));
save_system(src, f);          % in-memory copy, renamed (demo file untouched)

%% ---- remove the graded sensory layer ----
oldLines = { ...
    'Ia_biceps/1', 'SN_Ia_bi/1';  'Ib_biceps/1', 'SN_Ib_bi/1'; ...
    'Ia_triceps/1', 'SN_Ia_tri/1'; 'Ib_triceps/1', 'SN_Ib_tri/1'; ...
    'SN_Ia_bi/1',  'syn_IaBi_exc/1';       'SN_Ib_bi/1',  'syn_IbBi_inh/1'; ...
    'SN_Ia_tri/1', 'syn_IaTri_inh_onBi/1'; 'SN_Ia_tri/1', 'syn_IaTri_exc/1'; ...
    'SN_Ib_tri/1', 'syn_IbTri_inh/1';      'SN_Ia_bi/1',  'syn_IaBi_inh_onTri/1'; ...
    'syn_IaBi_exc/1',       'MN_biceps/2'; 'syn_IbBi_inh/1',       'MN_biceps/3'; ...
    'syn_IaTri_inh_onBi/1', 'MN_biceps/4'; 'syn_IaTri_exc/1',      'MN_triceps/2'; ...
    'syn_IbTri_inh/1',      'MN_triceps/3'; 'syn_IaBi_inh_onTri/1', 'MN_triceps/4'};
for k = 1:size(oldLines, 1)
    try, delete_line(mdl, oldLines{k,1}, oldLines{k,2}); catch, end
end
for b = {'SN_Ia_bi', 'SN_Ib_bi', 'SN_Ia_tri', 'SN_Ib_tri', ...
         'syn_IaBi_exc', 'syn_IbBi_inh', 'syn_IaTri_inh_onBi', ...
         'syn_IaTri_exc', 'syn_IbTri_inh', 'syn_IaBi_inh_onTri'}
    delete_block([mdl '/' b{1}]);
end

%% ---- spiking interneurons (baseline SN Gm values kept) ----
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ia_bi'],  'Position', [230 40 320 130], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.12', 'Cm', '2');
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ib_bi'],  'Position', [230 160 320 250], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.4', 'Cm', '2');
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ia_tri'], 'Position', [230 300 320 390], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.12', 'Cm', '2');
add_block('SNS_Library/SpikingNeuron', [mdl '/IN_Ib_tri'], 'Position', [230 420 320 510], ...
    'Vrest', '-52', 'Vth', '-46', 'Vreset', '-58', 'Gm', '0.4', 'Cm', '2');
add_line(mdl, 'Ia_biceps/1',  'IN_Ia_bi/1',  'autorouting', 'on');
add_line(mdl, 'Ib_biceps/1',  'IN_Ib_bi/1',  'autorouting', 'on');
add_line(mdl, 'Ia_triceps/1', 'IN_Ia_tri/1', 'autorouting', 'on');
add_line(mdl, 'Ib_triceps/1', 'IN_Ib_tri/1', 'autorouting', 'on');

%% ---- hybrid spiking synapses (ginc = 1.14 x baseline gmax, kReflex-scaled)
HS = {'gmax', '0.01', 'tau_syn', '10', 'ThrPre', '-47'};
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_IaBi_exc'],       'Position', [560 58 600 90], ...
    'ginc', '0.00057*kReflex', 'Esyn', '0',   HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_IbBi_inh'],       'Position', [560 100 600 132], ...
    'ginc', '0.00068*kReflex', 'Esyn', '-72', HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_IaTri_inh_onBi'], 'Position', [560 142 600 174], ...
    'ginc', '0.00091*kReflex', 'Esyn', '-72', HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_IaTri_exc'],      'Position', [560 318 600 350], ...
    'ginc', '0.00114*kReflex', 'Esyn', '0',   HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_IbTri_inh'],      'Position', [560 360 600 392], ...
    'ginc', '0.00046*kReflex', 'Esyn', '-72', HS{:});
add_block('SNS_Library/HybridSpikingSynapse', [mdl '/syn_IaBi_inh_onTri'], 'Position', [560 402 600 434], ...
    'ginc', '0.00285*kReflex', 'Esyn', '-72', HS{:});
add_line(mdl, 'IN_Ia_bi/1',  'syn_IaBi_exc/1',       'autorouting', 'on');
add_line(mdl, 'IN_Ib_bi/1',  'syn_IbBi_inh/1',       'autorouting', 'on');
add_line(mdl, 'IN_Ia_tri/1', 'syn_IaTri_inh_onBi/1', 'autorouting', 'on');
add_line(mdl, 'IN_Ia_tri/1', 'syn_IaTri_exc/1',      'autorouting', 'on');
add_line(mdl, 'IN_Ib_tri/1', 'syn_IbTri_inh/1',      'autorouting', 'on');
add_line(mdl, 'IN_Ia_bi/1',  'syn_IaBi_inh_onTri/1', 'autorouting', 'on');
add_line(mdl, 'syn_IaBi_exc/1',       'MN_biceps/2',  'autorouting', 'on');
add_line(mdl, 'syn_IbBi_inh/1',       'MN_biceps/3',  'autorouting', 'on');
add_line(mdl, 'syn_IaTri_inh_onBi/1', 'MN_biceps/4',  'autorouting', 'on');
add_line(mdl, 'syn_IaTri_exc/1',      'MN_triceps/2', 'autorouting', 'on');
add_line(mdl, 'syn_IbTri_inh/1',      'MN_triceps/3', 'autorouting', 'on');
add_line(mdl, 'syn_IaBi_inh_onTri/1', 'MN_triceps/4', 'autorouting', 'on');
for nm = {'syn_IaBi_exc', 'syn_IbBi_inh', 'syn_IaTri_inh_onBi', ...
          'syn_IaTri_exc', 'syn_IbTri_inh', 'syn_IaBi_inh_onTri'}
    try, set_param([mdl '/' nm{1}], 'ShowName', 'off'); catch, end
end

%% ---- spike-rate logging ----
add_block('simulink/Sinks/To Workspace', [mdl '/log_IN_Ia_bi_spk'], ...
    'VariableName', 'log_IN_Ia_bi_spk', 'SaveFormat', 'Timeseries', 'Position', [1180 980 1250 1010]);
add_line(mdl, 'IN_Ia_bi/2', 'log_IN_Ia_bi_spk/1', 'autorouting', 'on');

%% ---- annotation ----
try
    anno = Simulink.Annotation(mdl, ...
        'SPIKING beer-cup reflex: spiking Ia/Ib interneurons (rate coding) -> hybrid synapses -> NON-SPIKING MNs -> BPA_20mm pair (starts in the same empty-cup equilibrium)');
    anno.Position = [40 -130 1100 -90];
catch
end

save_system(mdl);
close_system(mdl, 0);
close_system(src, 0);
fprintf('%s.slx built (spiking interneuron layer on the BeerCupReflexDemo model).\n', mdl);
