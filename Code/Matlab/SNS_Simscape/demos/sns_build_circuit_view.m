%% sns_build_circuit_view.m — dissertation-ready Simulink VIEW of the
%% KneeReflexDemo neural circuit (2026-09-20).
%
% The full KneeReflexDemo.slx includes all plant arithmetic (torque sums,
% integrators, gains), which makes any print of it cluttered. This builder
% copies ONLY the neural + actuation blocks from the demo (mask values come
% along 1:1 — same gmax/Esyn/taus as the runnable model) and lays them out in
% the Rybak-style reading order of sns_draw_circuit.m:
%   afferent encoders -> sensory neurons -> synapses -> MNs -> activation -> BPA
% The Ib (force) feedback from the BPAs closes inside the view.
%
% Model is for PRINTING/export (sns_export_diagram), not simulation — ports
% are stubbed with In/Outport blocks.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                       % SNS_Simscape (SNS_Library)
cd(here);

src = 'KneeReflexDemo';
mdl = 'KneeReflexCircuit';
load_system('SNS_Library');
load_system(fullfile(fileparts(here), 'demos', [src '.slx']));
if bdIsLoaded(mdl), close_system(mdl, 0); end
f = [mdl '.slx'];
if exist(f, 'file'), delete(f); end
new_system(mdl); open_system(mdl);

% blocks to carry over from the demo, with fresh positions [x y x2 y2]
spec = { ...
 %  source block name          new name            position
 'Ia_ext',                    'Ia_ext',           [ 50  62 115 122];
 'Ia_flex',                   'Ia_flex',          [ 50 292 115 352];
 'Ib_ext',                    'Ib_ext',           [ 50 182 115 242];
 'Ib_flex',                   'Ib_flex',          [ 50 412 115 472];
 'SN_Ia_ext',                 'SN_Ia_ext',        [200  62 265 122];
 'SN_Ib_ext',                 'SN_Ib_ext',        [200 182 265 242];
 'SN_Ia_flex',                'SN_Ia_flex',       [200 292 265 352];
 'SN_Ib_flex',                'SN_Ib_flex',       [200 412 265 472];
 'desc_ext_c',                'desc_ext_c',       [520  -8 560   8];
 'desc_flex_c',               'desc_flex_c',      [520 222 560 238];
 'sumMN_ext',                 'sumMN_ext',        [520  70 545 100];
 'sumMN_flex',                'sumMN_flex',       [520 300 545 330];
 'MN_ext',                    'MN_ext',           [610  45 675 105];
 'MN_flex',                   'MN_flex',          [610 275 675 335];
 'syn_Iaext_exc',             'syn_Iaext_exc',    [360  60 410 110];
 'syn_Ibext_inh',             'syn_Ibext_inh',    [360 135 410 185];
 'syn_Iaflex_inh_on_ext',     'syn_Iaflex_inh_on_ext', [360 210 410 260];
 'syn_Iaflex_exc',            'syn_Iaflex_exc',   [360 310 410 360];
 'syn_Ibflex_inh',            'syn_Ibflex_inh',   [360 385 410 435];
 'syn_Iaext_inh_on_flex',     'syn_Iaext_inh_on_flex', [360 460 410 510];
 'Act_ext',                   'Act_ext',          [740  62 805 122];
 'Act_flex',                  'Act_flex',         [740 292 805 352];
 'BPA_ext',                   'BPA_ext',          [880  62 945 122];
 'BPA_flex',                  'BPA_flex',         [880 292 945 352];
 };
for k = 1:size(spec, 1)
    add_block([src '/' spec{k,1}], [mdl '/' spec{k,2}], 'Position', spec{k,3}, ...
        'MakeNameUnique', 'on');
end

% boundary ports (view stubs)
ports = { ...
 'simulink/Sources/In1', 'u_stretch_ext',  [ -80  72  -45  88],  '1', 'th\_norm (ext)';
 'simulink/Sources/In1', 'u_vel_ext',      [ -80 102  -45 118],  '2', 'thd\_norm';
 'simulink/Sources/In1', 'u_stretch_flex', [ -80 302  -45 318],  '3', 'stretch (flex)';
 'simulink/Sources/In1', 'u_vel_flex',     [ -80 332  -45 348],  '4', 'vel (flex)';
 'simulink/Sources/In1', 'u_strain_ext',   [ 800 172  835 188],  '5', 'strain\_ext';
 'simulink/Sources/In1', 'u_strain_flex',  [ 800 402  835 418],  '6', 'strain\_flex';
 'simulink/Sinks/Out1',  'F_ext',          [1010  72 1045  88],  '1', 'F\_ext';
 'simulink/Sinks/Out1',  'F_flex',         [1010 302 1045 318],  '2', 'F\_flex';
 };
for k = 1:size(ports, 1)
    add_block(ports{k,1}, [mdl '/' ports{k,2}], 'Position', ports{k,3}, ...
        'Port', ports{k,4});
end

% wiring — identical topology to the demo's neural section
L = { ...
 'u_stretch_ext/1',  'Ia_ext/1';    'u_vel_ext/1',      'Ia_ext/2';
 'u_stretch_flex/1', 'Ia_flex/1';   'u_vel_flex/1',     'Ia_flex/2';
 'Ia_ext/1',         'SN_Ia_ext/1'; 'Ia_flex/1',        'SN_Ia_flex/1';
 'BPA_ext/1',        'Ib_ext/1';    'BPA_flex/1',       'Ib_flex/1';
 'Ib_ext/1',         'SN_Ib_ext/1'; 'Ib_flex/1',        'SN_Ib_flex/1';
 'desc_ext_c/1',     'sumMN_ext/1'; 'desc_flex_c/1',    'sumMN_flex/1';
 'sumMN_ext/1',      'MN_ext/1';    'sumMN_flex/1',     'MN_flex/1';
 % presynaptic: sensory neurons -> synapses
 'SN_Ia_ext/1',      'syn_Iaext_exc/1';        'SN_Ib_ext/1', 'syn_Ibext_inh/1';
 'SN_Ia_flex/1',     'syn_Iaflex_inh_on_ext/1';'SN_Ia_flex/1','syn_Iaflex_exc/1';
 'SN_Ib_flex/1',     'syn_Ibflex_inh/1';       'SN_Ia_ext/1', 'syn_Iaext_inh_on_flex/1';
 % postsynaptic voltage: MN membrane -> Vpost
 'MN_ext/1',         'syn_Iaext_exc/2';        'MN_ext/1',    'syn_Ibext_inh/2';
 'MN_ext/1',         'syn_Iaflex_inh_on_ext/2';
 'MN_flex/1',        'syn_Iaflex_exc/2';       'MN_flex/1',   'syn_Ibflex_inh/2';
 'MN_flex/1',        'syn_Iaext_inh_on_flex/2';
 % synapse currents -> MN summing junctions
 'syn_Iaext_exc/1',      'sumMN_ext/2';  'syn_Ibext_inh/1',        'sumMN_ext/3';
 'syn_Iaflex_inh_on_ext/1','sumMN_ext/4'; 'syn_Iaflex_exc/1',      'sumMN_flex/2';
 'syn_Ibflex_inh/1',     'sumMN_flex/3'; 'syn_Iaext_inh_on_flex/1','sumMN_flex/4';
 % motoneuron drive -> activation -> BPA
 'MN_ext/2',         'Act_ext/1';    'MN_flex/2',         'Act_flex/1';
 'Act_ext/1',        'BPA_ext/1';    'Act_flex/1',        'BPA_flex/1';
 'u_strain_ext/1',   'BPA_ext/2';    'u_strain_flex/1',   'BPA_flex/2';
 'BPA_ext/1',        'F_ext/1';      'BPA_flex/1',        'F_flex/1';
 };
for k = 1:size(L, 1)
    add_line(mdl, L{k,1}, L{k,2}, 'autorouting', 'on');
end

% presentation: 12 pt everywhere (Ben: 10-12 pt), white canvas; synapse
% names hidden (export also enforces this, set here so the GUI view matches)
set_param(mdl, 'ScreenColor', 'white');
blks = find_system(mdl, 'Type', 'Block');
for b = 1:numel(blks)
    try, set_param(blks{b}, 'FontSize', '12'); catch, end
end
cand = find_system(mdl, 'Type', 'Block');
syns = {};
for b = 1:numel(cand)
    mt = '';
    try, mt = get_param(cand{b}, 'MaskType'); catch, end
    if strcmp(mt, 'SNS NonSpiking Synapse'), syns{end+1} = cand{b}; end %#ok<SAGROW>
end
for b = 1:numel(syns)
    set_param(syns{b}, 'ShowName', 'off');
end

save_system(mdl);
close_system(mdl, 0);
close_system(src, 0);
fprintf('%s.slx built (dissertation circuit view of %s)\n', mdl, src);
