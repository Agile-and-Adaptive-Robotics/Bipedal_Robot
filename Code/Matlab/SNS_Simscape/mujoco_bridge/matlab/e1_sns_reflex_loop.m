function e1_sns_reflex_loop()
% E1: FIRST CLOSED-LOOP NEUROMECHANICAL SIM on this machine:
% SNS_Library reflex arc <-> MuJoCo gait2392 (via the proven bridge).
%
% Circuit (vas_med_r, real knee): plant sensors (muscle len/vel/force)
%   -> Ia spindle encoder (SNS IaMuscleSpindle) + II length encoder (gain)
%      + Ib GTO encoder (SNS IbGolgiTendon)
%   -> Ia/II/Ib afferent neurons (taus from the tuned network JSON)
%   -> synapses onto MN_vas_med_r with the TUNED conductances
%      (ia_to_mn 0.6, ii_to_mn 0.4, ib_to_mn_inh 0.35 uS)
%   -> MN S(V) -> ctrl(29) = vas_med_r -> muscle force -> sensors again.
% Constant 2 nA MN bias sets the operating point (V ~ 2 mV, S ~ 0.4).
%
% Encoder mappings (documented choices):
%   stretch_norm = clip(0.5 + (len-L0)/(2*(Lmax-L0)), 0, 1)   % 0.5 at rest
%   vel_norm     = clip(0.5 + vel/0.4, 0, 1)                  % +-0.2 m/s fs
%   force_norm   = clip(-frc/Fmax, 0, 1)        % frc negative in pull
% Runs the same model twice: reflex OFF (gmax 0) vs ON, 1.5 s, ode45.
% Outputs: exp_out\e1_reflex.mat + e1_reflex.png.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..'));               % SNS_Simscape (SNS_Library)
xml = strrep(fullfile(here, '..', '..', '..', '..', '..', 'Solid_Models', ...
    'OpenSim', 'Gait2392_Robotbody', 'mjc', 'gait2392_simbody', ...
    'gait2392_simbody_cvt3_simbridge.xml'), '\', '/');
C = load(fullfile(here, 'exp_const.mat'));         % L0,Lmax,Fmax,taus,g's

mdl = 'sns_bridge_reflex';
if bdIsLoaded(mdl), close_system(mdl, 0); end
outDir = fullfile(here, 'exp_out');
if ~isfolder(outDir), mkdir(outDir); end
if exist(fullfile(outDir, [mdl '.slx']), 'file')
    delete(fullfile(outDir, [mdl '.slx']));
end
new_system(mdl);

%% ---------------- plant ---------------------------------------------------
add_block('mjLib/MuJoCo Plant', [mdl '/Plant'], ...
    'xmlFileRel', xml, 'renderingType', 'None', ...
    'rgbOutOption', 'off', 'depthOutOption', 'off', ...
    'Position', [900 100 1060 240]);

%% ---------------- afferent encoders (from plant sensor bus) ---------------
add_block('simulink/Signal Routing/Bus Selector', [mdl '/sens'], ...
    'OutputSignals', 'knee_r_pos,vas_med_r_len,vas_med_r_vel,vas_med_r_frc', ...
    'Position', [1140 100 1160 190]);
add_line(mdl, 'Plant/1', 'sens/1', 'autorouting', 'on');

% stretch_norm in [0,1]: 0.5 at L0, 1.0 at Lmax
add_block('simulink/Math Operations/Gain', [mdl '/dL'], ...
    'Gain', num2str(1 / (2 * (C.Lmax - C.L0))), 'Position', [1220 60 1260 100]);
add_block('simulink/Math Operations/Bias', [mdl '/half'], ...
    'Bias', '0.5', 'Position', [1290 60 1320 100]);
add_block('simulink/Discontinuities/Saturation', [mdl '/satL'], ...
    'UpperLimit', '1', 'LowerLimit', '0', 'Position', [1350 60 1380 100]);
add_line(mdl, 'sens/2', 'dL/1', 'autorouting', 'on');
add_line(mdl, 'dL/1', 'half/1', 'autorouting', 'on');
add_line(mdl, 'half/1', 'satL/1', 'autorouting', 'on');

% vel_norm: +-0.2 m/s -> 0..1
add_block('simulink/Math Operations/Gain', [mdl '/vN'], ...
    'Gain', '2.5', 'Position', [1220 130 1260 170]);
add_block('simulink/Math Operations/Bias', [mdl '/halfV'], ...
    'Bias', '0.5', 'Position', [1290 130 1320 170]);
add_block('simulink/Discontinuities/Saturation', [mdl '/satV'], ...
    'UpperLimit', '1', 'LowerLimit', '0', 'Position', [1350 130 1380 170]);
add_line(mdl, 'sens/3', 'vN/1', 'autorouting', 'on');
add_line(mdl, 'vN/1', 'halfV/1', 'autorouting', 'on');
add_line(mdl, 'halfV/1', 'satV/1', 'autorouting', 'on');

% forceNorm = clip(-frc/Fmax, 0, 1)
add_block('simulink/Math Operations/Gain', [mdl '/fN'], ...
    'Gain', num2str(-1 / C.Fmax), 'Position', [1220 200 1260 240]);
add_block('simulink/Discontinuities/Saturation', [mdl '/satF'], ...
    'UpperLimit', '1', 'LowerLimit', '0', 'Position', [1290 200 1320 240]);
add_line(mdl, 'sens/4', 'fN/1', 'autorouting', 'on');
add_line(mdl, 'fN/1', 'satF/1', 'autorouting', 'on');

add_block('SNS_Library/IaMuscleSpindle', [mdl '/IaEnc'], ...
    'Imax', '10', 'Wl', '4', 'Wv', '8', 'Position', [1440 80 1540 180]);
add_line(mdl, 'satL/1', 'IaEnc/1', 'autorouting', 'on');
add_line(mdl, 'satV/1', 'IaEnc/2', 'autorouting', 'on');
add_block('simulink/Math Operations/Gain', [mdl '/IIEnc'], ...
    'Gain', '4', 'Position', [1440 220 1480 260]);
add_line(mdl, 'satL/1', 'IIEnc/1', 'autorouting', 'on');
add_block('SNS_Library/IbGolgiTendon', [mdl '/IbEnc'], ...
    'Imax', '10', 'Kf', '10', 'Position', [1440 300 1540 380]);
add_line(mdl, 'satF/1', 'IbEnc/1', 'autorouting', 'on');

%% ---------------- afferent neurons + synapses (tuned network values) -----
tauof = @(nm) get_tau(here, nm);
add_neu(mdl, 'Ia_n', tauof('Ia_vas_med_r'), [1600 80 1680 160]);
add_neu(mdl, 'II_n', tauof('II_vas_med_r'), [1600 220 1680 300]);
add_neu(mdl, 'Ib_n', tauof('Ib_vas_med_r'), [1600 340 1680 420]);
add_line(mdl, 'IaEnc/1', 'Ia_n/1', 'autorouting', 'on');
add_line(mdl, 'IIEnc/1', 'II_n/1', 'autorouting', 'on');
add_line(mdl, 'IbEnc/1', 'Ib_n/1', 'autorouting', 'on');

add_syn(mdl, 'Ia_syn', C.g_ia, C.esyn_exc, [1760 100 1860 180]);
add_syn(mdl, 'II_syn', C.g_ii, C.esyn_exc, [1760 250 1860 330]);
add_syn(mdl, 'Ib_syn', C.g_ib, C.esyn_inh, [1760 400 1860 480]);
affpairs = {'Ia_n', 'Ia_syn'; 'II_n', 'II_syn'; 'Ib_n', 'Ib_syn'};
for j = 1:size(affpairs, 1)
    add_line(mdl, [affpairs{j, 1} '/1'], [affpairs{j, 2} '/1'], ...
             'autorouting', 'on');
end

%% ---------------- MN + bias + output assembly -----------------------------
add_block('simulink/Sources/Constant', [mdl '/Ibias'], 'Value', '2', ...
          'Position', [1760 550 1810 590]);
add_block('simulink/Math Operations/Sum', [mdl '/mnSum'], 'Inputs', '++++', ...
          'Position', [1930 300 1960 360]);
cells3 = {'Ia', 'II', 'Ib'};
for j = 1:3
    add_line(mdl, sprintf('%s_syn/1', cells3{j}), ...
             sprintf('mnSum/%d', j), 'autorouting', 'on');
end
add_line(mdl, 'Ibias/1', 'mnSum/4', 'autorouting', 'on');

add_neu(mdl, 'MN', C.tau_mn, [2020 280 2100 360]);
add_line(mdl, 'mnSum/1', 'MN/1', 'autorouting', 'on');
for j = 1:3   % postsynaptic voltage feedback to each synapse
    add_line(mdl, 'MN/1', sprintf('%s_syn/2', cells3{j}), ...
             'autorouting', 'on');
end

add_mlfn(mdl, 'ctrlvec', [ ...
    'function u = ctrlvec(s)' newline ...
    '%#codegen' newline ...
    'u = zeros(92,1); u(29) = s;'], [2160 300 2260 360]);
add_line(mdl, 'MN/2', 'ctrlvec/1', 'autorouting', 'on');
add_line(mdl, 'ctrlvec/1', 'Plant/1', 'autorouting', 'on');

%% ---------------- logging + sim (OFF then ON) -----------------------------
logs = { {'MN', 1, 'V_mn'}, {'MN', 2, 'S'}, {'Ia_n', 1, 'V_ia'}, ...
         {'Ib_n', 1, 'V_ib'}, {'sens', 1, 'knee'}, {'sens', 2, 'len'}, ...
         {'sens', 3, 'vel'}, {'sens', 4, 'frc'} };
for k = 1:numel(logs)
    ph = get_param([mdl '/' logs{k}{1}], 'PortHandles');
    set_param(ph.Outport(logs{k}{2}), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', logs{k}{3});
end

set_param(mdl, 'Solver', 'ode45', 'RelTol', '1e-5', 'AbsTol', '1e-7', ...
          'StopTime', '1.5', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');
gmap = struct('Ia_syn', 'g_ia', 'II_syn', 'g_ii', 'Ib_syn', 'g_ib');
R = struct();
for mode = ["OFF", "ON"]
    s = fieldnames(gmap);
    for k = 1:numel(s)
        if mode == "ON", g = C.(gmap.(s{k})); else, g = 0; end
        set_param([mdl '/' s{k}], 'gmax', num2str(g, 12));
    end
    out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
    for f = ["V_mn", "S", "V_ia", "V_ib", "knee", "len", "vel", "frc"]
        R.(mode).(f) = out.sigs.get(char(f)).Values;
    end
end
save_system(mdl, fullfile(outDir, [mdl '.slx']));
close_system(mdl, 0);

%% ---------------- report + figure ----------------------------------------
kd = @(r) rad2deg(r.knee.Data(:));
for mode = ["OFF", "ON"]
    r = R.(mode);
    kv = kd(r);
    fprintf(['E1 %s: knee %.1f -> %.1f deg (min %.1f), |frc| peak %.0f N, ' ...
             'S range [%.2f, %.2f], V_Ia [%.1f, %.1f] mV, V_Ib [%.1f, %.1f] mV\n'], ...
            mode, kv(1), kv(end), min(kv), ...
            max(abs(r.frc.Data(:))), min(r.S.Data(:)), max(r.S.Data(:)), ...
            min(r.V_ia.Data(:)), max(r.V_ia.Data(:)), ...
            min(r.V_ib.Data(:)), max(r.V_ib.Data(:)));
end

fig = figure('Visible', 'off', 'Position', [100 100 900 760]);
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact');
fields = {'knee', 'knee angle (deg)'; 'frc', 'vas\_med\_r force (N)'; ...
          'V_ia', 'Ia afferent V (mV)'; 'S', 'MN drive S'};
for k = 1:size(fields, 1)
    nexttile(tl); hold on;
    for mode = ["OFF", "ON"]
        el = R.(mode).(fields{k, 1});
        if strcmp(fields{k, 1}, 'knee')
            y = rad2deg(el.Data(:));
        else
            y = el.Data(:);
        end
        plot(el.Time, y, 'LineWidth', 1.1, 'DisplayName', char(mode));
    end
    grid on; ylabel(fields{k, 2}); legend('Location', 'northeast');
end
exportgraphics(fig, fullfile(outDir, 'e1_reflex.png'), 'Resolution', 130);
save(fullfile(outDir, 'e1_reflex.mat'), 'R');
fprintf('saved %s and e1_reflex.png\n', fullfile(outDir, 'e1_reflex.mat'));
end

% -------------------------------------------------------------------------
function add_neu(mdl, nm, tau, pos)
add_block('SNS_Library/NonSpikingNeuron', [mdl '/' nm], ...
    'Vrest', '0', 'Gm', '1', 'Cm', num2str(1000 * tau), ...
    'Thr', '0', 'Slope', '5', 'Position', pos);
end

function add_syn(mdl, nm, g, esyn, pos)
add_block('SNS_Library/NonSpikingSynapse', [mdl '/' nm], ...
    'gmax', num2str(g, 12), 'Esyn', num2str(esyn), 'ThrPre', '0', ...
    'SlopePre', '5', 'Position', pos);
end

function tau = get_tau(here, nm)
J = jsondecode(fileread(fullfile(here, '..', '..', '..', '..', '..', ...
    'Code', 'MuJoCo_SNS', 'spinal', 'spinal_net_export.json')));
tau = 0.05;
for k = 1:numel(J.neurons)
    if strcmp(J.neurons(k).name, nm)
        tau = J.neurons(k).tau_s;
        return
    end
end
end

function add_mlfn(mdl, name, script, pos)
add_block('simulink/User-Defined Functions/MATLAB Function', ...
          [mdl '/' name], 'Position', pos);
ch = find(sfroot, '-isa', 'Stateflow.EMChart', 'Path', [mdl '/' name]);
ch.Script = script;
end
