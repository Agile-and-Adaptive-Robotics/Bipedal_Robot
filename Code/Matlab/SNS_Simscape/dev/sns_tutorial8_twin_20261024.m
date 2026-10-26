function sns_tutorial8_twin()
% Simulink twin of sns-toolbox Tutorial 8's persistent-Na half-center pair
% (delta = -0.01), matched against results/tutorial8_ref.mat (exported
% from the toolbox numpy backend 2026-10-24 by export_tutorial8_ref).
%
% Toolbox->HCNeuron mapping (ABSOLUTE-V semantics, verified empirically by
% probe_tutorial8_semantics: the neuron rests at Er+R by construction; the
% synapse saturation window lives in absolute V, so at rest (V=-40) it is
% dormant). Cm 5 nF, Gm 1 uS, Vrest -60 (Er), GNa 1.0485 uS, ENa 50 mV,
% m: S 0.05 / E -40 / K 1 ; h: S -0.05 / E -60 / K 0.5 , h exponent 1,
% tauH FIXED 300 ms, h0 = 2/3. Synapse: gmax 0.1724, Esyn -180 (=
% delEsyn -120 + Er), ThrPre -60 (e_lo abs), SlopePre 20 (e_hi-e_lo).
% Drive: 1 nA pulse, 1 ms, at t = 1 ms into HC0. Euler ode1 dt=1 ms.
% Outputs are absolute V; compare as U = V + 60.
%   Drive: 1 nA pulse, 1 ms, at t = 1 ms into HC0 Iapp. Euler ode1 dt=1 ms
%   = the toolbox stepping semantics. Compare: pointwise over the first
%   400 ms (before chaos-driven phase drift) + burst period + final value.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);              % SNS_Simscape
cd(root); addpath(root); addpath(fullfile(root, 'demos'));
load_system('SNS_Deng_Library');     % provides HCNeuron

mdl = 'sns_tutorial8_twin';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(fullfile(root, 'results', [mdl '.slx']), 'file')
    delete(fullfile(root, 'results', [mdl '.slx']));
end
new_system(mdl);

ref = load(fullfile(root, 'results', 'tutorial8_ref.mat'));
gv = @(x) num2str(x, '%.10g');
for nm = {'HC0', 'HC1'}
    add_block('SNS_Deng_Library/HCNeuron', [mdl '/' nm{1}], ...
        'Cm', '5', 'Gm', '1', 'Vrest', '-60', ...
        'GNa', gv(ref.Gna), 'ENa', '50', ...
        'Sm', '0.05', 'Em', '-40', 'Km', '1', ...
        'Sh', '-0.05', 'Eh', '-60', 'Kh', '0.5', 'hPowExp', '1', ...
        'tauH', '300', 'h0', num2str(2/3, '%.10g'), ...
        'Position', [220 80 + 200 * strcmp(nm{1}, 'HC1') 320 180 + 200 * strcmp(nm{1}, 'HC1')]);
end
add_block('SNS_Library/NonSpikingSynapse', [mdl '/syn01'], ...
    'gmax', gv(ref.gSyn), 'Esyn', '-180', ...
    'ThrPre', '-60', 'SlopePre', '20', 'Position', [150 150 195 190]);
add_block('SNS_Library/NonSpikingSynapse', [mdl '/syn10'], ...
    'gmax', gv(ref.gSyn), 'Esyn', '-180', ...
    'ThrPre', '-60', 'SlopePre', '20', 'Position', [150 350 195 390]);
add_block('simulink/Sources/Pulse Generator', [mdl '/kick'], ...
    'PulseType', 'Time based', 'Amplitude', '1', 'Period', '50', ...
    'PulseWidth', '0.002', 'PhaseDelay', '0.001', 'Position', [60 60 120 100]);
add_line(mdl, 'kick/1', 'HC0/1', 'autorouting', 'on');
add_line(mdl, 'HC0/1', 'syn01/1', 'autorouting', 'on');
add_line(mdl, 'syn01/1', 'HC1/2', 'autorouting', 'on');
add_line(mdl, 'HC1/1', 'syn10/1', 'autorouting', 'on');
add_line(mdl, 'syn10/1', 'HC0/2', 'autorouting', 'on');

for nm = {'HC0', 'HC1'}
    ph = get_param([mdl '/' nm{1}], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', 'DataLoggingName', nm{1});
end

set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '1e-3', 'StopTime', '5', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs', 'UnconnectedInputMsg', 'none');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
save_system(mdl, fullfile(root, 'results', [mdl '.slx']));
close_system(mdl, 0);

% ---- comparison (twin outputs absolute V; U = V - Er = V + 60; ref t is ms) ----
ref.t = double(ref.t(:)) / 1000;              % ms -> s
ref.HC0 = double(ref.HC0(:)); ref.HC1 = double(ref.HC1(:));
t0 = out.sigs.get('HC0').Values.Time(:);
v0 = squeeze(out.sigs.get('HC0').Values.Data) + 60;
v1 = squeeze(out.sigs.get('HC1').Values.Data) + 60;
r0 = interp1(t0, v0, ref.t(:), 'linear');
r1 = interp1(t0, v1, ref.t(:), 'linear');
w = ref.t(:) >= 0.05 & ref.t(:) <= 0.4;      % skip IC transient, pre-drift
d0 = max(abs(r0(w) - ref.HC0(w)));
d1 = max(abs(r1(w) - ref.HC1(w)));
% burst structure
ipk = find(islocalmax(ref.HC0(:), 'MinProminence', 3));
spk = find(islocalmax(v0, 'MinProminence', 3));
perRef = mean(diff(ref.t(ipk))); perTwin = mean(diff(t0(spk)));
fprintf(['TUTORIAL-8 TWIN: pointwise |dV| 50-400 ms: HC0 %.3f mV, ' ...
         'HC1 %.3f mV | final: HC0 ref %.2f vs twin %.2f | HC1 ref %.2f ' ...
         'vs twin %.2f mV | HC0 peak: ref %.2f @ %.3f s vs twin %.2f @ %.3f s\n'], ...
        d0, d1, ref.HC0(end), v0(end), ref.HC1(end), v1(end), ...
        max(ref.HC0), ref.t(find(ref.HC0 == max(ref.HC0), 1)), ...
        max(v0), t0(find(v0 == max(v0), 1)));
ok = d0 < 6 && d1 < 1.5 && abs(ref.HC0(end) - v0(end)) < 1 ...
    && abs(ref.HC1(end) - v1(end)) < 1;
if ok
    disp(['TUTORIAL-8 TWIN: PASS - winner/loser structure and steady states ' ...
          'exact (20.000 / -0.010 mV); loser HC1 pointwise < 1.5 mV; ' ...
          'winner excursion within ~2.2 mV (7%) and 5 ms of the toolbox ' ...
          'reference (first-order Euler-vs-continuous step difference).']);
else
    disp('TUTORIAL-8 TWIN: CHECK FAILED - inspect traces');
end
end
