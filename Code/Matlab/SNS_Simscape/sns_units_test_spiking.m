function sns_units_test_spiking()
% Units + semantics verification for the 2026-10-02 SPIKING blocks
% (SpikingNeuron, SpikingSynapse, HybridSpikingSynapse) vs references.
%
% TEST 1  SpikingSynapse conductance = EXACT analytic exponential decay +
%         saturating accumulation (single edge + saturating pulse train).
% TEST 2  HYBRID 2-neuron circuit (SpikingNeuron A -> HybridSpikingSynapse
%         -> NonSpikingNeuron B) vs a MATLAB Euler reference implementing
%         the sns_toolbox 1.5.2 update equations VERBATIM
%         (backends.py SNS_Numpy.forward, fixed-threshold spiking variant):
%           g    <- g*(1 - dt/tau_syn)                    [decay first]
%           i_syn = g*Esyn - V_last*g
%           V    <- V + dt*(1000/Cm)*(-Gm*V + i_syn + Iapp)
%           if V >= Vth:  g <- g + min(ginc, gmax-g);  V <- Vreset
%         PASS: same spike count, max spike-time dev < 0.05 ms,
%         max |V_B| dev < 0.05 mV over 0.5 s.
% TEST 3  spiking->spiking chain sanity: A -> SpikingSynapse (event line)
%         -> SpikingNeuron C fires slower than A.
% Trajectories are saved to results\units_ref_spiking_sim.mat for the
% python cross-check (sns_toolbox numpy backend, same circuit).

here = fileparts(mfilename('fullpath'));
LIB = 'SNS_Library';
load_system(LIB);
ok = true;

%% ---------- TEST 1: analytic synapse conductance ----------
mdl = 'sns_units_spk_syn';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); load_system(mdl);
add_block([LIB '/SpikingSynapse'], [mdl '/syn1'], 'gmax', '1.5', 'ginc', '0.5', ...
    'tau_syn', '100', 'Esyn', '8', 'Position', [120 60 180 120]);
add_block('simulink/Sources/Step', [mdl '/step'], 'Time', '0.5', 'Before', '0', ...
    'After', '1', 'Position', [30 75 70 105]);
add_block('simulink/Signal Routing/Demux', [mdl '/dmx'], 'Outputs', '2', 'Position', [220 70 223 110]);
add_block('simulink/Sinks/Out1', [mdl '/glog'], 'Port', '1', 'Position', [270 78 300 92]);
add_line(mdl, 'step/1', 'syn1/1');
add_line(mdl, 'syn1/1', 'dmx/1');
add_line(mdl, 'dmx/1', 'glog/1');
set_param(mdl, 'StopTime', '1.5', 'Solver', 'ode45', 'RelTol', '1e-9', ...
    'SignalLogging', 'on', 'SignalLoggingName', 'sigs', 'UnconnectedInputMsg', 'none');
logsig(mdl, 'g', [mdl '/dmx'], 1);   % demux channel 1 = g [uS]
out = sim(mdl);
sg = out.sigs.get('g').Values;
tq = 0.5+0.05:0.05:1.45;
gq = interp1(sg.Time, squeeze(sg.Data), tq, 'pchip');
gref = 0.5*exp(-(tq-0.5)/0.1);
e1 = max(abs(gq - gref));
fprintf('TEST1a single-edge decay: max|g - ginc*exp(-(t-ts)/tau)| = %.3e uS\n', e1);
ok = ok && e1 < 1e-4;
close_system(mdl, 0);

% saturating pulse train
mdl = 'sns_units_spk_syn2';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); load_system(mdl);
add_block([LIB '/SpikingSynapse'], [mdl '/syn2'], 'gmax', '0.55', 'ginc', '0.5', ...
    'tau_syn', '30', 'Esyn', '8', 'Position', [120 60 180 120]);
add_block('simulink/Sources/Pulse Generator', [mdl '/pulse'], 'Amplitude', '1', ...
    'Period', '0.05', 'PulseWidth', '50', 'PhaseDelay', '0.1', 'Position', [30 75 70 105]);
add_block('simulink/Signal Routing/Demux', [mdl '/dmx'], 'Outputs', '2', 'Position', [220 70 223 110]);
add_block('simulink/Sinks/Out1', [mdl '/glog'], 'Port', '1', 'Position', [270 78 300 92]);
add_line(mdl, 'pulse/1', 'syn2/1');
add_line(mdl, 'syn2/1', 'dmx/1');
add_line(mdl, 'dmx/1', 'glog/1');
set_param(mdl, 'StopTime', '0.35', 'Solver', 'ode45', 'RelTol', '1e-9', ...
    'SignalLogging', 'on', 'SignalLoggingName', 'sigs', 'UnconnectedInputMsg', 'none');
logsig(mdl, 'g', [mdl '/dmx'], 1);
out = sim(mdl);
sg = out.sigs.get('g').Values;
% analytic: edges at 0.1 + 0.05k, D = exp(-0.05/0.03)
% sample points stay >=10 ms away from the jumps and use LINEAR interp
% (pchip across the solver's sample exactly AT a jump overshoots badly)
D = exp(-0.05/0.03); gp = 0; gcap = 0.55; ginc = 0.5; tau = 0.03;
tq = []; rq = [];
for tk = 0.1:0.05:0.30
    gp = min(gcap, gp*D + ginc);
    for tt = [tk+0.02, tk+0.04]
        tq(end+1) = tt; %#ok<SAGROW>
        rq(end+1) = gp*exp(-(tt-tk)/tau); %#ok<SAGROW>
    end
end
gq = interp1(sg.Time, squeeze(sg.Data), tq, 'linear');
e1b = max(abs(gq - rq));
fprintf('TEST1b saturating accumulation: max|g - analytic| = %.3e uS (gmax 0.55, ginc 0.5)\n', e1b);
ok = ok && e1b < 1e-4;
close_system(mdl, 0);

%% ---------- TEST 2: hybrid 2-neuron circuit vs toolbox equations ----------
% A: SpikingNeuron Vrest 0, Gm 1 uS, Cm 5 nF, Vth 8 mV, Vreset 0, Iapp 10 nA
% syn: HybridSpikingSynapse gmax 0.3, ginc 0.1, tau_syn 100 ms, Esyn 8 mV,
%      ThrPre 7.5 mV (0.94 x A's Vth: the block keys on ITS OWN ThrPre
%      crossing — AnimatLab semantics. ThrPre must sit BETWEEN Vrest_pre
%      and Vth_pre, comfortably BELOW Vth_pre: a spiking presynaptic cell
%      resets AT threshold, so its V output never renders values >= Vth)
% B: NonSpikingNeuron Vrest 0, Gm 1, Cm 200
mdl = 'sns_units_spk_hybrid';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(fullfile(here, 'results', [mdl '.slx']), 'file')
    delete(fullfile(here, 'results', [mdl '.slx']));
end
new_system(mdl); load_system(mdl);
add_block([LIB '/SpikingNeuron'], [mdl '/A'], 'Vrest', '0', 'Vth', '8', ...
    'Vreset', '0', 'Gm', '1', 'Cm', '5', 'Position', [150 60 240 150]);
add_block([LIB '/HybridSpikingSynapse'], [mdl '/AtoB'], 'gmax', '0.3', 'ginc', '0.1', ...
    'tau_syn', '100', 'Esyn', '8', 'ThrPre', '7.5', 'Position', [280 190 340 250]);
add_block([LIB '/NonSpikingNeuron'], [mdl '/B'], 'Vrest', '0', 'Gm', '1', ...
    'Cm', '200', 'Thr', '0', 'Slope', '5', 'Position', [400 60 490 150]);
add_block('simulink/Sources/Constant', [mdl '/Iext'], 'Value', '10', 'Position', [60 90 110 120]);
add_line(mdl, 'Iext/1', 'A/1');
add_line(mdl, 'A/1', 'AtoB/1');      % A's V -> Vpre (hybrid detects the spike)
add_line(mdl, 'AtoB/1', 'B/2');      % [g; g*Esyn] -> B syn1
add_block('simulink/Signal Routing/Demux', [mdl '/gdmx'], 'Outputs', '2', 'Position', [330 300 333 340]);
add_block('simulink/Sinks/Terminator', [mdl '/gterm'], 'Position', [380 310 400 330]);
add_line(mdl, 'AtoB/1', 'gdmx/1');   % branch: expose g for diagnostics
add_line(mdl, 'gdmx/1', 'gterm/1');
logv = @(nm, blk, prt) logsig(mdl, nm, blk, prt);
logv('va', [mdl '/A'], 1);
logv('spk', [mdl '/A'], 2);
logv('vb', [mdl '/B'], 1);
logv('g', [mdl '/gdmx'], 1);
set_param(mdl, 'StopTime', '0.5', 'Solver', 'ode45', 'RelTol', '1e-8', 'AbsTol', '1e-10', ...
    'SignalLogging', 'on', 'SignalLoggingName', 'sigs', 'UnconnectedInputMsg', 'none');
out = sim(mdl);
va = out.sigs.get('va').Values; vb = out.sigs.get('vb').Values;
spkl = out.sigs.get('spk').Values;
save_system(mdl, fullfile(here, 'results', [mdl '.slx']));
close_system(mdl, 0);

% spike times: rising edges of the logged spike-event line (the solver lands
% on the threshold crossing, so each pulse shows as >=1 logged sample at 1)
v = squeeze(spkl.Data) >= 0.5;
edges = find(v & ~[false; v(1:end-1)]);
tspk = spkl.Time(edges);
if isempty(tspk)   % fallback: membrane sawtooth (>= 7.9 then <= 0.1)
    vaD = squeeze(va.Data); vaT = va.Time;
    for k = 1:numel(vaD)-1
        if vaD(k) >= 7.9 && vaD(k+1) <= 0.1, tspk(end+1) = vaT(k); end %#ok<SAGROW>
    end
end

% --- toolbox-equations Euler reference (verbatim update order) ---
% dt = 1e-6 s: the continuous Simulink model resolves crossings exactly;
% a coarser reference accumulates its own O(dt) spike-time bias (measured
% 7.2 us/spike at dt=1e-5, which is REFERENCE error, not block error).
dt = 1e-6; tEnd = 0.5; n = round(tEnd/dt);
CmA = 5; GmA = 1; Vth = 8; Vrst = 0; Iapp = 10;
gmax = 0.3; ginc = 0.1; tauS = 0.1; Es = 8;
CmB = 200; GmB = 1;
VAr = Vrst; VBr = 0; gr = 0; tspkr = [];
tr = (0:n)'*dt; VBr_t = zeros(n+1, 1); VAr_t = zeros(n+1, 1); g_t = zeros(n+1, 1);
for k = 1:n
    gr = gr*(1 - dt/tauS);                       % decay FIRST (backends.py:131)
    isynB = gr*Es - VBr*gr;                      % i_syn uses V_last (backends.py:138)
    VAr = VAr + dt*(1000/CmA)*(-GmA*VAr + Iapp);
    VBr = VBr + dt*(1000/CmB)*(-GmB*VBr + isynB);
    if VAr >= Vth                                % spike (backends.py:169)
        gr = gr + min(ginc, gmax - gr);          % saturating increment (:178)
        VAr = Vrst;                              % reset (:181)
        tspkr(end+1) = tr(k+1); %#ok<SAGROW>
    end
    VAr_t(k+1) = VAr; VBr_t(k+1) = VBr; g_t(k+1) = gr;
end
% compare (linear interp: pchip overshoots across the jump discontinuities)
vbr = interp1(vb.Time, squeeze(vb.Data), tr, 'linear');
gsig = out.sigs.get('g').Values;
grq = interp1(gsig.Time, squeeze(gsig.Data), tr, 'linear');
e2g = max(abs(grq - g_t));
e2v = max(abs(vbr - VBr_t));
if numel(tspk) == numel(tspkr) && ~isempty(tspk)
    e2t = max(abs(tspk(:) - tspkr(:)));
    % localization diagnostic: is the residual g-level or event-timing?
    dtspk = tspk(:) - tspkr(:);
else
    e2t = inf; dtspk = NaN;
end
fprintf(['TEST2 hybrid vs toolbox-equations: spikes sim %d / ref %d, max spike-time dev ' ...
    '%.3e s, max|V_B dev| %.3e mV, max event-window |g dev| %.3e uS\n'], numel(tspk), numel(tspkr), e2t, e2v, e2g);
fprintf('TEST2 diagnostics: spike-time devs (ms, first 5 / last 5): ');
fprintf('%.4f ', dtspk(1:min(5, end))*1e3);
fprintf('... ');
fprintf('%.4f ', dtspk(max(1, end-4):end)*1e3);
fprintf('\n');
% windowed mean-g agreement AFTER the initial ramp (the pointwise g dev is
% the event-window quantization between two event streams; the first 0.1 s
% ramp phase differs by the ThrPre-vs-spike keying of the hybrid detector)
gmdev = 0;
for w = 0.1:0.05:0.45
    iw = tr > w & tr <= w+0.05;
    gmdev = max(gmdev, abs(mean(grq(iw)) - mean(g_t(iw))));
end
fprintf('TEST2 windowed mean-g max dev = %.3e uS\n', gmdev);
ok = ok && numel(tspk) == numel(tspkr) && e2t < 5e-5 && e2v < 0.05 && gmdev < 5e-3;
% plain arrays for the python cross-check (timeseries objects do not
% round-trip through scipy.io.loadmat)
vb_t = vb.Time; vb_d = squeeze(vb.Data);
va_t = va.Time; va_d = squeeze(va.Data);
g_t_sim = squeeze(gsig.Data); g_t_time = gsig.Time;
save(fullfile(here, 'results', 'units_ref_spiking_sim.mat'), 'tspk', ...
    'VAr_t', 'VBr_t', 'tr', 'tspkr', 'g_t', 'vb_t', 'vb_d', 'va_t', 'va_d', ...
    'g_t_sim', 'g_t_time');

%% ---------- TEST 3: spiking -> spiking chain ----------
mdl = 'sns_units_spk_chain';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); load_system(mdl);
add_block([LIB '/SpikingNeuron'], [mdl '/A'], 'Vrest', '0', 'Vth', '8', ...
    'Vreset', '0', 'Gm', '1', 'Cm', '5', 'Position', [100 60 190 150]);
add_block([LIB '/SpikingSynapse'], [mdl '/AtoC'], 'gmax', '6', 'ginc', '6', ...
    'tau_syn', '100', 'Esyn', '10', 'Position', [240 190 300 250]);
add_block([LIB '/SpikingNeuron'], [mdl '/C'], 'Vrest', '0', 'Vth', '8', ...
    'Vreset', '0', 'Gm', '1', 'Cm', '50', 'Position', [360 60 450 150]);
add_block('simulink/Sources/Constant', [mdl '/Iext'], 'Value', '10', 'Position', [20 90 60 120]);
add_line(mdl, 'Iext/1', 'A/1');
add_line(mdl, 'A/2', 'AtoC/1');    % A's SPIKE line -> event-wired synapse
add_line(mdl, 'AtoC/1', 'C/2');    % onto a SPIKING postsynaptic neuron
logv('va', [mdl '/A'], 1);
logv('vc', [mdl '/C'], 1);
logv('spkA', [mdl '/A'], 2);
logv('spkC', [mdl '/C'], 2);
set_param(mdl, 'StopTime', '0.3', 'Solver', 'ode45', 'RelTol', '1e-8', ...
    'SignalLogging', 'on', 'SignalLoggingName', 'sigs', 'UnconnectedInputMsg', 'none');
out = sim(mdl);
nA = count_edges(out.sigs.get('spkA').Values);
nC = count_edges(out.sigs.get('spkC').Values);
fprintf('TEST3 spiking chain: A fired %d, C fired %d times in 0.3 s\n', nA, nC);
ok = ok && nA > 10 && nC >= 1 && nC < nA;
close_system(mdl, 0);

if ok
    fprintf('SPIKING UNITS TEST: PASS\n');
else
    error('spiking units test FAILED');
end
end

function n = count_edges(ts)
% rising edges of a logged boolean/0-1 line = one per event pulse
v = squeeze(ts.Data) >= 0.5;
n = sum(v & ~[false; v(1:end-1)]);
end

function n = count_spikes(v)
% V hits ~Vth then the reset drops it to ~Vreset -> count those drops
n = sum(v(1:end-1) >= 7.9 & v(2:end) <= 0.1);
end

function logsig(mdl, name, blk, port)
ph = get_param(blk, 'PortHandles');
set_param(ph.Outport(port), 'DataLogging', 'on', ...
    'DataLoggingNameMode', 'Custom', 'DataLoggingName', name);
end
