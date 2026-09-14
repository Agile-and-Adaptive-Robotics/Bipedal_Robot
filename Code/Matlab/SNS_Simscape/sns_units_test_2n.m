function sns_units_test_2n()
% Units-mapping verification: SNS-Toolbox (numpy) vs SNS_Library Simulink.
%
% Circuit (identical to spinal/export_units_ref_2n.py):
%   5 nA constant -> neuron A (tau 0.2 s) -E synapse (g 0.3 uS, Esyn 8 mV,
%   ThrPre 0, SlopePre 5 mV)-> neuron B (tau 0.5 s)
% Mapping under test (toolbox -> SNS_Library):
%   C = tau uF  ->  Cm = 1000*tau nF ; Gm = 1 uS ; Vrest = 0 mV
%   e_lo/e_hi = 0/5 mV  ->  ThrPre 0 / SlopePre 5
% PASS: max |V_sim - V_ref| < 0.01 mV on both neurons over 3 s.

here = fileparts(mfilename('fullpath'));
mdl = 'sns_units_2n';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(fullfile(here, 'results', [mdl '.slx']), 'file')
    delete(fullfile(here, 'results', [mdl '.slx']));
end
new_system(mdl);

LIB = 'SNS_Library';   % this machine's copy loads on R2025a (verified)
add_block([LIB '/NonSpikingNeuron'], [mdl '/A'], ...
    'Vrest', '0', 'Gm', '1', 'Cm', '200', 'Thr', '0', 'Slope', '5', ...
    'Position', [200 80 300 180]);
add_block([LIB '/NonSpikingNeuron'], [mdl '/B'], ...
    'Vrest', '0', 'Gm', '1', 'Cm', '500', 'Thr', '0', 'Slope', '5', ...
    'Position', [200 300 300 400]);
add_block([LIB '/NonSpikingSynapse'], [mdl '/AtoB'], ...
    'gmax', '0.3', 'Esyn', '8', 'ThrPre', '0', 'SlopePre', '5', ...
    'Position', [60 180 160 280]);
add_block('simulink/Sources/Constant', [mdl '/Iext'], 'Value', '5', ...
    'Position', [60 60 120 100]);

add_line(mdl, 'Iext/1', 'A/1', 'autorouting', 'on');
add_line(mdl, 'A/1', 'AtoB/1', 'autorouting', 'on');    % V -> Vpre
add_line(mdl, 'B/1', 'AtoB/2', 'autorouting', 'on');    % V -> Vpost
add_line(mdl, 'AtoB/1', 'B/1', 'autorouting', 'on');    % Isyn -> B

logv = @(nm, blk, prt) logport(mdl, nm, blk, prt);
logv('va', [mdl '/A'], 1);
logv('vb', [mdl '/B'], 1);

set_param(mdl, 'Solver', 'ode45', 'RelTol', '1e-6', 'AbsTol', '1e-8', ...
          'StopTime', '3', 'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
va = squeeze(out.sigs.get('va').Values.Data);
ta = out.sigs.get('va').Values.Time;
vb = squeeze(out.sigs.get('vb').Values.Data);
tb = out.sigs.get('vb').Values.Time;

save_system(mdl, fullfile(here, 'results', [mdl '.slx']));
close_system(mdl, 0);

ref = load(fullfile(here, 'results', 'units_ref_2n.mat'), 't', 'va', 'vb');
var = interp1(ta, va, ref.t, 'pchip');
vbr = interp1(tb, vb, ref.t, 'pchip');
dva = max(abs(var - ref.va));
dvb = max(abs(vbr - ref.vb));
pass = dva < 0.01 && dvb < 0.01;
fprintf(['UNITS TEST 2-NEURON: %s | max|dV_A| %.2e mV, max|dV_B| %.2e mV ' ...
         '(tol 1e-2 mV); Simulink V_A(end) %.5f, V_B(end) %.5f mV\n'], ...
        tf(pass), dva, dvb, va(end), vb(end));
if ~pass, error('units test FAILED'); end
end

function logport(mdl, name, blk, port)
ph = get_param(blk, 'PortHandles');
set_param(ph.Outport(port), 'DataLogging', 'on', ...
          'DataLoggingNameMode', 'Custom', 'DataLoggingName', name);
end

function s = tf(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
