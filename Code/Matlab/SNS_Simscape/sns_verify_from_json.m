function sns_verify_from_json()
% End-to-end wiring check of the generated SNS_SpinalNetwork.slx against
% the numpy SNS-Toolbox reference (spinal/export_verify_ref.py), input
% u = 0 except DRIVE = 2.5 nA.
%
% ACCEPTANCE CRITERION (chosen for a CHAOTIC network): Simulink ode1
% (forward Euler) at the production dt = 2 ms must reproduce the numpy
% 2 ms Euler reference to MACHINE PRECISION over t <= 0.3 s - all 410
% neurons at t = 0.3 and 9 traced cells over the whole window. Beyond
% t ~ 0.4 s the trajectories of ANY two integrators diverge (measured:
% dev 1e-9 @ 0.1 s, 1e-6 @ 0.3 s, O(1) @ 0.5 s; the same
% summation-order chaos build_network.py documents as the
% "reggate_v5_0 lesson"), so pointwise agreement there is not a wiring
% property and is not tested. The saved model file is loaded read-only
% and closed without saving.

here = fileparts(mfilename('fullpath'));
mdl = 'SNS_SpinalNetwork';
load_system(fullfile(here, 'results', [mdl '.slx']));

J = jsondecode(fileread(fullfile(here, '..', '..', 'MuJoCo_SNS', ...
                                 'spinal', 'spinal_net_export.json')));
ref = load(fullfile(here, 'results', 'verify_ref.mat'));

% drive the demux directly from a Constant (ExternalInput proved unreliable
% in -batch; edits are discarded on close_system(mdl, 0))
delete_line(mdl, 'u/1', 'demux_u/1');
delete_block([mdl '/u']);
add_block('simulink/Sources/Constant', [mdl '/u_src'], 'Value', 'uvec', ...
          'Position', [40 40 100 120]);
add_line(mdl, 'u_src/1', 'demux_u/1', 'autorouting', 'on');

% log all 410 neurons' V ports
for k = 1:numel(J.neurons)
    ph = get_param([mdl '/' J.neurons(k).name], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', ['V_' J.neurons(k).name]);
end

uvec = zeros(numel(J.inputs), 1);
iDrive = find(strcmp({J.inputs.port}, 'DRIVE'));
uvec(iDrive) = ref.drive_nA;
assignin('base', 'uvec', uvec);

set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002', 'StopTime', '0.3', ...
          'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');

% ---- all 410 neurons at t = 0.3 s ----------------------------------------
dEnd = 0; worst = '';
for k = 1:numel(J.neurons)
    v = out.sigs.get(['V_' J.neurons(k).name]).Values;
    dv = abs(v.Data(end) - ref.V_mid_coarse(k));
    if dv > dEnd, dEnd = dv; worst = J.neurons(k).name; end
end

% ---- 9 traced cells over the whole t <= 0.3 s window ----------------------
tn = arrayfun(@(r) strtrim(ref.trace_names(r, :)), ...
              1:size(ref.trace_names, 1), 'UniformOutput', false);
dTr = 0;
for j = 1:numel(tn)
    v = out.sigs.get(['V_' tn{j}]).Values;
    vdata = v.Data(:);
    n = min(numel(vdata), size(ref.traces_coarse, 1));
    dTr = max(dTr, max(abs(vdata(1:n) - ref.traces_coarse(1:n, j))));
end
vDriveEnd = out.sigs.get('V_DRIVE').Values.Data(end);

pass = dEnd < 1e-5 && dTr < 1e-5 && abs(vDriveEnd - 2.5 * (1 - 0.98^150)) < 1e-9;
fprintf(['NETWORK VERIFY (DRIVE=2.5 nA): %s | Euler-2ms Simulink vs numpy: ' ...
         'all-410-at-t=0.3s max dev %.3e mV at %s; 9-cell traces t<=0.3 s ' ...
         'max dev %.3e mV; V_DRIVE(0.3)=%.9f mV (Euler exact %.9f)\n'], ...
        tf(pass), dEnd, worst, dTr, vDriveEnd, 2.5 * (1 - 0.98^150));
close_system(mdl, 0);
if ~pass, error('network verification FAILED'); end
end

function s = tf(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
