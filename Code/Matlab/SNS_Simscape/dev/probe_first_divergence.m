function probe_first_divergence()
% Which cell departs first: Simulink ode1@2ms vs numpy coarse 2ms Euler?
here = fileparts(mfilename('fullpath'));
mdl = 'SNS_SpinalNetwork';
load_system(fullfile(here, 'results', [mdl '.slx']));
ref = load(fullfile(here, 'results', 'verify_ref.mat'));

delete_line(mdl, 'u/1', 'demux_u/1');
delete_block([mdl '/u']);
add_block('simulink/Sources/Constant', [mdl '/u_src'], 'Value', 'uvec', ...
          'Position', [40 40 100 120]);
add_line(mdl, 'u_src/1', 'demux_u/1', 'autorouting', 'on');

tn = arrayfun(@(r) strtrim(ref.trace_names(r, :)), ...
              1:size(ref.trace_names, 1), 'UniformOutput', false);
for j = 1:numel(tn)
    ph = get_param([mdl '/' tn{j}], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', ['V_' tn{j}]);
end

uvec = zeros(376, 1); uvec(1) = ref.drive_nA;
assignin('base', 'uvec', uvec);

set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002', 'StopTime', '2', ...
          'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');

% Simulink row k (1-based, t=(k-1)*dt) should equal numpy trc row k.
fprintf('  cell         t=0.10        0.30        0.50        0.70        1.00        1.50        2.00\n');
for j = 1:numel(tn)
    v = out.sigs.get(['V_' tn{j}]).Values;
    vs = v.Data(:);
    rr = ref.traces_coarse(:, j);
    n = min(numel(vs), numel(rr));
    dev = abs(vs(1:n) - rr(1:n));
    tq = [51 151 251 351 501 751 1001];
    fprintf('%-12s %s\n', tn{j}, join(string(compose('%.1e', ...
          arrayfun(@(k) dev(min(k, n)), tq))), '  '));
end
close_system(mdl, 0);
end
