function probe_euler_match()
% Focused probe: does the generated model under ode1 @ 2 ms reproduce the
% numpy coarse (2 ms Euler) reference? Prints sample values side by side.
here = fileparts(mfilename('fullpath'));
mdl = 'SNS_SpinalNetwork';
load_system(fullfile(here, 'results', [mdl '.slx']));
ref = load(fullfile(here, 'results', 'verify_ref.mat'));

delete_line(mdl, 'u/1', 'demux_u/1');
delete_block([mdl '/u']);
add_block('simulink/Sources/Constant', [mdl '/u_src'], 'Value', 'uvec', ...
          'Position', [40 40 100 120]);
add_line(mdl, 'u_src/1', 'demux_u/1', 'autorouting', 'on');

for nm = ["DRIVE", "RG_E_r", "PF_E1_r"]
    ph = get_param([mdl '/' char(nm)], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', ['V_' char(nm)]);
end

uvec = zeros(376, 1); uvec(1) = ref.drive_nA;
assignin('base', 'uvec', uvec);

sp = get_param(mdl, 'SolverType');
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '0.002', 'StopTime', '2', ...
          'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
fprintf('solver now: %s / %s / step %s\n', get_param(mdl, 'SolverType'), ...
        get_param(mdl, 'Solver'), get_param(mdl, 'FixedStep'));
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');

tn = arrayfun(@(r) strtrim(ref.trace_names(r, :)), ...
              1:size(ref.trace_names, 1), 'UniformOutput', false);
jRG = find(strcmp(tn, 'RG_E_r')); jPF = find(strcmp(tn, 'PF_E1_r'));
vD = out.sigs.get('V_DRIVE').Values;
vR = out.sigs.get('V_RG_E_r').Values;
vP = out.sigs.get('V_PF_E1_r').Values;
fprintf('  t     V_DRIVE(sim/ref)              V_RG_E_r(sim/ref)            V_PF_E1_r(sim/ref)\n');
for q = [1 26 51 101 201 401 801 1001]
    tq = vD.Time(q);
    rD = interp1(ref.t_coarse, 2.5 * ones(size(ref.t_coarse)), tq);
    rR = interp1(ref.t_coarse, ref.traces_coarse(:, jRG), tq, 'previous');
    rP = interp1(ref.t_coarse, ref.traces_coarse(:, jPF), tq, 'previous');
    fprintf('%.3f  %+.6f / %+.6f   %+ .6f / %+.6f   %+ .6f / %+.6f\n', ...
            tq, vD.Data(q), rD, vR.Data(q), rR, vP.Data(q), rP);
end
close_system(mdl, 0);
end
