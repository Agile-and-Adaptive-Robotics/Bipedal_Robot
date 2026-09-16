function probe_deng_min()
% Minimal HCNeuron test: Constant 0 -> HCNeuron, 10 ms. Bisect: GNa=0 vs
% full, tauH huge vs 350. Print V trajectory samples.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
load_system('SNS_Deng_Library');

mdl = 'probe_hc';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
add_block('simulink/Sources/Constant', [mdl '/zero'], 'Value', '0', ...
          'Position', [40 60 90 100]);
add_block('SNS_Deng_Library/HCNeuron', [mdl '/hc'], 'Position', [180 40 280 140]);
add_line(mdl, 'zero/1', 'hc/1');
ph = get_param([mdl '/hc'], 'PortHandles');
set_param(ph.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', ...
          'Custom', 'DataLoggingName', 'v');
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '1e-4', 'StopTime', '5', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');

for cfg = {'full', 'noNa', 'noM'}
    switch cfg{1}
        case 'full',  set_param([mdl '/hc'], 'GNa', '1.5', 'Sm', '0.2', 'tauH', '350');
        case 'noNa',  set_param([mdl '/hc'], 'GNa', '0');
        case 'noM',   set_param([mdl '/hc'], 'GNa', '1.5', 'Sm', '0');
    end
    try
        out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
        vts = out.sigs.get('v').Values; t = vts.Time(:); y = vts.Data(:);
        
        mk = @(ts) interp1(t, y, ts, 'linear', 'extrap');
        fprintf('%-5s: V(0)=%.2f V(0.5)=%.2f V(1)=%.2f V(2)=%.2f V(3)=%.2f V(5)=%.2f  range [%.1f, %.1f]\n', ...
            cfg{1}, y(1), mk(0.5), mk(1), mk(2), mk(3), mk(5), min(y), max(y));
    catch e
        fprintf('%-5s: FAILED %s\n', cfg{1}, strtok(e.message, newline));
    end
end
close_system(mdl, 0);
end
