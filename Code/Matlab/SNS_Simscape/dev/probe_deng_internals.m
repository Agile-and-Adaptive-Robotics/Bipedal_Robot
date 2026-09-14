function probe_deng_internals()
% Tap inside HCNeuron (in-memory library edit, not saved): log m, hInf,
% hint, naDrive, naProd while a 3 nA test current drives the neuron.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
load_system('SNS_Deng_Library');
set_param('SNS_Deng_Library', 'Lock', 'off');
blk = 'SNS_Deng_Library/HCNeuron';
% extra outputs on internal lines
taps = {'mInv', 'm'; 'hint', 'h'; 'naDrive', 'ndrv'; 'naProd', 'nprod'; ...
        'Vint', 'V'};
for k = 1:size(taps, 1)
    src = taps{k, 1};
    if any(strcmp(src, {'mInv', 'hint', 'naDrive', 'naProd', 'Vint'}))
        % these are block outputs; tap their port 1 via a new Outport
        add_block('simulink/Sinks/Out1', [blk '/TAP_' taps{k, 2}], ...
                  'Port', num2str(k + 1), 'Position', [800 40 + 60 * k 830 60 + 60 * k]);
        add_line('SNS_Deng_Library', [blk '/' src '/1'], ...
                 [blk '/TAP_' taps{k, 2} '/1']);
    end
end

mdl = 'probe_hc2';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
add_block('simulink/Sources/Constant', [mdl '/iin'], 'Value', '3', ...
          'Position', [40 60 90 100]);
add_block('SNS_Deng_Library/HCNeuron', [mdl '/hc'], 'Position', [180 40 280 140]);
add_line(mdl, 'iin/1', 'hc/1');
for k = 1:size(taps, 1)
    ph = get_param([mdl '/hc'], 'PortHandles');
    set_param(ph.Outport(k), 'DataLogging', 'on', 'DataLoggingNameMode', ...
              'Custom', 'DataLoggingName', taps{k, 2});
end
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'ode1', ...
          'FixedStep', '1e-4', 'StopTime', '2', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
for k = 1:size(taps, 1)
    ts = out.sigs.get(taps{k, 2}).Values;
    fprintf('%-6s: t=0 %.4f  t=1 %.4f  t=2 %.4f\n', taps{k, 2}, ...
        ts.Data(1), ts.Data(round(numel(ts.Data) / 2)), ts.Data(end));
end
close_system(mdl, 0);
close_system('SNS_Deng_Library', 0);   % discard taps
end
