function probe_deng_inside()
% Log INTERNAL signals of the linked HCNeuron inside SNS_Deng_RG (zero
% stimulus). Disables the library link on the copy (model closed w/o save).
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
load_system(fullfile(here, 'SNS_Deng_RG.slx'));

blk = 'SNS_Deng_RG/pair1_HC_ext';
set_param(blk, 'LinkStatus', 'inactive');
% zero the stimulus
delete_line('SNS_Deng_RG', 'I_stim/1', 'pair1_sum_ext/2');
add_block('simulink/Sources/Constant', 'SNS_Deng_RG/z', 'Value', '0', ...
          'Position', [480 660 520 700]);
add_line('SNS_Deng_RG', 'z/1', 'pair1_sum_ext/2');

taps = {'mInv', 'm'; 'hInf', 'hinf'; 'hint', 'h'; 'naDrive', 'ndrv'; ...
        'dVsum', 'dV'; 'Vint', 'V'};
for k = 1:size(taps, 1)
    ph = get_param([blk '/' taps{k, 1}], 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', ['tap_' taps{k, 2}]);
end
set_param('SNS_Deng_RG', 'StopTime', '2', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');
out = sim('SNS_Deng_RG', 'ReturnWorkspaceOutputs', 'on');
for k = 1:size(taps, 1)
    ts = out.sigs.get(['tap_' taps{k, 2}]).Values;
    fprintf('%-5s: t=0 %+.4f  t=1 %+.4f  t=2 %+.4f\n', taps{k, 2}, ...
        ts.Data(1), ts.Data(round(numel(ts.Data) / 2)), ts.Data(end));
end
close_system('SNS_Deng_RG', 0);
end
