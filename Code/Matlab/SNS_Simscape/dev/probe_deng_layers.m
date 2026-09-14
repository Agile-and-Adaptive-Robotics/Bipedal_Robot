function probe_deng_layers()
% Isolate: simulate SNS_Deng_RG alone (zero stimulus) and SNS_Deng_PF alone
% (grounded RG inports), 0.05 s, print V ranges to find the blow-up.
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

% --- RG alone ---
load_system(fullfile(here, 'SNS_Deng_RG.slx'));
delete_line('SNS_Deng_RG', 'I_stim/1', 'pair1_sum_ext/2');
add_block('simulink/Sources/Constant', 'SNS_Deng_RG/z', 'Value', '0', ...
          'Position', [480 660 520 700]);
add_line('SNS_Deng_RG', 'z/1', 'pair1_sum_ext/2');
ph = get_param('SNS_Deng_RG/pair1_HC_ext', 'PortHandles');
set_param(ph.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', ...
          'Custom', 'DataLoggingName', 'vE');
ph = get_param('SNS_Deng_RG/pair1_HC_flx', 'PortHandles');
set_param(ph.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', ...
          'Custom', 'DataLoggingName', 'vF');
set_param('SNS_Deng_RG', 'StopTime', '0.05');
try
    out = sim('SNS_Deng_RG', 'ReturnWorkspaceOutputs', 'on');
    vE = out.sigs.get('vE'); vF = out.sigs.get('vF');
    fprintf('RG alone: HC_ext [%.2f, %.2f] mV, HC_flx [%.2f, %.2f] mV\n', ...
        min(vE.Data), max(vE.Data), min(vF.Data), max(vF.Data));
catch e
    fprintf('RG alone FAILED: %s\n', e.message);
end
close_system('SNS_Deng_RG', 0);

% --- PF alone (RG inports unconnected -> grounded) ---
load_system(fullfile(here, 'SNS_Deng_PF.slx'));
phE = get_param('SNS_Deng_PF/hip_HC_ext', 'PortHandles');
set_param(phE.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', ...
          'Custom', 'DataLoggingName', 'vH');
set_param('SNS_Deng_PF', 'StopTime', '0.05');
try
    out = sim('SNS_Deng_PF', 'ReturnWorkspaceOutputs', 'on');
    vH = out.sigs.get('vH');
    fprintf('PF alone: hip HC_ext [%.2f, %.2f] mV\n', ...
        min(vH.Data), max(vH.Data));
catch e
    fprintf('PF alone FAILED: %s\n', e.message);
end
close_system('SNS_Deng_PF', 0);
end
