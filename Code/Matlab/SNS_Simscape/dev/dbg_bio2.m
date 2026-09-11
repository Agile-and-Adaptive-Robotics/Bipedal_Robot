% dbg_bio2.m — fully instrumented BioMuscle check, no try/catch hiding anything
here = fileparts(mfilename('fullpath'));
cd(here); addpath(here);
sns_build_actuators;
load_system('SNS_Library');

mdl = 'dbg2';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '0.01');

add_block('simulink/Sources/Constant', [mdl '/A_c'], 'Value', '1', 'Position', [20 20 60 50]);
add_block('simulink/Sources/Constant', [mdl '/L_c'], 'Value', '0.15', 'Position', [20 80 60 110]);
add_block('simulink/Sources/Constant', [mdl '/V_c'], 'Value', '0', 'Position', [20 140 60 170]);
add_block('SNS_Library/BioMuscle', [mdl '/bio'], 'Position', [120 40 260 180]);
fprintf('ReferenceBlock before: %s\n', get_param([mdl '/bio'], 'ReferenceBlock'));
set_param([mdl '/bio'], 'LinkStatus', 'inactive');
n = numel(find_system([mdl '/bio'], 'SearchDepth', 1, 'LookUnderMasks', 'all', 'Type', 'Block'));
fprintf('blocks after inactive: %d\n', n);

add_block('simulink/Sinks/To Workspace', [mdl '/logF'], 'VariableName', 'Fbio', 'SaveFormat', 'Timeseries', 'Position', [320 90 400 120]);
add_line(mdl, 'A_c/1', 'bio/1', 'autorouting', 'on');
add_line(mdl, 'L_c/1', 'bio/2', 'autorouting', 'on');
add_line(mdl, 'V_c/1', 'bio/3', 'autorouting', 'on');
add_line(mdl, 'bio/1', 'logF/1', 'autorouting', 'on');

taps = {'lam', 'flExp', 'vSwitch', 'concA', 'pSat', 'pAdd', 'FmaxG', 'cosMul', ...
        'lsub', 'lmax', 'vDiv', 'vNeg', 'bSum', 'bAb', 'concDiv', 'cosA', 'lm1', 'plm1', 'pSub', 'pNorm'};
for k = 1:numel(taps)
    lbl = ['tap_' taps{k}];
    bioSys = [mdl '/bio'];
    add_block('simulink/Sinks/To Workspace', [bioSys '/' lbl], 'VariableName', lbl, ...
        'SaveFormat', 'Timeseries', 'Position', [500 40+40*k 580 70+40*k]);
    sph = get_param([bioSys '/' taps{k}], 'PortHandles').Outport(1);
    dph2 = get_param([bioSys '/' lbl], 'PortHandles').Inport(1);
    add_line(bioSys, sph, dph2, 'autorouting', 'on');
end
out = sim(mdl);
fprintf('F = %.3f (expect 1000)\n', out.Fbio.Data(end));
for k = 1:numel(taps)
    fprintf('  %-8s = %.5g\n', taps{k}, out.(['tap_' taps{k}]).Data(end));
end
close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
fprintf('DBG2 DONE\n');
