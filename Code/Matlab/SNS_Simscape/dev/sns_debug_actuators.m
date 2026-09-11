%% sns_debug_actuators.m — signal-level debug of BPA_10mm + BioMuscle internals
% Instances the blocks at fixed inputs, prints block outputs vs closed-form
% references, and taps BioMuscle internal signals.
cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);
addpath(fullfile(cdto, '..', 'Functions'));

load_system('SNS_Library');
mdl = 'dbgAct';

%% ---- BPA_10mm at three constant lengths ----
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '0.01');

Lvals = [0.20 0.19 0.18];
for k = 1:numel(Lvals)
    pc = sprintf('%s/P_c%d', mdl, k);
    lc = sprintf('%s/L_c%d', mdl, k);
    bb = sprintf('%s/b%d', mdl, k);
    lg = sprintf('%s/logF%d', mdl, k);
    add_block('simulink/Sources/Constant', pc, 'Value', '620', 'Position', [20 20+120*k 60 50+120*k]);
    add_block('simulink/Sources/Constant', lc, 'Value', num2str(Lvals(k)), 'Position', [20 80+120*k 60 110+120*k]);
    add_block('SNS_Library/BPA_10mm', bb, 'Position', [120 40+120*k 220 140+120*k]);
    add_block('simulink/Sinks/To Workspace', lg, 'VariableName', sprintf('F10_%d', k), 'SaveFormat', 'Timeseries', 'Position', [280 70+120*k 360 100+120*k]);
    add_line(mdl, sprintf('P_c%d/1', k), sprintf('b%d/1', k), 'autorouting', 'on');
    add_line(mdl, sprintf('L_c%d/1', k), sprintf('b%d/2', k), 'autorouting', 'on');
    add_line(mdl, sprintf('b%d/1', k), sprintf('logF%d/1', k), 'autorouting', 'on');
end
out = sim(mdl);
rest = 0.20; kmax = 0.165; KMAX = (rest-kmax)/rest;
for k = 1:numel(Lvals)
    Fb = out.(['F10_' num2str(k)]).Data(end);
    rel = (rest - Lvals(k))/rest/KMAX;
    Fcf = (0.568207874671*(exp(-4.25442545542*rel)-1) + (620/620)*exp(-0.55972777762*rel^2)) * ...
          (620*(0.4895*atan(0.03068*(rest-0.0075)*620)));
    Fref = festo4(10, rel, 620)*(620*(0.4895*atan(0.03068*(rest-0.0075)*620)));
    fprintf('BPA10 L=%.3f: block=%9.2f  closed-form=%9.2f  festo4ref=%9.2f  rel=%.3f\n', Lvals(k), Fb, Fcf, Fref);
end
close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end

%% ---- BioMuscle isometric + internal taps ----
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '0.01');
add_block('simulink/Sources/Constant', [mdl '/A_c'], 'Value', '1', 'Position', [20 20 60 50]);
add_block('simulink/Sources/Constant', [mdl '/L_c'], 'Value', '0.15', 'Position', [20 80 60 110]);
add_block('simulink/Sources/Constant', [mdl '/V_c'], 'Value', '0', 'Position', [20 140 60 170]);
add_block('SNS_Library/BioMuscle', [mdl '/bio'], 'Position', [120 40 260 180]);
% break the library link so internal taps can attach
set_param([mdl '/bio'], 'LinkStatus', 'inactive');
add_block('simulink/Sinks/To Workspace', [mdl '/logF'], 'VariableName', 'Fbio', 'SaveFormat', 'Timeseries', 'Position', [320 90 400 120]);
add_line(mdl, 'A_c/1', 'bio/1', 'autorouting', 'on');
add_line(mdl, 'L_c/1', 'bio/2', 'autorouting', 'on');
add_line(mdl, 'V_c/1', 'bio/3', 'autorouting', 'on');
add_line(mdl, 'bio/1', 'logF/1', 'autorouting', 'on');
% tap internal signals for logging
taps = {'lam', 'flExp', 'vSwitch', 'concA', 'pSat', 'pAdd', 'FmaxG', 'cosMul', 'lsub', 'lmax', 'vDiv', 'vNeg', 'bSum', 'bAb', 'concDiv'};
for k = 1:numel(taps)
    lbl = ['tap_' taps{k}];
    try
        add_block('simulink/Sinks/To Workspace', [mdl '/' lbl], 'VariableName', lbl, 'SaveFormat', 'Timeseries', 'Position', [500 40+60*k 580 70+60*k]);
        add_line(mdl, sprintf('bio/%s/1', taps{k}), [lbl '/1'], 'autorouting', 'on');
    catch ME
        fprintf('tap %s failed: %s\n', taps{k}, ME.message);
    end
end
out = sim(mdl);
fprintf('BioMuscle isometric: F = %.2f N (expect 1000)\n', out.Fbio.Data(end));
for k = 1:numel(taps)
    v = ['tap_' taps{k}];
    if isfield(out, v)
        fprintf('  %-8s = %.4g\n', taps{k}, out.(v).Data(end));
    end
end
close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
fprintf('DEBUG DONE\n');
