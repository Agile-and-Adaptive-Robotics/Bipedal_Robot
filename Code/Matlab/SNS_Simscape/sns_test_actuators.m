%% sns_test_actuators.m — validate BPA_10/20/40mm + BioMuscle blocks against Ben's equations
%
% 1) BPA blocks: ramp actuator length at 620 kPa, compare block force with
%    festo4(dia, rel, P) * Fmax reference (Functions\festo4.m, on path).
% 2) BPA blocks: pressure sweep at fixed length, same reference.
% 3) BioMuscle: constant-input spot checks vs inline Thelen equations.
cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);
addpath(fullfile(cdto, '..', 'Functions'));   % festo4.m + FestoLookup.mat

sns_build_actuators;
load_system('SNS_Library');

mdl = 'actTest';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '1');

% length ramp: Rest -> 0.5*Rest over 1 s (rel: 0 -> 1 at 620 kPa)
add_block('simulink/Sources/Ramp', [mdl '/Lramp'], 'slope', '-0.5', 'start', '0', 'Position', [50 100 80 130]);
add_block('simulink/Sources/Constant', [mdl '/Rest_c'], 'Value', '0.20', 'Position', [50 160 80 190]);
add_block('simulink/Math Operations/Sum', [mdl '/Lsum'], 'Inputs', '+-', 'Position', [120 120 150 160]);
add_block('simulink/Sources/Constant', [mdl '/P620'], 'Value', '620', 'Position', [50 250 80 280]);
rest = 0.20; kmax = 0.165; tend = 0; fit = 0;
cfg = {'BPA_10mm', 10; 'BPA_20mm', 20; 'BPA_40mm', 40};
refFmax = [ ...
    620*(0.4895*atan(0.03068*(rest-0.0075)*620)); ...
    620*(1.4877*atan(0.0248*(rest-0.0075)*620)); ...
    6000];
y = 100;
for k = 1:3
    nm = cfg{k, 1};
    add_block(['SNS_Library/' nm], [mdl '/' nm], 'Position', [220 y 320 y+80]);
    % per-block Rest/Kmax via mask defaults 0.20/0.165 already
    add_line(mdl, 'Lsum/1', [nm '/2'], 'autorouting', 'on');
    add_line(mdl, 'P620/1', [nm '/1'], 'autorouting', 'on');
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' nm], 'VariableName', ['F_' nm], 'SaveFormat', 'Timeseries', 'Position', [360 y 430 y+30]);
    add_line(mdl, [nm '/1'], ['log_' nm '/1'], 'autorouting', 'on');
    y = y + 120;
end
add_line(mdl, 'Lramp/1', 'Lsum/1', 'autorouting', 'on');
add_line(mdl, 'Rest_c/1', 'Lsum/2', 'autorouting', 'on');
add_block('simulink/Sinks/To Workspace', [mdl '/logL'], 'VariableName', 'log_L', 'SaveFormat', 'Timeseries', 'Position', [360 60 430 90]);
add_line(mdl, 'Lsum/1', 'logL/1', 'autorouting', 'on');
out = sim(mdl);

L = out.log_L.Data;
tmaxErr = 0;
names = {'BPA_10mm', 'BPA_20mm', 'BPA_40mm'};
for k = 1:3
    Fblk = out.(['F_' names{k}]).Data;
    dia = cfg{k, 2};
    fmax = refFmax(k);
    KMAX = (rest - kmax)/rest;
    rel = (rest - L)/rest/KMAX;
    Fref = festo4(dia, rel, 620) * fmax;
    err = max(abs(Fblk - Fref));
    fprintf('%s: length-ramp max |err| = %.4g N  (Fmax ref %.1f N)\n', names{k}, err, fmax);
    tmaxErr = max(tmaxErr, err);
end

% pressure sweep at fixed rel = 0.5: P from 0 to 620 over 1 s
close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '1');
add_block('simulink/Sources/Ramp', [mdl '/Pramp'], 'slope', '620', 'Position', [50 100 80 130]);
add_block('simulink/Sources/Constant', [mdl '/Lc'], 'Value', '0.5*(0.20+0.165)', 'Position', [50 160 80 190]);
y = 100;
for k = 1:3
    nm = cfg{k, 1};
    add_block(['SNS_Library/' nm], [mdl '/' nm], 'Position', [220 y 320 y+80]);
    add_line(mdl, 'Pramp/1', [nm '/1'], 'autorouting', 'on');
    add_line(mdl, 'Lc/1', [nm '/2'], 'autorouting', 'on');
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' nm], 'VariableName', ['F_' nm], 'SaveFormat', 'Timeseries', 'Position', [360 y 430 y+30]);
    add_line(mdl, [nm '/1'], ['log_' nm '/1'], 'autorouting', 'on');
    y = y + 120;
end
add_block('simulink/Sinks/To Workspace', [mdl '/logP'], 'VariableName', 'log_P', 'SaveFormat', 'Timeseries', 'Position', [360 60 430 90]);
add_line(mdl, 'Pramp/1', 'logP/1', 'autorouting', 'on');
out = sim(mdl);
P = out.log_P.Data;
relMid = (rest - 0.5*(rest+kmax))/rest/((rest-kmax)/rest);
for k = 1:3
    Fblk = out.(['F_' names{k}]).Data;
    fmax = refFmax(k);
    Fref = festo4(cfg{k, 2}, relMid*ones(size(P)), P) * fmax;
    err = max(abs(Fblk - Fref));
    fprintf('%s: pressure-sweep max |err| = %.4g N\n', names{k}, err);
    tmaxErr = max(tmaxErr, err);
end

% ---- BioMuscle spot checks ----
close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode45', 'StopTime', '0.5');
add_block('simulink/Sources/Constant', [mdl '/A_c'], 'Value', '1', 'Position', [50 100 80 130]);
add_block('simulink/Sources/Constant', [mdl '/L_c'], 'Value', '0.15', 'Position', [50 160 80 190]);
add_block('simulink/Sources/Constant', [mdl '/V_c'], 'Value', '0', 'Position', [50 220 80 250]);
add_block('SNS_Library/BioMuscle', [mdl '/BioMuscle'], 'Position', [220 110 340 210], ...
    'Fmax', '1000', 'l_opt', '0.10', 'l_slack', '0.05', 'alpha', '0', 'v_max', '10');
add_line(mdl, 'A_c/1', 'BioMuscle/1', 'autorouting', 'on');
add_line(mdl, 'L_c/1', 'BioMuscle/2', 'autorouting', 'on');
add_line(mdl, 'V_c/1', 'BioMuscle/3', 'autorouting', 'on');
add_block('simulink/Sinks/To Workspace', [mdl '/logF'], 'VariableName', 'log_F', 'SaveFormat', 'Timeseries', 'Position', [400 120 470 150]);
add_line(mdl, 'BioMuscle/1', 'logF/1', 'autorouting', 'on');
out = sim(mdl);
Fbio = out.log_F.Data(end);
% reference at Lmt=0.15, l_slack=0.05 -> l_ce = 0.10 = l_opt, lambda = 1
% fl=1, fv=1, fpas=0 -> F = Fmax
fprintf('BioMuscle isometric (lambda=1, A=1): F = %.2f N (expect 1000)\n', Fbio);
bmaxErr = abs(Fbio - 1000);

% shortening check: v = -0.5 m/s (= -5 l_opt/s, below v_max=10)
set_param([mdl '/V_c'], 'Value', '-0.5');
out = sim(mdl);
Fshort = out.log_F.Data(end);
sHat = 5; A = 0.25; b = A*10;
fv = b*(1+A)/(sHat + b) - A;
Fref = 1000*fv;
fprintf('BioMuscle concentric (v=5 l_opt/s): F = %.2f N (expect %.2f)\n', Fshort, Fref);
bmaxErr = max(bmaxErr, abs(Fshort - Fref));

% passive stretch check: Lmt = 0.05 + 1.3*0.10 = 0.18 -> lambda = 1.3
set_param([mdl '/V_c'], 'Value', '0');
set_param([mdl '/L_c'], 'Value', '0.18');
set_param([mdl '/A_c'], 'Value', '0');
out = sim(mdl);
Fpas = out.log_F.Data(end);
lam = 1.3; kpe = 4; epas = 0.6;
Fref = 1000*(exp(kpe*(lam-1)/epas)-1)/(exp(kpe)-1);
fprintf('BioMuscle passive (lambda=1.3, A=0): F = %.2f N (expect %.2f)\n', Fpas, Fref);
bmaxErr = max(bmaxErr, abs(Fpas - Fref));

close_system(mdl, 0);
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
fprintf('DONE. Overall BPA max err %.4g N; BioMuscle max err %.4g N\n', tmaxErr, bmaxErr);
assert(tmaxErr < 1.0 && bmaxErr < 2.0, 'actuator validation FAILED');
fprintf('ACTUATOR VALIDATION PASSED\n');
