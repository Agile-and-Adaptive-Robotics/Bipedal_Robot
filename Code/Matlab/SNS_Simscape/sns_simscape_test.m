%% sns_simscape_test.m — prove the Simscape SNS blocks: RC membrane + E/I synapses
%
% Network (one electrical network, one shared reference rail):
%   DC Current Source (4 nA) -> PRE membrane
%   PRE.m --[NonSpikingSynapse Esyn=  0 V ]--> EXC.m   (excitatory)
%   PRE.m --[NonSpikingSynapse Esyn=-72 mV]--> INH.m   (inhibitory)
%
% Analytics (steady state):
%   PRE: Vrest + I/Gm = -52 mV + 4nA/0.4uS = -42 mV  (Sat drive = 1, fully on)
%   EXC: V = (Gm*Vrest + gmax*Esyn)/(Gm+gmax) with Esyn=0    -> ~ -26 mV (depolarizes)
%   INH: same with Esyn=-72 mV                                -> ~ -62 mV (hyperpolarizes)
% The depol/hyperpol pair IS the E/I distinction.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
cd(thisDir);
addpath(fullfile(thisDir, 'simscape_sources'));   % where SNS_lib.slx lives
load_system('SNS_lib');

mdl = 'sns_simscape_EI_test';
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end
new_system(mdl); load_system(mdl);
set_param(mdl, 'Solver', 'ode23t', 'StopTime', '0.3');

%% blocks
add_block('nesl_utility/Solver Configuration', [mdl '/SolverCfg'], 'Position', [40 320 100 360]);
add_block('fl_lib/Electrical/Electrical Elements/Electrical Reference', [mdl '/GND'], 'Position', [40 380 70 410]);
add_block('fl_lib/Electrical/Electrical Sources/DC Current Source', [mdl '/Stim'], 'Position', [40 60 90 110], 'i0', '4e-9');

% resolve library block paths at runtime (ssc_build names may wrap)
libBlocks = find_system('SNS_lib', 'SearchDepth', 1, 'Type', 'Block');
libPath = containers.Map();
for h = 1:numel(libBlocks)
    nm = strrep(strrep(get_param(libBlocks{h}, 'Name'), newline, ' '), ' ', '');
    if contains(nm, 'NonSpikingNeuron'),     libPath('neuron')  = libBlocks{h}; end
    if contains(nm, 'NonSpikingSynapse'),    libPath('synapse') = libBlocks{h}; end
end
disp(['lib neuron block:  ' strrep(libPath('neuron'), newline, '|')]);
disp(['lib synapse block: ' strrep(libPath('synapse'), newline, '|')]);

add_block(libPath('neuron'), [mdl '/PRE'], 'Position', [200 40 300 140], 'Gm', '4e-7');
add_block(libPath('neuron'), [mdl '/EXC'], 'Position', [470 40 570 140]);
add_block(libPath('neuron'), [mdl '/INH'], 'Position', [470 240 570 340]);

add_block(libPath('synapse'), [mdl '/SynEXC'], 'Position', [360 60 410 110], 'Esyn', '0',      'gmax', '1e-7');
add_block(libPath('synapse'), [mdl '/SynINH'], 'Position', [360 260 410 310], 'Esyn', '-0.072', 'gmax', '1e-7');

%% physical wiring by port handles (order = node declaration order m,r / p,n)
ph  = struct();
for nm = {'PRE','EXC','INH','SynEXC','SynINH','Stim','GND','SolverCfg'}
    ph.(nm{1}) = get_param([mdl '/' nm{1}], 'PortHandles');
    fprintf('%-8s LConn=%d RConn=%d Out=%d\n', nm{1}, numel(ph.(nm{1}).LConn), numel(ph.(nm{1}).RConn), numel(ph.(nm{1}).Outport));
end

% neuron conserving ports: m = first conn port, r = second (check both sides)
neuronConns = struct();
for nm = {'PRE','EXC','INH'}
    h = ph.(nm{1});
    neuronConns.(nm{1}) = [h.LConn h.RConn];   % declaration order
end
synConns = struct();
for nm = {'SynEXC','SynINH'}
    h = ph.(nm{1});
    synConns.(nm{1}) = [h.LConn h.RConn];      % [p, n]
end
stimConns = [ph.Stim.LConn ph.Stim.RConn];     % [+, -]
gndConn   = [ph.GND.LConn ph.GND.RConn];       % [V]

% tie all r ports to the shared rail (they'll be connected through GND node)
% Strategy: single Electrical Reference connects to PRE.r; r's of all neurons tie together.
add_line(mdl, stimConns(1), neuronConns.PRE(1), 'autorouting', 'on');   % Stim + -> PRE.m
add_line(mdl, stimConns(2), neuronConns.PRE(2),  'autorouting', 'on');  % Stim - -> PRE.r
add_line(mdl, neuronConns.EXC(2), neuronConns.PRE(2), 'autorouting', 'on');  % EXC.r == rail
add_line(mdl, neuronConns.INH(2), neuronConns.PRE(2), 'autorouting', 'on');  % INH.r == rail
add_line(mdl, gndConn(1), neuronConns.PRE(2), 'autorouting', 'on');     % GND on rail
add_line(mdl, ph.SolverCfg.RConn(1), neuronConns.PRE(2), 'autorouting', 'on'); % solver on network
% synapses: p senses PRE.m, n injects into post membrane
add_line(mdl, synConns.SynEXC(1), neuronConns.PRE(1), 'autorouting', 'on');
add_line(mdl, synConns.SynEXC(2), neuronConns.EXC(1), 'autorouting', 'on');
add_line(mdl, synConns.SynINH(1), neuronConns.PRE(1), 'autorouting', 'on');
add_line(mdl, synConns.SynINH(2), neuronConns.INH(1), 'autorouting', 'on');

save_system(mdl);

%% voltage sensors -> PS-Simulink converters -> To Workspace (bulletproof logging)
targets = {'PRE', 'EXC', 'INH'};
y0 = [200 470 470];
for t = 1:numel(targets)
    nm = targets{t};
    add_block('fl_lib/Electrical/Electrical Sensors/Voltage Sensor', [mdl '/Vs_' nm], 'Position', [180 y0(t)+120 230 y0(t)+160]);
    add_block('nesl_utility/PS-Simulink Converter', [mdl '/ps_' nm], 'Position', [280 y0(t)+125 330 y0(t)+155]);
    add_block('simulink/Sinks/To Workspace', [mdl '/log_' nm], 'VariableName', ['log_' nm], 'SaveFormat', 'Timeseries', 'Position', [370 y0(t)+128 440 y0(t)+152]);
    phv = get_param([mdl '/Vs_' nm], 'PortHandles');
    php = get_param([mdl '/ps_' nm], 'PortHandles');
    disp(['Vs_' nm ': sensor L=' num2str(numel(phv.LConn)) ' R=' num2str(numel(phv.RConn)) ...
        ' | conv L=' num2str(numel(php.LConn)) ' R=' num2str(numel(php.RConn)) ' O=' num2str(numel(php.Outport))]);
    add_line(mdl, neuronConns.(nm)(1), phv.LConn(1), 'autorouting', 'on');   % + senses membrane
    add_line(mdl, neuronConns.(nm)(2), phv.RConn(2), 'autorouting', 'on');   % - on reference rail (RConn2 = n; RConn1 is the PS output)
    add_line(mdl, phv.RConn(1), php.LConn(1), 'autorouting', 'on');          % PS V out -> converter
    add_line(mdl, php.Outport(1), get_param([mdl '/log_' nm], 'PortHandles').Inport(1), 'autorouting', 'on');
end

%% simulate
out = sim(mdl);
vPre = out.log_PRE.Data(end)*1e3; vExc = out.log_EXC.Data(end)*1e3; vInh = out.log_INH.Data(end)*1e3;
fprintf('\nsteady-state membrane voltages:\n');
fprintf('  PRE (stimulated) : %7.2f mV   (analytic -42.0)\n', vPre);
fprintf('  EXC (Esyn=  0 mV): %7.2f mV   (analytic  -26.0)  -> depolarized\n', vExc);
fprintf('  INH (Esyn=-72 mV): %7.2f mV   (analytic  -62.0)  -> hyperpolarized\n', vInh);

% polarity assertions vs rest (-52 mV)
depol   = vExc > -45;
hyper   = vInh < -58;
correct = depol && hyper;
fprintf('E/I polarity check: EXC depolarizes=%d, INH hyperpolarizes=%d -> %s\n', ...
    depol, hyper, ternary(correct, 'PASS', 'FAIL'));

fig = figure('Visible', 'off', 'Position', [80 80 900 420]);
plot(out.log_PRE.Time, out.log_PRE.Data*1e3, 'LineWidth', 1.4); hold on;
plot(out.log_EXC.Time, out.log_EXC.Data*1e3, 'LineWidth', 1.4);
plot(out.log_INH.Time, out.log_INH.Data*1e3, 'LineWidth', 1.4);
yline(-52, '--', 'V_{rest}', 'HandleVisibility', 'off');
grid on; xlabel('time (s)'); ylabel('membrane potential (mV)');
legend('PRE (4 nA stim)', 'postsyn EXC (Esyn = 0)', 'postsyn INH (Esyn = -72 mV)', 'Location', 'southeast');
title('Simscape SNS proof: RC-membrane neurons + E/I synapses');
exportgraphics(fig, fullfile(thisDir, 'sns_simscape_EI_test.png'), 'Resolution', 150);
fprintf('plot saved: sns_simscape_EI_test.png\n');

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
