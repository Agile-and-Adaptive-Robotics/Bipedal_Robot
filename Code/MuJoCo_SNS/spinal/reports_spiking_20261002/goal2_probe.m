function goal2_probe()
% Quick Simulink probes for the spiking-block design:
%  P1: does Sum('+' ) with one port pass a vector through or sum elements?
%     (the NonSpikingNeuron Iapp front end claims element-summing)
%  P2: Integrator port layout with external initial condition + rising reset
%     (needed for the SpikingSynapse conductance state g)
rep = fileparts(mfilename('fullpath'));
mdl = 'goal2_probe';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); load_system(mdl);

add_block('simulink/Sources/Constant', [mdl '/cv'], 'Value', '[1 2 3]', 'Position', [40 30 80 60]);
add_block('simulink/Math Operations/Sum', [mdl '/s1'], 'Inputs', '+', 'Position', [120 30 150 60]);
add_block('simulink/Sinks/Out1', [mdl '/o1'], 'Port', '1', 'Position', [200 30 230 44]);
add_line(mdl, 'cv/1', 's1/1');
add_line(mdl, 's1/1', 'o1/1');
ph = get_param([mdl '/s1'], 'PortHandles');
fprintf('P1: Sum ''+'' port handles: in=%d out=%d\n', numel(ph.Inport), numel(ph.Outport));

add_block('simulink/Continuous/Integrator', [mdl '/integ'], ...
    'InitialConditionSource', 'external', 'ExternalReset', 'rising', ...
    'Position', [120 120 160 170]);
ph2 = get_param([mdl '/integ'], 'PortHandles');
fn = fieldnames(ph2);
fprintf('P2: Integrator(extIC,rising) PortHandles fields:\n');
for k = 1:numel(fn)
    fprintf('   %-18s n=%d\n', fn{k}, numel(ph2.(fn{k})));
end
close_system(mdl, 0);
fprintf('PROBE DONE\n');
end
