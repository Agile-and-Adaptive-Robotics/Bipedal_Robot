function goal2_probe4()
% Where does the Integrator STATE port land when combined with reset/IC ports?
mdl = 'goal2_probe4';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); load_system(mdl);
add_block('simulink/Continuous/Integrator', [mdl '/a'], ...
    'InitialCondition', '0', 'ExternalReset', 'rising', 'Position', [100 100 140 150]);
try
    set_param([mdl '/a'], 'ShowStatePort', 'on');
    fprintf('set_param ShowStatePort ok\n');
catch ME
    fprintf('set_param ShowStatePort FAILED: %s\n', ME.message);
end
ph = get_param([mdl '/a'], 'PortHandles');
fprintf('reset+state: In=%d Out=%d Reset=%d State=%d\n', numel(ph.Inport), ...
    numel(ph.Outport), numel(ph.Reset), numel(ph.State));

add_block('simulink/Continuous/Integrator', [mdl '/b'], ...
    'InitialConditionSource', 'external', 'ExternalReset', 'rising', 'Position', [100 200 140 250]);
set_param([mdl '/b'], 'ShowStatePort', 'on');
ph2 = get_param([mdl '/b'], 'PortHandles');
fprintf('extIC+reset+state: In=%d Out=%d Reset=%d State=%d\n', numel(ph2.Inport), ...
    numel(ph2.Outport), numel(ph2.Reset), numel(ph2.State));
% which numeric outport is the state port? connect a terminator and try
try
    add_block('simulink/Sinks/Terminator', [mdl '/t1'], 'Position', [200 210 220 230]);
    add_line(mdl, 'b/2', 't1/1');
    fprintf('state port of b IS outport 2\n');
catch
    try
        add_line(mdl, 'b/1', 't1/1');
        fprintf('state port of b is NOT 2; only outport 1 exists\n');
    catch ME2
        fprintf('confusing: %s\n', ME2.message);
    end
end
close_system(mdl, 0);
fprintf('PROBE4 DONE\n');
end
