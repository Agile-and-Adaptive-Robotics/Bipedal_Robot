function goal2_probe2()
% Determine the Integrator port ORDER with external IC + rising reset.
% Circuit: y' = -y  (decay); external IC constant 5; reset = step 0->1 at t=0.2.
% Hypothesis A (xdot=1, IC=2, reset=3): y(0)=5, decays toward 0, at t=0.2
%   jumps to IC-port value (5) then decays again.
% If IC/reset are swapped, the sim either errors (dim) or y resets to the
% step value (1) -> y jumps to 1.
mdl = 'goal2_probe2';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl); load_system(mdl);
add_block('simulink/Continuous/Integrator', [mdl '/y'], ...
    'InitialConditionSource', 'external', 'ExternalReset', 'rising', 'Position', [150 100 190 150]);
add_block('simulink/Sources/Constant', [mdl '/ic'], 'Value', '5', 'Position', [40 160 80 190]);
add_block('simulink/Sources/Step', [mdl '/rst'], 'Time', '0.2', 'Before', '0', 'After', '1', 'Position', [40 220 80 250]);
add_block('simulink/Math Operations/Gain', [mdl '/neg'], 'Gain', '-1', 'Position', [250 110 280 140]);
add_line(mdl, 'y/1', 'neg/1');
add_line(mdl, 'neg/1', 'y/1');          % xdot  -> Inport 1 (hypothesis A)
add_line(mdl, 'ic/1', 'y/2');           % IC     -> Inport 2
add_line(mdl, 'rst/1', 'y/3');          % reset  -> Inport 3
ph = get_param([mdl '/y'], 'PortConnectivity');
for k = 1:numel(ph)
    fprintf('port %d: type=%s\n', ph(k).Type, num2str(ph(k).Type));
end
set_param(mdl, 'StopTime', '0.4', 'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
pho = get_param([mdl '/y'], 'PortHandles');
set_param(pho.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'y');
out = sim(mdl);
y = out.sigs.get('y').Values;
fprintf('y(0)=%.3f  y(0.15)=%.3f  y(0.21)=%.3f  y(0.35)=%.3f\n', ...
    y.Data(1), interp1(y.Time, y.Data, 0.15), interp1(y.Time, y.Data, 0.21), interp1(y.Time, y.Data, 0.35));
close_system(mdl, 0);
fprintf('PROBE2 DONE\n');
end
