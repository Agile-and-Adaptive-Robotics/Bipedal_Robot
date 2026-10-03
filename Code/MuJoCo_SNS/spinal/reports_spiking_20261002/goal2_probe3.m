function goal2_probe3()
% Try all 6 permutations of connecting {xdot, IC, reset} to Integrator
% inports 1..3 and identify the correct one from the trajectory signature.
% Circuit: y' = -y, IC constant 5, reset = step 0->1 at t=0.2.
% Correct wiring => y(0)=5, y(0.15)~4.30, y(0.21)~4.95, y(0.35)~4.30.
sig = {'neg','ic','rst'};
P = perms(1:3);
mdl0 = 'goal2_probe3';
for p = 1:size(P, 1)
    mdl = sprintf('%s_%d', mdl0, p);
    if bdIsLoaded(mdl), close_system(mdl, 0); end
    new_system(mdl); load_system(mdl);
    add_block('simulink/Continuous/Integrator', [mdl '/y'], ...
        'InitialConditionSource', 'external', 'ExternalReset', 'rising', 'Position', [150 100 190 150]);
    add_block('simulink/Sources/Constant', [mdl '/ic'], 'Value', '5', 'Position', [40 160 80 190]);
    add_block('simulink/Sources/Step', [mdl '/rst'], 'Time', '0.2', 'Before', '0', 'After', '1', 'Position', [40 220 80 250]);
    add_block('simulink/Math Operations/Gain', [mdl '/neg'], 'Gain', '-1', 'Position', [250 110 280 140]);
    add_line(mdl, 'y/1', 'neg/1');
    for j = 1:3
        add_line(mdl, [sig{j} '/1'], sprintf('y/%d', P(p, j)));
    end
    set_param(mdl, 'StopTime', '0.4', 'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
    pho = get_param([mdl '/y'], 'PortHandles');
    set_param(pho.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', 'DataLoggingName', 'y');
    try
        out = sim(mdl);
        y = out.sigs.get('y').Values;
        fprintf('perm xdot=%d ic=%d rst=%d : y(0)=%7.3f y(0.15)=%7.3f y(0.21)=%7.3f y(0.35)=%7.3f\n', ...
            P(p,1), P(p,2), P(p,3), y.Data(1), ...
            interp1(y.Time, y.Data, 0.15), interp1(y.Time, y.Data, 0.21), ...
            interp1(y.Time, y.Data, 0.35));
    catch ME
        fprintf('perm xdot=%d ic=%d rst=%d : ERROR %s\n', P(p,1), P(p,2), P(p,3), ME.message);
    end
    close_system(mdl, 0);
end
fprintf('PROBE3 DONE (correct = ports xdot=1, ic=2, rst=3 if row 1 shows 5/4.3/4.95/4.3)\n');
end
