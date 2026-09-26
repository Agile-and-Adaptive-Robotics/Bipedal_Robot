function probe_if_ports_20260926()
% discovery: Internal Force port roles + Transform Sensor output order
mdl = 'probe_if_scratch';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
load_system('sm_lib');

add_block('sm_lib/Forces and Torques/Internal Force', [mdl '/IF'], 'Position', [100 100 160 160]);
add_block('sm_lib/Frames and Transforms/World Frame', [mdl '/W1'], 'Position', [30 80 60 110]);
add_block('nesl_utility/Solver Configuration', [mdl '/Sol'], 'Position', [30 240 90 290]);

ph = get_param([mdl '/IF'], 'PortHandles');
fn = fieldnames(ph);
for k = 1:numel(fn)
    for j = 1:numel(ph.(fn{k}))
        fprintf('IF %s(%d)\n', fn{k}, j);
    end
end

wph = get_param([mdl '/W1'], 'PortHandles');
wfn = fieldnames(wph);
for k = 1:numel(wfn)
    fprintf('W1 %s: %d\n', wfn{k}, numel(wph.(wfn{k})));
end
wport = [wph.RConn wph.LConn];
ifh = get_param([mdl '/IF'], 'PortHandles');
iports = [ifh.LConn ifh.RConn];
for j = 1:numel(iports)
    try
        add_line(mdl, wport(1), iports(j));
        fprintf('WORLD -> IF port %d OK (frame port)\n', j);
    catch ME
        fprintf('WORLD -> IF port %d FAILED (%s)\n', j, ME.message(1:min(end, 80)));
    end
end
close_system(mdl, 0);
end
