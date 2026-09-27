function probe_conv_20260926()
mdl = 'probe_conv_scratch';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
add_block('nesl_utility/PS-Simulink Converter', [mdl '/C1'], 'Position', [50 50 100 80]);
ph = get_param([mdl '/C1'], 'PortHandles');
fn = fieldnames(ph);
for k = 1:numel(fn)
    for j = 1:numel(ph.(fn{k}))
        fprintf('C1 %s(%d)\n', fn{k}, j);
    end
end
add_block('simulink/Sinks/To Workspace', [mdl '/W1'], 'Position', [150 50 210 80]);
try
    add_line(mdl, 'C1/1', 'W1/1');
    fprintf('C1/1 -> W1/1 OK\n');
catch ME
    fprintf('C1/1 -> W1/1 FAILED: %s\n', ME.message(1:min(end, 90)));
end
% and the suspicious pattern: second identical pair
add_block('nesl_utility/PS-Simulink Converter', [mdl '/C2'], 'Position', [50 130 100 160]);
add_block('simulink/Sinks/To Workspace', [mdl '/W2'], 'Position', [150 130 210 160]);
try
    add_line(mdl, 'C2/1', 'W2/1');
    fprintf('C2/1 -> W2/1 OK\n');
catch ME
    fprintf('C2/1 -> W2/1 FAILED: %s\n', ME.message(1:min(end, 90)));
end
% string concat check
cn = 'C2';
fprintf('name built: %s\n', [cn 'w/1']);
close_system(mdl, 0);
end
