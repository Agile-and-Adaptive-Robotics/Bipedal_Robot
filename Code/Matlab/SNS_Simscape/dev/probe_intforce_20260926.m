function probe_intforce_20260926()
load_system('sm_lib');
blk = 'sm_lib/Forces and Torques/Internal Force';
try
    dp = fieldnames(get_param(blk, 'DialogParameters'));
    for k = 1:numel(dp)
        v = '';
        try, v = get_param(blk, dp{k}); catch, end
        if ischar(v) && ~isempty(v), v = [' = ' v]; else, v = ''; end
        fprintf('IF.%s%s\n', dp{k}, v);
    end
catch ME
    fprintf('IF probe: %s\n', ME.message);
end
ph = get_param(blk, 'PortHandles');
fprintf('IF ports: RConn=%d LConn=%d In=%d Out=%d\n', numel(ph.RConn), numel(ph.LConn), numel(ph.In), numel(ph.Out));

ts = 'sm_lib/Frames and Transforms/Transform Sensor';
try
    load_system(ts);
    dp = fieldnames(get_param(ts, 'DialogParameters'));
    fprintf('TS params: %s\n', strjoin(dp', ', '));
    ph = get_param(ts, 'PortHandles');
    fprintf('TS ports: RConn=%d LConn=%d Out=%d\n', numel(ph.RConn), numel(ph.LConn), numel(ph.Out));
catch ME
    fprintf('TS probe: %s\n', ME.message);
end
% PS-Simulink / Simulink-PS converter paths (nesl_utility)
load_system('nesl_utility');
l = find_system('nesl_utility', 'SearchDepth', 1, 'Type', 'Block');
for k = 1:numel(l), fprintf('NESL: %s\n', l{k}); end
end
