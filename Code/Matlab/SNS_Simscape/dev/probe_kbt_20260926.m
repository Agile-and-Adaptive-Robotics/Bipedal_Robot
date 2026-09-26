function probe_kbt_20260926()
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
sns = fileparts(fileparts(mfilename('fullpath')));
load_system(fullfile(sns, [mdl '.slx']));
ph = get_param([sub '/x04_02_KB_R_001_1_RIGID'], 'PortHandles');
ports = [ph.RConn ph.LConn];
fprintf('KB_R ports: %d\n', numel(ports));
for j = 1:numel(ports)
    ln = get_param(ports(j), 'Line');
    if ln < 0
        fprintf('port%d: free\n', j);
        continue;
    end
    phs = [];
    try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
    names = {};
    for q = 1:numel(phs)
        h = phs(q);
        if isscalar(h) && h > 0 && h ~= ports(j)
            try
                names{end+1} = get_param(get_param(h, 'Parent'), 'Name'); %#ok<SAGROW>
            catch
            end
        end
    end
    fprintf('port%d: %s\n', j, strjoin(names, ' | '));
end
close_system(mdl, 0);
end
