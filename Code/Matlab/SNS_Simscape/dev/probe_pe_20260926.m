function probe_pe_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

d1 = find_system(sub, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
fprintf('depth-1 blocks in sub: %d\n', numel(d1));
rig = d1(contains(d1, '_RIGID'));
fprintf('RIGID bodies at depth 1: %d\n', numel(rig));

pe = [sub '/x02_01_PE_001_1_RIGID'];
fprintf('PE exists: %d\n', getSimulinkBlockHandle(pe) > 0);
ph = get_param(pe, 'PortHandles');
ports = [ph.RConn ph.LConn];
fprintf('PE ports: %d (RConn %d, LConn %d)\n', numel(ports), numel(ph.RConn), numel(ph.LConn));
for pp = 1:min(numel(ports), 8)
    ln = -1; try, ln = get_param(ports(pp), 'Line'); catch, end
    d1w = port_parent(ports(pp));
    fprintf('  port%d line=%d  parent=%s\n', pp, ln, d1w);
end
close_system(mdl, 0);
end

function s = port_parent(h)
s = '<n/a>';
try, s = get_param(get_param(h, 'Parent'), 'Name'); catch, end
end

%% appended: far_ends internals for PE ports 3..8
function probe2()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));
pe = [sub '/x02_01_PE_001_1_RIGID'];
ph = get_param(pe, 'PortHandles');
ports = [ph.RConn ph.LConn];
for pp = 3:8
    porth = ports(pp);
    ln = get_param(porth, 'Line');
    fprintf('port%d line=%.0f class', pp, ln);
    try, fprintf(' SrcPH=%.0f', get_param(ln, 'SrcPortHandle')); catch, fprintf(' SrcPH=ERR'); end
    try, fprintf(' DstPH=%.0f', get_param(ln, 'DstPortHandle')); catch, fprintf(' DstPH=ERR'); end
    try, fprintf(' BRH=%.0f', get_param(ln, 'BranchHandles')); catch, end
    fprintf('\n');
    % try the far end via line's ports
    for fld = ["SrcPortHandle", "DstPortHandle"]
        try
            h = get_param(ln, char(fld));
            if isscalar(h) && h > 0 && h ~= porth
                fprintf('   far via %s: parent=%s\n', char(fld), get_param(get_param(h,'Parent'),'Name'));
            end
        catch
        end
    end
end
close_system(mdl, 0);
end
