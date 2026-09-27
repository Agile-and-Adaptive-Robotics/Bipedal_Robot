function probe_foot_path_20260926()
% exhaustive weld-path BFS foot_R -> KB_R (joint blocks excluded from paths)
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
sns = fileparts(fileparts(mfilename('fullpath')));
load_system(fullfile(sns, [mdl '.slx']));

adj = containers.Map('KeyType', 'char', 'ValueType', 'any');
bodies = find_system(sub, 'SearchDepth', 1, 'Type', 'Block');
for k = 1:numel(bodies)
    bn = get_param(bodies{k}, 'Name');
    if ~endsWith(bn, '_RIGID'), continue; end
    ph = get_param(bodies{k}, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    for pp = 1:numel(ports)
        neigh = all_far_bodies(ports(pp), bn);
        for q = 1:numel(neigh)
            key = bn;
            if ~isKey(adj, key), adj(key) = {}; end
            lst = adj(key);
            if ~any(strcmp(lst, neigh{q}))
                lst{end+1} = neigh{q};
                adj(key) = lst;
            end
        end
    end
end
fn = keys(adj);
fprintf('rigid adjacency:\n');
for k = 1:numel(fn)
    lst = adj(fn{k});
    fprintf('  %s <-> %s\n', fn{k}, strjoin(lst, ' | '));
end
close_system(mdl, 0);
end

function nbrs = all_far_bodies(porth, selfName)
% walk the whole connection net: line Src/Dst, port Branch, line BranchHandles,
% and any branch segment's own BranchHandles
nbrs = {};
phs = [porth];
seenL = [];
seenP = [porth];
while ~isempty(phs)
    h = phs(1);
    phs(1) = [];
    if ~isscalar(h) || h <= 0, continue; end
    if any(seenP == h), continue; end
    seenP(end+1) = h; %#ok<AGROW>
    % is this port on a RIGID body?
    try
        pn = get_param(get_param(h, 'Parent'), 'Name');
        if endsWith(pn, '_RIGID') && ~strcmp(pn, selfName)
            if ~any(strcmp(nbrs, pn))
                nbrs{end+1} = pn; %#ok<AGROW>
            end
            continue;   % crossed into another body; do not walk past it
        end
    catch
    end
    % walk the port's/segment's line
    ln = -1;
    try, ln = get_param(h, 'Line'); catch, end
    if isscalar(ln) && ln > 0 && ~any(seenL == ln)
        seenL(end+1) = ln; %#ok<AGROW>
        phs = [phs, line_ports(ln)]; %#ok<AGROW>
    end
    % branch structures at port or line level
    for obj = {h, ln}
        o = obj{1};
        if ~isscalar(o) || o <= 0, continue; end
        br = -1;
        try, br = get_param(o, 'Branch'); catch, end
        if isscalar(br) && br > 0
            try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
        end
        try
            bb = get_param(o, 'BranchHandles');
            if ~isempty(bb), phs = [phs bb]; end %#ok<AGROW>
        catch
        end
    end
end
end

function phs = line_ports(ln)
phs = [];
try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
end
