function probe_rigid_path_20260926()
% BFS over direct (weld) lines between RIGID bodies inside the sub
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
sns = fileparts(fileparts(mfilename('fullpath')));
load_system(fullfile(sns, [mdl '.slx']));

start = [sub '/x04_01_KT_R_001_1_RIGID'];
target = [sub '/x04_02_KB_R_001_1_RIGID'];

adj = containers.Map('KeyType', 'char', 'ValueType', 'any');
bodies = find_system(sub, 'SearchDepth', 1, 'Type', 'Block');
for k = 1:numel(bodies)
    bn = get_param(bodies{k}, 'Name');
    if ~endsWith(bn, '_RIGID'), continue; end
    ph = get_param(bodies{k}, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    for pp = 1:numel(ports)
        ln = get_param(ports(pp), 'Line');
        if ln < 0, continue; end
        phs = [];
        try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
        br = -1;
        try, br = get_param(ports(pp), 'Branch'); catch, end
        if isscalar(br) && br > 0
            try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
        end
        for q = 1:numel(phs)
            h = phs(q);
            if ~isscalar(h) || h <= 0 || h == ports(pp), continue; end
            try
                pn = get_param(get_param(h, 'Parent'), 'Name');
            catch
                continue;
            end
            if endsWith(pn, '_RIGID') && ~strcmp(pn, bn)
                key = bn;
                if ~isKey(adj, key), adj(key) = {}; end
                lst = adj(key);
                if ~any(strcmp(lst, pn))
                    lst{end+1} = pn;
                    adj(key) = lst;
                end
            end
        end
    end
end

% BFS
queue = {start(1:strlength(start)-6)};  % strip '_RIGID'
startN = queue{1};
visited = {startN};
prev = containers.Map('KeyType', 'char', 'ValueType', 'char');
found = false;
while ~isempty(queue)
    cur = queue{1};
    queue(1) = [];
    if strcmp(cur, target(1:strlength(target)-6))
        found = true;
        break;
    end
    if ~isKey(adj, cur), continue; end
    nbrs = adj(cur);
    for k = 1:numel(nbrs)
        if ~any(strcmp(visited, nbrs{k}))
            visited{end+1} = nbrs{k}; %#ok<SAGROW>
            prev(nbrs{k}) = cur;
            queue{end+1} = nbrs{k}; %#ok<SAGROW>
        end
    end
end
fprintf('rigid path KT_R -> KB_R exists: %d\n', found);
if found
    path = target(1:strlength(target)-6);
    cur = path;
    while ~strcmp(cur, startN)
        cur = prev(cur);
        path = [cur ' -> ' path]; %#ok<AGROW>
    end
    fprintf('%s\n', path);
end
% also list KB_R weld neighbours
if isKey(adj, target(1:strlength(target)-6))
    fprintf('KB_R weld neighbours: %s\n', strjoin(adj(target(1:strlength(target)-6)), ' | '));
end
close_system(mdl, 0);
end
