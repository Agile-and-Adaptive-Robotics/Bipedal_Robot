% probe4b_jointgraph.m  (AH001_Simscape_Repair, 2026-11-01)
% Joint graph of ah001_full_import with newline-normalized ReferenceBlocks.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'probe4b_jointgraph_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_full_import';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));

blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try
        r = string(get_param(blks{k}, 'ReferenceBlock'));
        r = strrep(r, newline, " ");
        r = strtrim(r);
        ref(k) = r;
    catch
        ref(k) = "";
    end
end

isJoint  = startsWith(ref, "sm_lib/Joints/");
isConstr = startsWith(ref, "sm_lib/Constraints/");
sel = find(isJoint | isConstr);
p(fid, 'selected joint/constraint blocks: %d\n', numel(sel));

for kk = 1:numel(sel)
    k = sel(kk);
    b = blks{k};
    nm = get_param(b, 'Name');
    jtype = "joint";
    if isConstr(k), jtype = "constr"; end
    sideB = hop(neigh(b, 1));
    sideF = hop(neigh(b, 2));
    p(fid, '%-32s %-6s  B-> %-44s  F-> %s\n', nm, jtype, char(sideB), char(sideF));
end

p(fid, '\n=== totals ===\n');
uj = unique(ref(sel));
for t = uj'
    p(fid, '  %2d x %s\n', nnz(ref(sel) == t), t);
end

try, close_system(mdl, 0); catch, end
fclose(fid);

function h = neigh(b, which)
    h = -1;
    pc = get_param(b, 'PortConnectivity');
    c = 0;
    for j = 1:numel(pc)
        t = pc(j).Type;
        if strcmp(t, 'LConn') || strcmp(t, 'RConn')
            c = c + 1;
            if c == which
                h = pc(j).Block;
                return;
            end
        end
    end
end

function nm = hop(h0)
    if h0 == -1 || isempty(h0)
        nm = "<unconnected>";
        return;
    end
    prev = -1;
    cur = h0;
    for it = 1:20
        try
            r = string(get_param(cur, 'ReferenceBlock'));
            r = strrep(r, newline, " ");
        catch
            r = "";
        end
        if ~(contains(r, "Rigid") && contains(r, "Transform"))
            nm = string(get_param(cur, 'Name'));
            return;
        end
        pc = get_param(cur, 'PortConnectivity');
        nxt = -1;
        for j = 1:numel(pc)
            t = pc(j).Type;
            if (strcmp(t, 'LConn') || strcmp(t, 'RConn')) && pc(j).Block ~= -1 ...
                    && pc(j).Block ~= prev && pc(j).Block ~= cur
                nxt = pc(j).Block;
                break;
            end
        end
        if nxt == -1
            nm = string(get_param(cur, 'Name')) + " <deadend>";
            return;
        end
        prev = cur;
        cur = nxt;
    end
    nm = "<depth>";
end

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
