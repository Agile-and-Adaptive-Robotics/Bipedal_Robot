% probe4_jointgraph.m  (AH001_Simscape_Repair, 2026-11-01)
% Dump the joint graph of ah001_full_import: for every joint/constraint block,
% resolve both conserving sides through Rigid Transform chains to the owning
% RIGID group (or World/other block), via block PortConnectivity.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'probe4_jointgraph_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_full_import';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));

blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try, ref(k) = strtrim(string(get_param(blks{k}, 'ReferenceBlock'))); catch, ref(k) = ""; end
end

isJoint = startsWith(ref, "sm_lib/Joints/");
isConstr = contains(ref, "sm_lib/Constraints/");
p(fid, '=== joints / constraints ===\n');
for k = find(isJoint | isConstr)
    b = blks{k};
    nm = get_param(b, 'Name');
    jtype = "joint";
    if isConstr(k), jtype = "constr"; end
    sideB = hop(neigh(b, 1));
    sideF = hop(neigh(b, 2));
    p(fid, '%-30s %-6s  B-> %-38s  F-> %s\n', nm, jtype, sideB, sideF);
end

p(fid, '\n=== totals ===\n');
p(fid, 'revolute=%d spherical=%d cartesian=%d cylindrical=%d planar=%d sixdof=%d constraint=%d\n', ...
    nnz(ref == "sm_lib/Joints/Revolute Joint"), nnz(ref == "sm_lib/Joints/Spherical Joint"), ...
    nnz(ref == "sm_lib/Joints/Cartesian Joint"), nnz(ref == "sm_lib/Joints/Cylindrical Joint"), ...
    nnz(ref == "sm_lib/Joints/Planar Joint"), nnz(ref == "sm_lib/Joints/6-DOF Joint"), nnz(isConstr));

try, close_system(mdl, 0); catch, end
fclose(fid);

function h = neigh(b, which)
% h = handle of the block connected to b's which-th conserving port
    h = -1;
    pc = get_param(b, 'PortConnectivity');
    c = 0;
    for j = 1:numel(pc)
        t = pc(j).Type;
        if strcmp(t, 'LConn') || strcmp(t, 'RConn')
            c = c + 1;
            if c == which
                h = pc(j).Block;   % -1 if unconnected
                return;
            end
        end
    end
end

function nm = hop(h0)
% walk from h0 across Rigid Transform chains to the first non-transform block
    if h0 == -1 || isempty(h0)
        nm = "<unconnected>";
        return;
    end
    prev = -1;
    cur = h0;
    for it = 1:20
        try
            r = string(get_param(cur, 'ReferenceBlock'));
        catch
            r = "";
        end
        if ~(contains(r, "Rigid") && contains(r, "Transform"))
            nm = string(get_param(cur, 'Name'));
            par = get_param(cur, 'Parent');
            if ~strcmp(par, cur)
                nm = nm + " @" + string(get_param(par, 'Name'));
            end
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
