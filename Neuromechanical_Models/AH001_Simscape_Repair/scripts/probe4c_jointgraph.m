% probe4c_jointgraph.m  (AH001_Simscape_Repair, 2026-11-01)
% Joint graph of ah001_full_import using PortHandles + Line handles.

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'probe4c_jointgraph_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_full_import';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));

blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try
        r = strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " ");
        ref(k) = strtrim(r);
    catch
        ref(k) = "";
    end
end

sel = find(startsWith(ref, "sm_lib/Joints/") | startsWith(ref, "sm_lib/Constraints/"));
p(fid, 'joints+constraints: %d\n', numel(sel));
for kk = 1:numel(sel)
    k = sel(kk);
    b = blks{k};
    nm = get_param(b, 'Name');
    jtype = "joint"; if startsWith(ref(k), "sm_lib/Constraints"), jtype = "constr"; end
    sideB = hopside(b, 'LConn');
    sideF = hopside(b, 'RConn');
    p(fid, '%-28s %-6s  B-> %-42s  F-> %s\n', nm, jtype, sideB, sideF);
end

try, close_system(mdl, 0); catch, end
fclose(fid);

function nm = hopside(b, which)
% resolve the conserving neighbor of block b on port group which, then walk
% through Rigid Transform chains to the owning non-transform block
    nm = "<none>";
    ph = get_param(b, 'PortHandles');
    ports = ph.(which);
    if isempty(ports), nm = "<no" + which + ">"; return; end
    h = -1;
    for ii = 1:numel(ports)
        l = get_param(ports(ii), 'Line');
        if l == -1, continue; end
        src = get_param(l, 'SrcBlockHandle');
        dst = get_param(l, 'DstBlockHandle');
        cand = [];
        for x = [src dst]
            if x ~= -1 && x ~= get_param(b, 'Handle')
                cand(end+1) = x; %#ok<AGROW>
            end
        end
        if ~isempty(cand)
            h = cand(1);
            break;
        end
    end
    if h == -1, nm = "<unconnected>"; return; end
    % walk transforms
    prev = -1; cur = h;
    for it = 1:20
        try
            r = strrep(string(get_param(cur, 'ReferenceBlock')), newline, " ");
        catch
            r = "";
        end
        if ~(contains(r, "Rigid") && contains(r, "Transform"))
            owner = string(get_param(cur, 'Name'));
            % climb to named subsystem if inside one
            par = string(get_param(cur, 'Parent'));
            nm = owner + " [" + par + "]";
            return;
        end
        php = get_param(cur, 'PortHandles');
        allp = [php.LConn php.RConn];
        nxt = -1;
        for ii = 1:numel(allp)
            l = get_param(allp(ii), 'Line');
            if l == -1, continue; end
            src = get_param(l, 'SrcBlockHandle');
            dst = get_param(l, 'DstBlockHandle');
            for x = [src dst]
                if x ~= -1 && x ~= prev && x ~= get_param(cur, 'Handle')
                    nxt = x;
                    break;
                end
            end
            if nxt ~= -1, break; end
        end
        if nxt == -1
            nm = string(get_param(cur, 'Name')) + " <deadend>";
            return;
        end
        prev = cur; cur = nxt;
    end
    nm = "<depth>";
end

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
