% probe5_adjacency.m  (AH001_Simscape_Repair, 2026-11-01)
% Adjacency of the knee-region RIGID groups in ah001_full_import: for each
% group of interest, list every conserving connection (neighbor block path,
% through one hop). Also resolve both sides of the two Cartesian joints and of
% every Rigid Transform that directly couples the KT and KB groups (the
% smimport-added rigid connections).

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
logFile = fullfile(repairRoot, 'logs', 'probe5_adjacency_20261101.txt');
fid = fopen(logFile, 'w');

mdl = 'ah001_full_import';
load_system(fullfile(repairRoot, 'working', [mdl '.slx']));

blks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
ref = strings(numel(blks),1);
for k = 1:numel(blks)
    try
        ref(k) = strtrim(strrep(string(get_param(blks{k}, 'ReferenceBlock')), newline, " "));
    catch
        ref(k) = "";
    end
end

groups = ["x04_01_KT_R_001_1_RIGID", "x04_02_KB_R_001_1_RIGID", "x04_03_KT_L_001_1_RIGID", ...
          "x04_04_KB_L_001_1_RIGID", "x04_05_BL_001_1_RIGID", "x04_05_BL_001_2_RIGID", ...
          "x04_06_FL_001_1_RIGID", "x04_06_FL_001_2_RIGID", ...
          "x03_01_FB_R_001_1_RIGID", "x03_03_FB_L_001_1_RIGID", ...
          "x03_02_FH_R_001_1_RIGID", "x03_04_FH_L_001_1_RIGID", ...
          "x05_01_TI_R_001_1_RIGID", "x05_02_TI_001_1_RIGID", ...
          "x02_01_PE_001_1_RIGID"];

for g = groups
    b = find_system(mdl, 'Name', char(g), 'Type', 'Block');
    if isempty(b), p(fid, 'group %s NOT FOUND\n', g); continue; end
    b = b{1};
    p(fid, '\n[%s]\n', g);
    ph = get_param(b, 'PortHandles');
    cons = [ph.LConn ph.RConn];
    tags = [repmat("L",1,numel(ph.LConn)), repmat("R",1,numel(ph.RConn))];
    for ii = 1:numel(cons)
        l = get_param(cons(ii), 'Line');
        if l == -1
            p(fid, '  %s%d: <unconnected>\n', tags(ii), ii);
            continue;
        end
        src = get_param(l, 'SrcBlockHandle');
        dst = get_param(l, 'DstBlockHandle');
        others = [];
        for x = [src dst]
            if x ~= -1 && x ~= get_param(b, 'Handle')
                others(end+1) = x; %#ok<AGROW>
            end
        end
        for jj = 1:numel(others)
            p(fid, '  %s%d -> %s\n', tags(ii), ii, blkpath(others(jj)));
        end
    end
end

% Cartesian joints: full side resolution (print the transform chain)
p(fid, '\n=== Cartesian joint chains ===\n');
for nm = ["Cartesian", "Cartesian1"]
    b = find_system(mdl, 'Name', char(nm), 'Type', 'Block', 'ReferenceBlock', 'sm_lib/Joints/Cartesian Joint');
    if isempty(b), continue; end
    b = b{1};
    p(fid, '\n[%s @ %s]\n', nm, get_param(b, 'Parent'));
    for which = ["LConn" "RConn"]
        ph = get_param(b, 'PortHandles');
        ports = ph.(which);
        for ii = 1:numel(ports)
            l = get_param(ports(ii), 'Line');
            if l == -1, continue; end
            src = get_param(l, 'SrcBlockHandle');
            dst = get_param(l, 'DstBlockHandle');
            for x = [src dst]
                if x ~= -1 && x ~= get_param(b, 'Handle')
                    p(fid, '  %s: chain=%s\n', which, chain(x, get_param(b,'Handle')));
                end
            end
        end
    end
end

try, close_system(mdl, 0); catch, end
fclose(fid);

function s = blkpath(h)
    s = "<?>";
    try
        s = string(get_param(h, 'Name'));
        par = string(get_param(h, 'Parent'));
        s = s + " @" + extractAfter(par, string(get_param(0,'CurrentModel')) + "/");
    catch
    end
end

function s = chain(h, stopAt)
% printable chain: walk conserving connections up to 6 hops
    s = "";
    cur = h;
    for it = 1:6
        if cur == stopAt, s = s + " <loop>"; return; end
        nm = "?";
        try, nm = string(get_param(cur, 'Name')); catch, end
        s = s + nm;
        try
            r = strtrim(strrep(string(get_param(cur, 'ReferenceBlock')), newline, " "));
        catch
            r = "";
        end
        if contains(r, "Transform"), s = s + "->"; else, return; end
        ph = get_param(cur, 'PortHandles');
        allp = [ph.LConn ph.RConn];
        nxt = -1;
        for ii = 1:numel(allp)
            l = get_param(allp(ii), 'Line');
            if l == -1, continue; end
            src = get_param(l, 'SrcBlockHandle');
            dst = get_param(l, 'DstBlockHandle');
            for x = [src dst]
                if x ~= -1 && x ~= cur
                    nxt = x; break;
                end
            end
            if nxt ~= -1, break; end
        end
        if nxt == -1, s = s + "<dead>"; return; end
        cur = nxt;
    end
    s = s + "<deep>";
end

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
