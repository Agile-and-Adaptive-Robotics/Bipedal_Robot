% probe3_audit.m  (AH001_Simscape_Repair, 2026-11-01)
% Deep audit of the copied imported model:
%   - Mechanism/Solver config, model solver
%   - per-joint: primitive targets (fixture state), connectivity (which RIGID group /
%     block each port leads to, walking through Rigid Transforms)
%   - per-Solid: type, file, mass/density params -> mass ledger
%   - foot / spine / lump block identities and params
%   - BPA channel + neuron + synapse payload parameters (evidence record)

repairRoot = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Neuromechanical_Models/AH001_Simscape_Repair';
workingDir = fullfile(repairRoot, 'working');
snsDir     = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape';
logFile    = fullfile(repairRoot, 'logs', 'probe3_audit_20261101.txt');
fid = fopen(logFile, 'w');
addpath(fullfile(repairRoot, 'lib'));   % private SNS_Library copy; avoids shadowing originals

mdl = 'mdl_humanoid_lower_ah001_imported';
load_system(fullfile(workingDir, [mdl '.slx']));
top = [mdl '/x09_BA_001_1'];

% enumerate once, filter by ReferenceBlock string (find_system value-pair on
% ReferenceBlock is unreliable here)
allBlks = find_system(mdl, 'LookUnderMasks','all','FollowLinks','on','Type','Block');
allRefs = strings(numel(allBlks),1);
for k = 1:numel(allBlks)
    try, allRefs(k) = strtrim(string(get_param(allBlks{k}, 'ReferenceBlock'))); catch, allRefs(k) = ""; end
end

%% Model-level configuration
p(fid, '=== configuration ===\n');
p(fid, 'solver: %s\n', get_param(mdl, 'Solver'));
p(fid, 'SolverType: %s\n', get_param(mdl, 'SolverType'));
mc = find_system(mdl,'LookUnderMasks','all','FollowLinks','on','ReferenceBlock','sm_lib/Utilities/Mechanism Configuration');
for k = 1:numel(mc)
    dp = get_param(mc{k}, 'DialogParameters');
    fn = fieldnames(dp);
    p(fid, 'Mechanism Configuration: %s\n', mc{k});
    for j = 1:numel(fn)
        try
            val = get_param(mc{k}, fn{j});
            if ischar(val) || isstring(val)
                p(fid, '   %-22s = %s\n', fn{j}, val);
            end
        catch
        end
    end
end

%% Joint blocks: primitive targets + connectivity
p(fid, '\n=== joints ===\n');
jointMask = allRefs == "sm_lib/Joints/Revolute Joint" | ...
            allRefs == "sm_lib/Joints/Spherical Joint" | ...
            allRefs == "sm_lib/Joints/6-DOF Joint";
jblks = allBlks(jointMask);
for k = 1:numel(jblks)
    b = jblks{k};
    p(fid, '\n[%s]\n', b);
    dp = get_param(b, 'DialogParameters');
    fn = fieldnames(dp);
    for j = 1:numel(fn)
        % only interesting params: targets, axes, limits
        if contains(fn{j}, 'Target', 'IgnoreCase', true) || contains(fn{j}, 'Axis') || contains(fn{j}, 'Limit') ...
           || contains(fn{j}, 'ValueSource') || contains(fn{j}, 'Priority') || contains(fn{j}, 'Motion')
            try
                val = get_param(b, fn{j});
                if (ischar(val) || isstring(val)) && ~isempty(strtrim(char(val)))
                    p(fid, '   %-34s = %s\n', fn{j}, val);
                end
            catch
            end
        end
    end
    % connectivity: B and F ports
    ph = get_param(b, 'PortHandles');
    cons = [ph.RConn ph.LConn];
    tags = {'B','F'};
    for c = 1:numel(cons)
        try
            l = get_param(cons(c), 'Line');
            if l == -1
                p(fid, '   port %s : <unconnected>\n', tags{min(c,2)});
                continue;
            end
            dst = get_param(l, 'DstBlockHandle');
            src = get_param(l, 'SrcBlockHandle');
            other = setdiff([dst src], cons(c));
            names = '';
            for o = 1:numel(other)
                ob = other(o);
                % walk up through Rigid Transform chains: name the first non-transform block
                nm = walkToBody(fid, ob);
                names = [names char(nm) '; ']; %#ok<AGROW>
            end
            p(fid, '   port %s -> %s\n', tags{min(c,2)}, names);
        catch ME
            p(fid, '   port %s : err %s\n', tags{min(c,2)}, ME.message);
        end
    end
end

%% Solid blocks ledger
p(fid, '\n=== solid ledger ===\n');
solidMask = startsWith(allRefs, "sm_lib/Body Elements/");
sblks = allBlks(solidMask);
totMass = 0; totKnown = 0;
for k = 1:numel(sblks)
    b = sblks{k};
    dp = get_param(b, 'DialogParameters');
    if isfield(dp, 'FileName')
        fname = get_param(b, 'FileName');
        dens  = '?'; mtype = 'file';
        try, dens = get_param(b, 'Density'); catch, end
        p(fid, '  FILE  %-70s dens=%s\n', b, dens);
    else
        % generic solid: type, dimensions, mass/density
        st = '?'; try, st = get_param(b, 'Shape'); catch, end %#ok<NASGU>
        x = '';
        try, x = ['dims=' get_param(b, 'Dimensions')]; catch, end
        d = ''; try, d = ['dens=' get_param(b, 'Density')]; catch, end
        m = ''; try, m = ['mass=' get_param(b, 'Mass')]; catch, end
        r = ''; try, r = ['rad=' get_param(b, 'Radius')]; catch, end
        h = ''; try, h = ['h=' get_param(b, 'Height')]; catch, end
        p(fid, '  GEN   %-70s %s %s %s %s %s\n', b, x, r, h, d, m);
    end
end

%% Named blocks of interest (foot/spine/lump/pinT)
p(fid, '\n=== named aux blocks ===\n');
interest = find_system(top, 'SearchDepth', 1, 'Type', 'Block');
for k = 1:numel(interest)
    nm = get_param(interest{k}, 'Name');
    if any(startsWith(nm, {'foot','spine','lump','pinT','World_pelvis','Spherical','knee','pelvis_free','N_','S_','EXT_','FLX_'}))
        bt = get_param(interest{k}, 'BlockType');
        ref = '';
        try, ref = get_param(interest{k}, 'ReferenceBlock'); catch, end
        p(fid, '  %-18s BlockType=%-16s ref=%s\n', nm, bt, ref);
    end
end

%% BPA / neural payload parameters
p(fid, '\n=== BPA + neural payload ===\n');
pay = find_system(top, 'SearchDepth', 1, 'Type', 'Block', 'BlockType', 'SubSystem');
for k = 1:numel(pay)
    nm = get_param(pay{k}, 'Name');
    if contains(nm, 'BPA') || startsWith(nm, 'N_') || startsWith(nm, 'S_')
        dp = get_param(pay{k}, 'MaskNames');
        if ~isempty(dp)
            vals = get_param(pay{k}, 'MaskValues');
            p(fid, '  %s:\n', nm);
            for j = 1:numel(dp)
                p(fid, '     %-16s = %s\n', dp{j}, char(vals{j}));
            end
        end
    end
end

%% Model workspace smiData masses by ID (instance-independent ledger)
p(fid, '\n=== smiData.Solid mass ledger (per unique part) ===\n');
hws = get_param(mdl, 'ModelWorkspace');
sm = evalin(hws, 'smiData');
for k = 1:numel(sm.Solid)
    s = sm.Solid(k);
    try
        p(fid, '  %2d  %-58s mass=%.9g\n', k, char(string(s.ID)), s.mass);
    catch
    end
end

try, close_system(mdl, 0); catch, end
p(fid, '\nAUDIT DONE.\n');
fclose(fid);

function nm = walkToBody(fid, h)
% Walk from block handle h across Rigid Transform chains to find the first
% "body-ish" block (Solid subsystem / Reference Frame owner). Returns its path.
nm = '?';
cur = h;
for it = 1:12
    bt = get_param(cur, 'BlockType');
    try
        ref = get_param(cur, 'ReferenceBlock');
    catch
        ref = '';
    end
    if isempty(ref) || ~contains(ref, 'Rigid')
        nm = [get_param(cur,'Name') ' (' bt ')'];
        return;
    end
    % move across the other conserving port of the transform
    ph = get_param(cur, 'PortHandles');
    cons = [ph.RConn ph.LConn];
    nxt = -1;
    for c = 1:numel(cons)
        try
            l = get_param(cons(c), 'Line');
            if l ~= -1
                dst = get_param(l, 'DstBlockHandle');
                src = get_param(l, 'SrcBlockHandle');
                other = setdiff([dst src], cur);
                if ~isempty(other)
                    nxt = other(1);
                end
            end
        catch
        end
    end
    if nxt == -1
        nm = [get_param(cur,'Name') ' (transform, dead end)'];
        return;
    end
    cur = nxt;
end
nm = [get_param(cur,'Name') ' (depth limit)'];
end

function p(fid, fmt, varargin)
    s = sprintf(fmt, varargin{:});
    fprintf('%s', s);
    fprintf(fid, '%s', s);
end
