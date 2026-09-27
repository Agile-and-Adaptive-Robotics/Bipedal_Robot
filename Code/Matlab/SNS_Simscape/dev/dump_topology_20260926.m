function dump_topology_20260926()
%% dump_topology_20260926  Read-only topology recon of the two imported models:
%% every joint/constraint block with its two neighbours (via Simscape line
%% handles), plus the smiData mass inventory. Feeds the knee/pelvis surgery.

sns = fileparts(fileparts(mfilename('fullpath')));

models = {'mdl_leg_rig_ba003_imported', 'mdl_humanoid_lower_ah001_imported'};
datafiles = {'mdl_leg_rig_ba003_DataFile.m', 'mdl_humanoid_lower_ah001_DataFile.m'};

for mi = 1:numel(models)
    mdl = models{mi};
    fprintf('\n############ %s ############\n', mdl);

    % ---- mass inventory ----
    smiData = [];
    run(fullfile(sns, datafiles{mi}));   % defines smiData in this function ws
    fprintf('--- smiData.Solid inventory ---\n');
    tot = 0;
    for k = 1:numel(smiData.Solid)
        fprintf('  Solid(%2d) %10.4f kg  %s\n', k, smiData.Solid(k).mass, smiData.Solid(k).ID);
        tot = tot + smiData.Solid(k).mass;
    end
    fprintf('  SUM (per unique part, not per instance): %.4f kg\n', tot);

    % ---- connectivity ----
    fprintf('--- JOINT CONNECTIVITY (via line handles) ---\n');
    load_system(fullfile(sns, [mdl '.slx']));
    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    for k = 1:numel(blks)
        rb = '';
        try, rb = get_param(blks{k}, 'ReferenceBlock'); if ~ischar(rb), rb = ''; end, catch, end
        isJoint = ~isempty(strfind(rb, 'sm_lib/Joints/'));
        isWeld  = ~isempty(strfind(rb, 'sm_lib/Constraints/'));
        if ~isJoint && ~isWeld, continue; end
        nm = strrep(blks{k}, [mdl '/'], '');
        ph = get_param(blks{k}, 'PortHandles');
        ports = [ph.RConn ph.LConn];
        nbr = {};
        for pp = 1:numel(ports)
            others = neighbors_of_port(ports(pp));
            if isempty(others)
                nbr{end+1} = '<unconnected>'; %#ok<SAGROW>
            else
                lbl = '';
                for o = 1:numel(others)
                    on = strrep(get_param(others(o), 'Name'), [mdl '/'], '');
                    if strcmp(on, nm), continue; end
                    lbl = [lbl strrep(on, [mdl '/'], '') ' ']; %#ok<AGROW>
                end
                if isempty(lbl), lbl = '<self>'; end
                nbr{end+1} = lbl; %#ok<SAGROW>
            end
        end
        tag = 'JOINT';
        if isWeld, tag = 'WELD/CON'; end
        fprintf('%s %-42s : %s -- %s\n', tag, nm, nbr{1}, nbr{min(2, numel(nbr))});
    end
    close_system(mdl, 0);
end
fprintf('\n=== dump_topology DONE ===\n');
end

function out = neighbors_of_port(porth)
% all blocks at the far end of the Simscape connection attached to porth
out = [];
ln = -1;
try, ln = get_param(porth, 'Line'); catch, end
if ln < 0, return; end
phs = [];
try
    phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')];
catch
end
% branched connections: also walk branch handles
br = -1;
try, br = get_param(porth, 'Branch'); catch, end
if isscalar(br) && br > 0
    try
        phs = [phs, get_param(br, 'BranchHandles')]; %#ok<AGROW>
    catch
    end
end
for q = 1:numel(phs)
    h = phs(q);
    if h <= 0 || h == porth, continue; end
    try
        b = get_param(h, 'Parent');
        b = get_param(b, 'Handle');
    catch
        continue;
    end
    if isempty(out) || ~any(out == b)
        out(end+1) = b; %#ok<AGROW>
    end
end
end
