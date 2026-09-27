function recon_patch_20260926()
%% recon_patch_20260926  Patch the 2026-09-26 smimport results:
%%   1. gravity -> [0 -9.80665 0] (CAD convention is +y up, NOT +z)
%%   2. File Solid geometry -> absolute STEP paths (map from the source XMLs)
%%   3. dump joint<->body connectivity (for the knee/pelvis rewire decision)
%%   4. compile check; save in every case so the path patch survives.

here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
addpath(sns);
rawRig = jsondecode(fileread(fullfile(here, 'instance_step_map_rig.json')));
rawHum = jsondecode(fileread(fullfile(here, 'instance_step_map_humanoid.json')));

models = {'mdl_leg_rig_ba003_imported', 'mdl_humanoid_lower_ah001_imported'};
raws = {rawRig, rawHum};

for mi = 1:numel(models)
    mdl = models{mi};
    fprintf('\n############ %s ############\n', mdl);
    % normalized block-parent-name -> STEP path (parallel arrays; jsondecode
    % mangles non-identifier object keys, so never key by instance name)
    raw = raws{mi};
    lut = containers.Map('KeyType', 'char', 'ValueType', 'char');
    for k = 1:numel(raw.instances)
        norm = regexprep(regexprep(raw.instances{k}, '[ \-\.]', '_'), '_+', '_');
        lut(['x' norm '_RIGID']) = char(raw.steps{k});   % digit-leading names
        lut([norm '_RIGID'])     = char(raw.steps{k});   % letter-leading names (no x prefix)
    end

    load_system(fullfile(sns, [mdl '.slx']));

    % ---- 1. gravity: CAD +y up ----
    mc = find_system(mdl, 'LookUnderMasks', 'all', 'MaskType', 'Mechanism Configuration');
    if isempty(mc)
        fprintf('NO Mechanism Configuration block\n');
    else
        for g = 1:numel(mc)
            set_param(mc{g}, 'GravityVector', '[0 -9.80665 0]');
        end
        fprintf('gravity set to [0 -9.80665 0] on %d Mechanism Configuration block(s)\n', numel(mc));
    end

    % ---- 2. File Solid paths ----
    fs = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'MaskType', 'File Solid');
    fprintf('%d File Solid blocks\n', numel(fs));
    % learn the parameter schema from the first block
    dp = fieldnames(get_param(fs{1}, 'DialogParameters'));
    fprintf('File Solid dialog params: %s\n', strjoin(dp', ', '));
    for p = 1:numel(dp)
        v = '';
        try, v = get_param(fs{1}, dp{p}); catch, end
        if ischar(v), fprintf('  %-24s = %s\n', dp{p}, v); end
    end

    nOK = 0; nMiss = 0; missNames = {};
    for f = 1:numel(fs)
        parent = regexprep(strrep(fs{f}, [mdl '/'], ''), '/Solid$', '');
        [~, base] = fileparts(parent);
        if lut.isKey(base)
            pth = lut(base);
        else
            % sub-assembly nesting: take everything after the last '/' segment set
            toks = split(base, '/');
            base2 = toks{end};
            if lut.isKey(base2), pth = lut(base2); base = base2;
            else
                nMiss = nMiss + 1; missNames{end+1} = base; %#ok<SAGROW>
                continue;
            end
        end
        nSet = 0;
        for p = 1:numel(dp)
            v = '';
            try, v = get_param(fs{f}, dp{p}); catch, end
            if ischar(v) && ~isempty(v) && ~isempty(regexpi(v, '\.(step|stl)$|package://', 'once'))
                set_param(fs{f}, dp{p}, pth);
                nSet = nSet + 1;
            end
        end
        if nSet == 0
            % maybe the file param is EMPTY rather than a wrong path: set the
            % params whose name suggests a file
            for p = 1:numel(dp)
                if ~isempty(regexpi(dp{p}, 'file', 'once'))
                    try, set_param(fs{f}, dp{p}, pth); nSet = nSet + 1; catch, end
                end
            end
        end
        if nSet > 0, nOK = nOK + 1; else, nMiss = nMiss + 1; missNames{end+1} = base; end %#ok<SAGROW>
    end
    fprintf('File Solids patched: %d OK, %d unresolved\n', nOK, nMiss);
    if ~isempty(missNames)
        fprintf('unresolved (first 10): %s\n', strjoin(missNames(1:min(10, numel(missNames))), ' | '));
    end

    % ---- 3. connectivity dump ----
    blks = find_system(mdl, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
    fprintf('\n--- JOINT CONNECTIVITY ---\n');
    for k = 1:numel(blks)
        rb = '';
        try, rb = get_param(blks{k}, 'ReferenceBlock'); if ~ischar(rb), rb = ''; end, catch, end
        if isempty(strfind(rb, 'sm_lib/Joints/')), continue; end
        nm = strrep(blks{k}, [mdl '/'], '');
        try
            pc = get_param(blks{k}, 'PortConnectivity');
            nbrs = {};
            for pp = 1:numel(pc)
                h = [];
                try, h = pc(pp).SrcBlockHandle; catch, end
                if ~isscalar(h) || h <= 0
                    try, h = pc(pp).DstBlockHandle; catch, end
                end
                if isscalar(h) && h > 0
                    nbrs{end+1} = strrep(get_param(h, 'Name'), [mdl '/'], ''); %#ok<SAGROW>
                else
                    nbrs{end+1} = '<unconnected>'; %#ok<SAGROW>
                end
            end
            fprintf('JOINT %-46s : %s\n', nm, strjoin(nbrs, '  --  '));
        catch ME
            fprintf('JOINT %-46s : <pc failed: %s>\n', nm, ME.message);
        end
    end

    % ---- 4. compile + save ----
    save_system(mdl);   % keep the path patch regardless
    try
        set_param(mdl, 'StopTime', '0.01');
        set_param(mdl, 'SimulationCommand', 'update');
        fprintf('%s: UPDATE OK\n', mdl);
    catch ME
        fprintf('%s: UPDATE FAILED: %s\n', mdl, ME.message);
        for c = 1:min(numel(ME.cause), 8)
            fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 160)));
        end
    end
    save_system(mdl);
    close_system(mdl, 0);
end
fprintf('\n=== recon_patch DONE ===\n');
end
