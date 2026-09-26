function surgery1_20260926()
%% surgery1_20260926  Humanoid model fixes, pass 1:
%%   A. map ALL plain lines on key bodies (find hidden rigid welds)
%%   B. delete the two degenerate blocks (Cartesian, Parallel) in x09_BA_001_1
%%   C. delete the KT<->KB direct rigid weld lines (both knees)
%%   D. compile + save

sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
load_system(fullfile(sns, [mdl '.slx']));

H = containers.Map('KeyType', 'char', 'ValueType', 'any');  % body name -> handle
targets = {'x09_BA_001_1', ...
    'x09_BA_001_1/x02_01_PE_001_1_RIGID', 'x09_BA_001_1/x03_04_FH_L_001_1_RIGID', ...
    'x09_BA_001_1/x04_03_KT_L_001_1_RIGID', 'x09_BA_001_1/x04_04_KB_L_001_1_RIGID', ...
    'x09_BA_001_1/x04_01_KT_R_001_1_RIGID', 'x09_BA_001_1/x04_02_KB_R_001_1_RIGID'};
for t = 1:numel(targets)
    hb = getSimulinkBlockHandle([mdl '/' targets{t}]);
    if hb < 0
        [~, nm] = fileparts(targets{t});
        error('cannot resolve block %s', targets{t});
    end
    H(targets{t}) = hb;
    fprintf('%-54s handle=%d\n', targets{t}, hb);
end

%% ---- A. plain-line map -----------------------------------------------------
fprintf('\n--- PLAIN-LINE MAP (key bodies) ---\n');
for t = 1:numel(targets)
    body = targets{t};
    hb = H(body);
    ph = get_param(hb, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    for pp = 1:numel(ports)
        far = far_ends(ports(pp));
        if isempty(far), continue; end
        for f = 1:size(far, 1)
            fprintf('%-34s port%d  ->  %-38s (port %s)\n', body, pp, far{f, 1}, far{f, 2});
        end
    end
end

%% ---- B. delete degenerate blocks -------------------------------------------
for b = {'x09_BA_001_1/Cartesian', 'x09_BA_001_1/Parallel'}
    try
        delete_block([mdl '/' b{1}]);
        fprintf('deleted block %s\n', b{1});
    catch ME
        fprintf('could not delete %s: %s\n', b{1}, ME.message);
    end
end

%% ---- C. delete KT<->KB rigid weld lines ------------------------------------
for leg = {'L', 'R'}
    if strcmp(leg{1}, 'L')
        kt = 'x09_BA_001_1/x04_03_KT_L_001_1_RIGID'; kb = 'x09_BA_001_1/x04_04_KB_L_001_1_RIGID';
    else
        kt = 'x09_BA_001_1/x04_01_KT_R_001_1_RIGID'; kb = 'x09_BA_001_1/x04_02_KB_R_001_1_RIGID';
    end
    nCut = 0;
    ph = get_param(H(kt), 'PortHandles');
    ports = [ph.RConn ph.LConn];
    for pp = 1:numel(ports)
        far = far_ends(ports(pp));
        for f = 1:size(far, 1)
            if isnumeric(far{f, 3}) && far{f, 3} == H(kb)   % direct line to KB (no intermediate block)
                try
                    delete_line(mdl, far{f, 4}, far{f, 5});
                    nCut = nCut + 1;
                    fprintf('cut knee weld %s port%d <-> %s port %s\n', kt, pp, kb, far{f, 2});
                catch ME
                    fprintf('delete_line failed: %s\n', ME.message);
                end
            end
        end
    end
    fprintf('%s knee: cut %d rigid line(s)\n', leg{1}, nCut);
end

%% ---- D. compile + save -----------------------------------------------------
save_system(mdl);
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
fprintf('=== surgery1 DONE ===\n');
end

function far = far_ends(porth)
% far ends of the connection at porth: {blockName, portPath, blockHandle, myPortPath?...}
% returns rows: {farBlockName, farPortPath, farBlockHandle, thisPortPath, farPortPath}
far = {};
ln = -1;
try, ln = get_param(porth, 'Line'); catch, end
if ln < 0, return; end
thisPath = '';
try, thisPath = get_param(porth, 'PortPath'); catch, end
if isempty(thisPath)
    try, thisPath = port_ref_to_path(porth); catch, end
end
phs = [];
try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
br = -1;
try, br = get_param(porth, 'Branch'); catch, end
if isscalar(br) && br > 0
    try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
end
seen = [];
for q = 1:numel(phs)
    h = phs(q);
    if ~isscalar(h) || h <= 0 || h == porth, continue; end
    if ~isempty(seen) && any(seen == h), continue; end
    seen(end+1) = h; %#ok<AGROW>
    fp = '';
    try, fp = get_param(h, 'PortPath'); catch, end
    if isempty(fp)
        try, fp = port_ref_to_path(h); catch, end
    end
    try
        bh = get_param(get_param(h, 'Parent'), 'Handle');
        bn = get_param(bh, 'Name');
    catch
        continue;
    end
    far(end+1, :) = {bn, fp, bh, thisPath, fp}; %#ok<AGROW>
end
end

function p = port_ref_to_path(h)
% build 'blkPath/portName#k' for a port handle when PortPath is unavailable
par = get_param(h, 'Parent');
p = '';
try
    pt = get_param(h, 'PortType');
    ph = get_param(par, 'PortHandles');
    fn = fieldnames(ph);
    for k = 1:numel(fn)
        v = ph.(fn{k});
        idx = find(v == h, 1);
        if ~isempty(idx)
            p = sprintf('%s/%s#%d', par, fn{k}, idx);
            return;
        end
    end
catch
end
if isempty(p)
    error('cannot resolve port path for handle %d', h);
end
end
