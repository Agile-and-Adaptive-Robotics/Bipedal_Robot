function surgery2_20260926()
%% surgery2_20260926  Humanoid model fixes, pass 2 — free the biped:
%%   1. delete degenerate Cartesian1 (KT_R <-> Transform3)
%%   2. hip chains: delete stud/knob Revolutes 6,7,8,9, rejoin their far ports
%%      with rigid lines (stud+knob are BOLTED to the femur head in hardware;
%%      the only real hip DOF is the ball socket = Spherical PE<->knob)
%%   3. cut every plain line from a part body to an F#/Transform# frame block
%%      (the SW "grounded" assembly ties that weld the biped to the root)
%%   4. add pelvis_free 6-DOF joint World -> PE
%%   5. compile + save

sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

%% ---- 1. Cartesian1 ---------------------------------------------------------
try
    delete_block([sub '/Cartesian1']);
    fprintf('deleted Cartesian1\n');
catch ME
    fprintf('Cartesian1: %s\n', ME.message);
end

%% ---- 2. hip chains -> rigid ------------------------------------------------
hipRevs = {'Revolute6', 'Revolute7', 'Revolute8', 'Revolute9'};
for hr = 1:numel(hipRevs)
    blk = [sub '/' hipRevs{hr}];
    ph = get_param(blk, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    farH = [];
    for pp = 1:numel(ports)
        far = far_ends(ports(pp));
        for f = 1:size(far, 1)
            farH(end+1) = far{f, 6}; %#ok<SAGROW>  % far PORT HANDLE
        end
    end
    if numel(farH) ~= 2
        fprintf('%s: expected 2 far ports, got %d - SKIPPED\n', hipRevs{hr}, numel(farH));
        continue;
    end
    nmA = get_param(get_param(farH(1), 'Parent'), 'Name');
    nmB = get_param(get_param(farH(2), 'Parent'), 'Name');
    delete_block(blk);
    add_line(sub, farH(1), farH(2));
    fprintf('%s: welded  %s  <->  %s\n', hipRevs{hr}, nmA, nmB);
end

%% ---- 3. cut body-to-frame world ties ---------------------------------------
bodies = find_system(sub, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
nCut = 0;
for k = 1:numel(bodies)
    bn = get_param(bodies{k}, 'Name');
    if ~endsWith(bn, '_RIGID', 'IgnoreCase', true), continue; end   % part bodies only
    ph = get_param(bodies{k}, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    for pp = 1:numel(ports)
        far = far_ends(ports(pp));
        for f = 1:size(far, 1)
            fb = far{f, 1};
            if isempty(regexpi(fb, '^(F\d+|Transform\d*)$', 'once')), continue; end
            % plain line to a frame/boundary block -> world tie, cut it
            try
                ln = get_param(ports(pp), 'Line');
                delete_line(ln);
                nCut = nCut + 1;
                fprintf('cut tie: %s  ->  %s\n', bn, fb);
            catch ME
                fprintf('cut FAILED %s -> %s: %s\n', bn, fb, ME.message);
            end
            break;   % one cut per line; move to next port
        end
    end
end
fprintf('cut %d world tie line(s)\n', nCut);

%% ---- 4. pelvis_free 6-DOF --------------------------------------------------
% find a free frame port on the pelvis body
pe = [sub '/x02_01_PE_001_1_RIGID'];
ph = get_param(pe, 'PortHandles');
ports = [ph.RConn ph.LConn];
freePort = [];
for pp = 1:numel(ports)
    if get_param(ports(pp), 'Line') < 0
        freePort = ports(pp);
        break;
    end
end
assert(~isempty(freePort), 'no free frame port on pelvis body');

% world frame: the smimport-generated block at root level
allRoot = find_system(mdl, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
wf = {};
for k = 1:numel(allRoot)
    nm = get_param(allRoot{k}, 'Name');
    if contains(nm, 'World', 'IgnoreCase', true)
        wf{end+1} = allRoot{k}; %#ok<SAGROW>
        fprintf('World candidate: %s (ref: %s)\n', allRoot{k}, get_param(allRoot{k}, 'ReferenceBlock'));
    end
end
assert(~isempty(wf), 'no World Frame block found at root');
wfph = get_param(wf{1}, 'PortHandles');
wfPort = [wfph.RConn wfph.LConn];

add_block('sm_lib/Joints/6-DOF Joint', [sub '/pelvis_free'], 'Position', [200 380 260 440]);
wph = get_param([sub '/pelvis_free'], 'PortHandles');
wp = [wph.RConn wph.LConn];   % 1 = base (B), 2 = follower (F)

% World block inside the subassembly (all World blocks = the same inertial
% frame; avoids a cross-hierarchy connection through a boundary port that
% also carries the loose hardware welds)
add_block('sm_lib/Frames and Transforms/World Frame', [sub '/World_pelvis'], 'Position', [80 400 110 430]);
wiph = get_param([sub '/World_pelvis'], 'PortHandles');
wiPort = [wiph.RConn wiph.LConn];

add_line(sub, wiPort(1), wp(1));
add_line(sub, freePort, wp(2));
fprintf('pelvis_free 6-DOF wired (World_pelvis + pelvis)\n');

%% ---- 5. compile + save ------------------------------------------------------
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
fprintf('=== surgery2 DONE ===\n');
end

function far = far_ends(porth)
% far ends of the connection at porth:
% {farBlockName, farPortPath, farBlockHandle, thisPortPath, farPortPath}
far = {};
ln = -1;
try, ln = get_param(porth, 'Line'); catch, end
if ln < 0, return; end
thisPath = port_ref_to_path(porth);
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
    fp = port_ref_to_path(h);
    try
        bh = get_param(get_param(h, 'Parent'), 'Handle');
        bn = get_param(bh, 'Name');
    catch
        continue;
    end
    far(end+1, :) = {bn, fp, bh, thisPath, fp, h}; %#ok<AGROW>
end
end

function p = port_ref_to_path(h)
par = get_param(h, 'Parent');
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
error('cannot resolve port path for handle %d', h);
end
