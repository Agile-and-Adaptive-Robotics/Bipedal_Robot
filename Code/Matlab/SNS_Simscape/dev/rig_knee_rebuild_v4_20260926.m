function rig_knee_rebuild_v4_20260926()
% v4 — explicit reduction, all numbers precomputed from the XML mates:
%   keep ports: KB-side of old Revolute (=Concentric16 frame on KB),
%               KT-side of old Cylindrical (=Concentric19 frame on KT),
%               FL_002-side of old Revolute2 (=Concentric13 frame on FL_002)
%   welds:      BL<->KB, FL_001<->KT, FL_002<->KT (at its crank pin)
%   hinge:      knee_hinge KT<->KB at crank midpoint (z axis, part-frame axes)
%   muscles:    re-anchored to the pin frames, then recalibrated separately
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

bodyOf = @(h) get_param(get_param(h, 'Parent'), 'Name');
joints = {'Revolute', 'Cylindrical1', 'Revolute2', 'Cylindrical', 'Revolute1'};
far = struct();
for w = 1:numel(joints)
    jb = [mdl '/' joints{w}];
    ph = get_param(jb, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    fh = zeros(0, 1);
    for pp = 1:numel(ports)
        fh = [fh; far_ports_of(ports(pp))]; %#ok<AGROW>
    end
    far.(joints{w}) = fh;
end
h_KB   = pick_single(far.Revolute,    'x04_02_KB_R_003_1_RIGID');
h_BL   = pick_single(far.Revolute,    'x04_05_BL_001_1_RIGID');
h_KT   = pick_single(far.Cylindrical, 'x04_01_KT_R_003_1_RIGID');
h_FL1  = pick_single(far.Cylindrical1, 'x04_06_FL_001_1_RIGID');
h_FL1k = pick_single(far.Cylindrical1, 'x04_01_KT_R_003_1_RIGID');
h_FL2  = pick_single(far.Revolute2,   'x04_06_FL_002_1_RIGID');

% ---- delete all five joints ---------------------------------------------------
for w = 1:numel(joints)
    delete_block([mdl '/' joints{w}]);
end
fprintf('5 joints deleted\n');

% BL -> KB (weld at their own hinge frames)
add_line(mdl, h_BL, h_KB);
% FL_001 -> KT
add_line(mdl, h_FL1, h_FL1k);

% FL_002 -> KT: frame on KT at the crank pin, welded to the FL_002 port
deal_frame('fl2T', h_KT, [-0.0178 0.0152 -0.0740], mdl);
phF = get_param([mdl '/fl2T'], 'PortHandles'); pF = [phF.RConn phF.LConn];
add_line(mdl, pF(2), h_FL2);

% knee hinge frames at the crank midpoint
deal_frame('pinT_KT', h_KT, [-0.0061 -0.0077 -0.0745], mdl);
deal_frame('pinT_KB', h_KB, [0.0281 0.0242 0.0005], mdl);
phK = get_param([mdl '/pinT_KT'], 'PortHandles'); pK = [phK.RConn phK.LConn];
phB = get_param([mdl '/pinT_KB'], 'PortHandles'); pB = [phB.RConn phB.LConn];
add_block('sm_lib/Joints/Revolute Joint', [mdl '/knee_hinge'], 'Position', [200 560 260 620]);
nh = get_param([mdl '/knee_hinge'], 'PortHandles'); np = [nh.RConn nh.LConn];
add_line(mdl, pK(2), np(1));
add_line(mdl, pB(2), np(2));
fprintf('knee_hinge added\n');

% ---- re-point muscle anchors ---------------------------------------------------
for pfx = {'EXT', 'FLX'}
    oph = get_param([mdl '/' pfx{1} '_origT'], 'PortHandles'); oP = [oph.RConn oph.LConn];
    iph = get_param([mdl '/' pfx{1} '_insT'], 'PortHandles'); iP = [iph.RConn iph.LConn];
    for pr = {oP(1), iP(1)}
        ln = get_param(pr{1}, 'Line');
        if ln > 0, delete_line(ln); end
    end
    add_line(mdl, pK(2), oP(1));
    add_line(mdl, pB(2), iP(1));
    fprintf('%s anchors re-pointed\n', pfx{1});
end

try
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
catch ME
    fprintf('UPDATE FAILED: %s\n', ME.message);
    for c = 1:min(numel(ME.cause), 6)
        fprintf('  CAUSE %d: %s\n', c, ME.cause{c}.message(1:min(end, 170)));
    end
    save_system(mdl);
    close_system(mdl, 0);
    return;
end
save_system(mdl);
close_system(mdl, 0);
fprintf('=== rig_knee_rebuild_v4 DONE ===\n');
end

function deal_frame(nm, srcPort, dT, mdl)
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/' nm], 'Position', [40 700 100 750]);
set_param([mdl '/' nm], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', mat2str(dT, 6), 'RotationMethod', 'None');
ph = get_param([mdl '/' nm], 'PortHandles');
add_line(mdl, srcPort, ph.RConn(1));
end

function fp = far_ports_of(porth)
fp = zeros(0, 1);
ln = get_param(porth, 'Line');
if ln < 0, return; end
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
    if any(seen == h), continue; end
    seen(end+1) = h; %#ok<AGROW>
    try
        pn = get_param(get_param(h, 'Parent'), 'Name');
    catch
        continue;
    end
    if ~isempty(regexpi(pn, '_RIGID$', 'once'))
        fp(end+1, 1) = h; %#ok<AGROW>
    end
end
end

function h = pick_single(fh, bodyA)
h = [];
for k = 1:size(fh, 1)
    pn = get_param(get_param(fh(k, 1), 'Parent'), 'Name');
    if strcmp(pn, bodyA), h = fh(k, 1); end
end
assert(~isempty(h), 'pick_single failed for %s', bodyA);
end
