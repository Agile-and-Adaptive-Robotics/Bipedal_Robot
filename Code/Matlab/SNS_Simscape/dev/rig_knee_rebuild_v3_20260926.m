function rig_knee_rebuild_v3_20260926()
% v3: record (port,parent) pairs including branch walks; weld brackets to the
% correct side; knee hinge at the crank midpoint; re-point muscle anchors.
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
    fh = [];
    for pp = 1:numel(ports)
        fh = [fh; far_ports_of(ports(pp))]; %#ok<AGROW>
    end
    far.(joints{w}) = fh;
    for k = 1:size(fh, 1)
        fprintf('%s far: %s\n', joints{w}, bodyOf(fh(k, 1)));
    end
end

% picks
h_BL_KB   = pick_pair(far.Revolute,    'x04_05_BL_001_1_RIGID', 'x04_02_KB_R_003_1_RIGID');
h_KT_FL1  = pick_pair(far.Cylindrical1, 'x04_01_KT_R_003_1_RIGID', 'x04_06_FL_001_1_RIGID');
h_KT_FL2  = pick_pair(far.Revolute2,   'x04_01_KT_R_003_1_RIGID', 'x04_06_FL_002_1_RIGID');
kbPin     = pick_single(far.Revolute1, 'x04_02_KB_R_003_1_RIGID');

% welds
delete_block([mdl '/Revolute']);    add_line(mdl, h_BL_KB(1),  h_BL_KB(2));
delete_block([mdl '/Cylindrical1']); add_line(mdl, h_KT_FL1(1), h_KT_FL1(2));
delete_block([mdl '/Revolute2']);   add_line(mdl, h_KT_FL2(1), h_KT_FL2(2));
delete_block([mdl '/Cylindrical']);
delete_block([mdl '/Revolute1']);
fprintf('joints rebuilt\n');

% knee hinge
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/pinT_KT'], 'Position', [40 520 100 570]);
set_param([mdl '/pinT_KT'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[0.0036 -0.0083 -0.0375]', 'RotationMethod', 'None');
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/pinT_KB'], 'Position', [40 600 100 650]);
set_param([mdl '/pinT_KB'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[0.0122 0.0331 -0.0374]', 'RotationMethod', 'None');
add_block('sm_lib/Joints/Revolute Joint', [mdl '/knee_hinge'], 'Position', [200 560 260 620]);
phK = get_param([mdl '/pinT_KT'], 'PortHandles'); pK = [phK.RConn phK.LConn];
phB = get_param([mdl '/pinT_KB'], 'PortHandles'); pB = [phB.RConn phB.LConn];
nh = get_param([mdl '/knee_hinge'], 'PortHandles'); np = [nh.RConn nh.LConn];
add_line(mdl, h_KT_FL2(1), pK(1));
add_line(mdl, kbPin, pB(1));
add_line(mdl, pK(2), np(1));
add_line(mdl, pB(2), np(2));
fprintf('knee_hinge added\n');

% re-point muscle anchors (their old branches died with the deleted joints)
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
fprintf('=== rig_knee_rebuild_v3 DONE ===\n');
end

function fp = far_ports_of(porth)
% all far PORT HANDLES on bodies (RIGID blocks) attached via lines+branches
fp = zeros(0, 1);
ln = get_param(porth, 'Line');
if ln < 0, return; end
phs = [];
try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
% branches live at BOTH the line and the port level
br = -1;
try, br = get_param(porth, 'Branch'); catch, end
if isscalar(br) && br > 0
    try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
end
br = -1;
try, br = get_param(ln, 'BranchHandles'); catch, end
if isscalar(br) && br > 0
    phs = [phs br]; %#ok<AGROW>
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
    fprintf('    far: %s\n', pn);
    if ~isempty(regexpi(pn, '_RIGID$', 'once'))
        fp(end+1, 1) = h; %#ok<AGROW>
    end
end
end

function pair = pick_pair(fh, bodyA, bodyB)
ia = []; ib = [];
for k = 1:size(fh, 1)
    pn = get_param(get_param(fh(k, 1), 'Parent'), 'Name');
    if strcmp(pn, bodyA), ia = fh(k, 1); end
    if strcmp(pn, bodyB), ib = fh(k, 1); end
end
assert(~isempty(ia) && ~isempty(ib), 'pick_pair failed for %s/%s', bodyA, bodyB);
pair = [ia ib];
end

function h = pick_single(fh, bodyA)
h = [];
for k = 1:size(fh, 1)
    pn = get_param(get_param(fh(k, 1), 'Parent'), 'Name');
    if strcmp(pn, bodyA), h = fh(k, 1); end
end
assert(~isempty(h), 'pick_single failed for %s', bodyA);
end
