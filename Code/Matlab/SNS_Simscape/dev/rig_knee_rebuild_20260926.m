function rig_knee_rebuild_20260926()
% Reduce the rig knee to a single z-axis revolute at the FL_002 crank centre:
%   delete Cylindrical, Cylindrical1, Revolute, Revolute1, Revolute2
%   weld BL->KB, FL_001->KT, FL_002->KT
%   knee_hinge Revolute KT<->KB at the crank midpoint (z axis)
%   re-point muscle anchor transforms to the new pin frames
%   (calibration script then re-zeros the muscle geometry)
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

joints = {'Cylindrical', 'Cylindrical1', 'Revolute', 'Revolute1', 'Revolute2'};
% far ports per joint, in order: {bodyAfar, bodyBfar} -> weld pairs:
%   Cylindrical  KT<->BL   : weld BL -> KT?? (BL bolted set is to KB) -> weld BL to KB side
%   Cylindrical1 KT<->FL001: weld FL001 -> KT
%   Revolute     KB<->BL   : (BL already welded to KB by first) -> welds redundant; just cut
%   Revolute1    FL002<->KB: weld FL002 -> KT is wanted, not KB; KB-side port freed
%   Revolute2    FL002<->KT: weld FL002 -> KT
% simple rule: weld Cylindrical's far pair, Cylindrical1's far pair, Revolute1
% KB-side->FL002 is dropped (FL002 welded to KT via Revolute2), Revolute (BL<->KB) welded.
weldPairs = {'Cylindrical', 'Cylindrical1', 'Revolute', 'Revolute2'};
for w = 1:numel(weldPairs)
    jb = [mdl '/' weldPairs{w}];
    ph = get_param(jb, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    farH = [];
    for pp = 1:numel(ports)
        ln = get_param(ports(pp), 'Line');
        if ln < 0, continue; end
        phs = [];
        try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
        for q = 1:numel(phs)
            h = phs(q);
            if ~isscalar(h) || h <= 0 || h == ports(pp), continue; end
            farH(end+1) = h; %#ok<SAGROW>
        end
    end
    if numel(farH) == 2
        delete_block(jb);
        add_line(mdl, farH(1), farH(2));
        fprintf('%s: deleted + welded far ports\n', weldPairs{w});
    else
        fprintf('%s: %d far ports, skipped\n', weldPairs{w}, numel(farH));
    end
end
% Revolute (BL<->KB): BL already welded to KB through Cylindrical weld; delete
% the joint and weld its ports too (harmless duplicate rigidity within one body)
try
    jb = [mdl '/Revolute'];
    ph = get_param(jb, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    farH = [];
    for pp = 1:numel(ports)
        ln = get_param(ports(pp), 'Line');
        if ln < 0, continue; end
        phs = [];
        try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
        for q = 1:numel(phs)
            h = phs(q);
            if ~isscalar(h) || h <= 0 || h == ports(pp), continue; end
            farH(end+1) = h; %#ok<SAGROW>
        end
    end
    if numel(farH) == 2
        delete_block(jb);
        add_line(mdl, farH(1), farH(2));
        fprintf('Revolute: deleted + welded\n');
    end
catch ME
    fprintf('Revolute: %s\n', ME.message);
end

% ---- knee hinge at the FL_002 crank midpoint ---------------------------------
kt = [mdl '/x04_01_KT_R_003_1_RIGID'];
kb = [mdl '/x04_02_KB_R_003_1_RIGID'];
% branch frames: use the muscle anchor transforms' source ports is gone; branch
% off the welded lines: find any free frame port pair on KT and KB via new
% Rigid Transforms connected to the EXISTING weld lines (branching)
ktPort = line_partner_port([mdl '/EXT_origT'], kt);   % EXT_origT B is wired to a KT-side port
kbPort = line_partner_port([mdl '/EXT_insT'], kb);

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
add_line(mdl, ktPort, pK(1));
add_line(mdl, kbPort, pB(1));
add_line(mdl, pK(2), np(1));
add_line(mdl, pB(2), np(2));
fprintf('knee_hinge added (z axis)\n');

% ---- re-point muscle anchors to the new pin frames ---------------------------
for pfx = {'EXT', 'FLX'}
    oB = line_partner_port([mdl '/' pfx{1} '_origT'], kt);   % not used; we cut by line
    % cut origT/insT base lines and re-branch from the pin transforms
    oph = get_param([mdl '/' pfx{1} '_origT'], 'PortHandles'); oP = [oph.RConn oph.LConn];
    iph = get_param([mdl '/' pfx{1} '_insT'], 'PortHandles'); iP = [iph.RConn iph.LConn];
    ln = get_param(oP(1), 'Line');
    if ln > 0, delete_line(ln); end
    ln = get_param(iP(1), 'Line');
    if ln > 0, delete_line(ln); end
    add_line(mdl, pK(2), oP(1));
    add_line(mdl, pB(2), iP(1));
    fprintf('%s anchors re-pointed to pin frames\n', pfx{1});
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
fprintf('=== rig_knee_rebuild DONE (now run calibrate_muscles_20260926) ===\n');
end

function port = line_partner_port(jointBlk, bodyBlk)
ph = get_param(jointBlk, 'PortHandles');
ports = [ph.RConn ph.LConn];
hb = get_param(bodyBlk, 'Handle');
port = [];
for pp = 1:numel(ports)
    ln = get_param(ports(pp), 'Line');
    if ln < 0, continue; end
    phs = [];
    try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
    for q = 1:numel(phs)
        h = phs(q);
        if ~isscalar(h) || h <= 0 || h == ports(pp), continue; end
        try
            if get_param(get_param(h, 'Parent'), 'Handle') == hb
                port = ports(pp);
                return;
            end
        catch
        end
    end
end
assert(~isempty(port), 'no line from %s to %s', jointBlk, bodyBlk);
end
