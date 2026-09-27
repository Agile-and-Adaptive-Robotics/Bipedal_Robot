function rig_knee_rebuild_v2_20260926()
% v2: record far ports FIRST, then:
%   weld Revolute(BL<->KB)      -> BL rides with KB
%   weld Cylindrical1(KT<->FL1) -> FL_001 rides with KT
%   weld Revolute2(KT<->FL2)    -> FL_002 rides with KT
%   delete bare: Cylindrical(KT<->BL), Revolute1(FL2<->KB)
%   knee_hinge Revolute between pinT_KT / pinT_KB (branch the welded lines)
%   re-point EXT/FLX anchors to the pin frames
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

joints = {'Revolute', 'Cylindrical1', 'Revolute2', 'Cylindrical', 'Revolute1'};
far = struct();
for w = 1:numel(joints)
    jb = [mdl '/' joints{w}];
    ph = get_param(jb, 'PortHandles');
    ports = [ph.RConn ph.LConn];
    fh = [];
    for pp = 1:numel(ports)
        ln = get_param(ports(pp), 'Line');
        if ln < 0, continue; end
        phs = [];
        try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
        for q = 1:numel(phs)
            h = phs(q);
            if ~isscalar(h) || h <= 0 || h == ports(pp), continue; end
            if isempty(fh) || ~any(fh == h)
                fh(end+1) = h; %#ok<SAGROW>
            end
        end
    end
    far.(joints{w}) = fh;
    fprintf('%s: %d far ports\n', joints{w}, numel(fh));
end

% weld three joints (far-port pairs), delete the other two bare
for w = {'Revolute', 'Cylindrical1', 'Revolute2'}
    fh = far.(w{1});
    delete_block([mdl '/' w{1}]);
    if numel(fh) == 2
        add_line(mdl, fh(1), fh(2));
    end
    fprintf('%s: deleted + welded\n', w{1});
end
for w = {'Cylindrical', 'Revolute1'}
    delete_block([mdl '/' w{1}]);
    fprintf('%s: deleted bare\n', w{1});
end

% knee hinge at FL_002 crank midpoint (z axis; frames = part axes)
ktPort = far.Revolute2(1);   % KT-side port of the old KT<->FL_002 hinge
kbPort = far.Revolute1(1);   % KB-side port of the old FL_002<->KB hinge
% identify which recorded port is on which body
if ~strcmp(get_param(get_param(ktPort, 'Parent'), 'Name'), 'x04_01_KT_R_003_1_RIGID')
    ktPort = far.Revolute2(2);
end
if ~strcmp(get_param(get_param(kbPort, 'Parent'), 'Name'), 'x04_02_KB_R_003_1_RIGID')
    kbPort = far.Revolute1(2);
end
fprintf('hinge anchors: %s / %s\n', ...
    get_param(get_param(ktPort, 'Parent'), 'Name'), ...
    get_param(get_param(kbPort, 'Parent'), 'Name'));

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
fprintf('knee_hinge added\n');

% re-point muscle anchors: cut dangling base lines, branch from pin frames
for pfx = {'EXT', 'FLX'}
    oph = get_param([mdl '/' pfx{1} '_origT'], 'PortHandles'); oP = [oph.RConn oph.LConn];
    iph = get_param([mdl '/' pfx{1} '_insT'], 'PortHandles'); iP = [iph.RConn iph.LConn];
    for pp = {oP(1), iP(1)}
        ln = get_param(pp{1}, 'Line');
        if ln > 0
            hasFar = false;
            phs = [];
            try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
            for q = 1:numel(phs)
                h = phs(q);
                if isscalar(h) && h > 0 && h ~= pp{1}, hasFar = true; end
            end
            if hasFar
                fprintf('%s: base line still connected, leaving it\n', pfx{1});
            else
                delete_line(ln);
            end
        end
    end
    if get_param(oP(1), 'Line') < 0, add_line(mdl, pK(2), oP(1)); end
    if get_param(iP(1), 'Line') < 0, add_line(mdl, pB(2), iP(1)); end
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
fprintf('=== rig_knee_rebuild_v2 DONE ===\n');
end
