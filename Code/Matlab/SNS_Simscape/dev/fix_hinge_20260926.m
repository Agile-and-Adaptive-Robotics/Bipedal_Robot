function fix_hinge_20260926()
% re-anchor the knee hinge on IDENTITY-rotation source frames:
%   KT side: the old Revolute2 KT-side port (smiData entry 25, R=I, at
%            (-0.00806, 0.01461, -0.0370) in KT coords)
%   KB side: the KB<->TI_007 bolt-weld frame (entry 9, R=I, at
%            (-0.03992, 0, 0.00349) in KB coords)
% new pin frames carry pure translations; hinge z = part z. Muscles re-pointed
% and recalibrated afterwards.
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

kt = [mdl '/x04_01_KT_R_003_1_RIGID'];
kb = [mdl '/x04_02_KB_R_003_1_RIGID'];
ti7 = [mdl '/x05_01_TI_R_007_1_RIGID'];

% KT-side identity frame = a FREE port on KT (the old Revolute2 leftover).
% Identify: KT ports that are free AND not the fl2T/pinT branch sources.
phK = get_param(kt, 'PortHandles');
kPorts = [phK.RConn phK.LConn];
ktSrc = [];
for j = 1:numel(kPorts)
    if get_param(kPorts(j), 'Line') < 0
        ktSrc = kPorts(j);
        break;
    end
end
assert(~isempty(ktSrc), 'no free port on KT');
fprintf('KT free port found (idx %d)\n', j);

% KB-side identity frame = the KB port welded to TI_007
phB = get_param(kb, 'PortHandles');
bPorts = [phB.RConn phB.LConn];
kbSrc = [];
for j = 1:numel(bPorts)
    ln = get_param(bPorts(j), 'Line');
    if ln < 0, continue; end
    phs = [];
    try, phs = [get_param(ln, 'SrcPortHandle') get_param(ln, 'DstPortHandle')]; catch, end
    br = -1;
    try, br = get_param(bPorts(j), 'Branch'); catch, end
    if isscalar(br) && br > 0
        try, phs = [phs get_param(br, 'BranchHandles')]; catch, end
    end
    for q = 1:numel(phs)
        h = phs(q);
        if isscalar(h) && h > 0 && h ~= bPorts(j)
            if strcmp(get_param(get_param(h, 'Parent'), 'Name'), 'x05_01_TI_R_007_1_RIGID')
                kbSrc = bPorts(j);
            end
        end
    end
end
assert(~isempty(kbSrc), 'no KB<->TI_007 weld port found');
fprintf('KB<->TI_007 weld port found\n');

% new pin frames (pure translations, identity orientation)
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/pinT_KT2'], 'Position', [40 760 100 810]);
set_param([mdl '/pinT_KT2'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[0.01166 -0.02291 -0.0005]', 'RotationMethod', 'None');
add_block('sm_lib/Frames and Transforms/Rigid Transform', [mdl '/pinT_KB2'], 'Position', [40 840 100 890]);
set_param([mdl '/pinT_KB2'], 'TranslationMethod', 'Cartesian', ...
    'TranslationCartesianOffset', '[0.0521 0.0331 -0.0409]', 'RotationMethod', 'None');
phK2 = get_param([mdl '/pinT_KT2'], 'PortHandles'); pK2 = [phK2.RConn phK2.LConn];
phB2 = get_param([mdl '/pinT_KB2'], 'PortHandles'); pB2 = [phB2.RConn phB2.LConn];
add_line(mdl, ktSrc, pK2(1));
add_line(mdl, kbSrc, pB2(1));

% rewire the hinge to the new frames
nh = get_param([mdl '/knee_hinge'], 'PortHandles');
np = [nh.RConn nh.LConn];
for pp = {np(1), np(2)}
    ln = get_param(pp{1}, 'Line');
    if ln > 0, delete_line(ln); end
end
add_line(mdl, pK2(2), np(1));
add_line(mdl, pB2(2), np(2));

% re-point muscle anchors to the new pin frames
for pfx = {'EXT', 'FLX'}
    oph = get_param([mdl '/' pfx{1} '_origT'], 'PortHandles'); oP = [oph.RConn oph.LConn];
    iph = get_param([mdl '/' pfx{1} '_insT'], 'PortHandles'); iP = [iph.RConn iph.LConn];
    for pr = {oP(1), iP(1)}
        ln = get_param(pr{1}, 'Line');
        if ln > 0, delete_line(ln); end
    end
    add_line(mdl, pK2(2), oP(1));
    add_line(mdl, pB2(2), iP(1));
end
% retire the old rotated pin frames
for nm = {'pinT_KT', 'pinT_KB'}
    try, delete_block([mdl '/' nm{1}]); catch, end
end
fprintf('hinge re-anchored on identity frames\n');

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
fprintf('=== fix_hinge DONE ===\n');
end
