function rig_hinge_20260926()
% replace the rig's KT<->KB Parallel constraint with a Revolute Joint
% (the true knee hinge pin; frames come from the constraint's own ports)
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

par = [mdl '/Parallel'];
ph = get_param(par, 'PortHandles');
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
fprintf('Parallel far ports: %d\n', numel(farH));
if numel(farH) == 2
    nmA = get_param(get_param(farH(1), 'Parent'), 'Name');
    nmB = get_param(get_param(farH(2), 'Parent'), 'Name');
    delete_block(par);
    add_block('sm_lib/Joints/Revolute Joint', [mdl '/knee_hinge'], 'Position', [420 300 480 360]);
    nph = get_param([mdl '/knee_hinge'], 'PortHandles');
    np = [nph.RConn nph.LConn];
    add_line(mdl, farH(1), np(1));
    add_line(mdl, farH(2), np(2));
    fprintf('knee_hinge wired: %s <-> %s\n', nmA, nmB);
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
set_param(mdl, 'StopTime', '6');
out = sim(mdl);
LE = out.L_EXT.signals(1).values;
LF = out.L_FLX.signals(1).values;
FE = out.F_EXT.signals(1).values;
FF = out.F_FLX.signals(1).values;
fprintf('EXT length: %.4f -> %.4f..%.4f m | FLX: %.4f..%.4f m\n', LE(1), min(LE), max(LE), min(LF), max(LF));
fprintf('forces: EXT max %.1f N | FLX max %.1f N\n', max(abs(FE)), max(abs(FF)));
r = corrcoef(LE(:), LF(:));
fprintf('corr(L_EXT, L_FLX) = %.3f (want NEGATIVE for antagonism)\n', r(1, 2));
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== rig_hinge DONE ===\n');
end
