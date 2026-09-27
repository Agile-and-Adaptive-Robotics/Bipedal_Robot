function repair_mux2_20260926()
% delete dangling mux-in lines; wire product OUTPUTS -> mux inputs by handle
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

for pfx = {'EXT', 'FLX'}
    % 1. cut whatever currently feeds the mux inputs (dangling lines)
    mh = get_param([mdl '/' pfx{1} '_Mux'], 'PortHandles');
    pc = get_param([mdl '/' pfx{1} '_Mux'], 'PortConnectivity');
    for j = 1:numel(mh.Inport)
        ln = get_param(mh.Inport(j), 'Line');
        if ln > 0
            hasSrc = false;
            for q = 1:numel(pc)
                try
                    h = pc(q).SrcBlockHandle;
                    if isscalar(h) && h > 0, hasSrc = true; end
                catch, end
            end
            if ~hasSrc
                delete_line(ln);
                fprintf('%s: cut dangling mux-in %d\n', pfx{1}, j);
            end
        end
    end
    % 2. wire product outputs -> mux inputs (handles)
    muxh = get_param([mdl '/' pfx{1} '_Mux'], 'PortHandles');
    for c = 1:3
        prh = get_param([mdl '/' sprintf('%s_F%d', pfx{1}, c)], 'PortHandles');
        if get_param(muxh.Inport(c), 'Line') < 0
            add_line(mdl, prh.Outport(1), muxh.Inport(c));
            fprintf('%s: F%d out -> mux in%d\n', pfx{1}, c, c);
        end
    end
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
set_param(mdl, 'StopTime', '4');
out = sim(mdl);
L_E = out.L_EXT; L_F = out.L_FLX; F_E = out.F_EXT; F_F = out.F_FLX; VE = out.V_E;
v = L_E.signals(1).values;
fprintf('EXT length: t0 %.4f  min %.4f  max %.4f m\n', v(1), min(v), max(v));
v = L_F.signals(1).values;
fprintf('FLX length: t0 %.4f  min %.4f  max %.4f m\n', v(1), min(v), max(v));
fprintf('EXT force max %.1f N | FLX force max %.1f N\n', max(abs(F_E.signals(1).values)), max(abs(F_F.signals(1).values)));
fprintf('V_E range %.1f..%.1f mV (Vrest -52)\n', min(VE.signals(1).values), max(VE.signals(1).values));
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== repair_mux2 DONE ===\n');
end
