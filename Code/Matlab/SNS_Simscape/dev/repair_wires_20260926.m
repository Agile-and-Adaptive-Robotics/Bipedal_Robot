function repair_wires_20260926()
% wire the missing F and 1/dist inputs into the force-component products
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

for pfx = {'EXT', 'FLX'}
    for c = 1:3
        add_line(mdl, [pfx{1} '_BPA/1'], sprintf('%s_F%d/1', pfx{1}, c), 'autorouting', 'on');
        add_line(mdl, [pfx{1} '_invD/1'], sprintf('%s_F%d/3', pfx{1}, c), 'autorouting', 'on');
    end
    % sanity: log invD? skip. verify all inputs connected now
    for c = 1:3
        ph = get_param([mdl '/' sprintf('%s_F%d', pfx{1}, c)], 'PortHandles');
        ok = zeros(1, numel(ph.Inport));
        for j = 1:numel(ph.Inport), ok(j) = get_param(ph.Inport(j), 'Line') > 0; end
        fprintf('%s_F%d inputs: %s\n', pfx{1}, c, mat2str(ok));
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
fprintf('=== repair_wires DONE ===\n');
end
