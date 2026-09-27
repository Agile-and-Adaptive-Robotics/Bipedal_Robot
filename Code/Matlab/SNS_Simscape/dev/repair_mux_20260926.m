function repair_mux_20260926()
% repair: force Mux widths, verify wiring, compile, sim, report
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

for pfx = {'EXT', 'FLX'}
    set_param([mdl '/' pfx{1} '_Mux'], 'Inputs', '3');
    % verify the three product->mux lines exist
    mh = get_param([mdl '/' pfx{1} '_Mux'], 'PortHandles');
    fprintf('%s_Mux inputs connected: ', pfx{1});
    for j = 1:numel(mh.Inport)
        ln = get_param(mh.Inport(j), 'Line');
        fprintf('%d ', ln > 0);
    end
    fprintf('\n');
    % verify product inputs: F_c ports 1(F) 2(x) 3(y/z)
    for c = 1:3
        ph = get_param([mdl '/' sprintf('%s_F%d', pfx{1}, c)], 'PortHandles');
        n = numel(ph.Inport);
        ok = zeros(1, n);
        for j = 1:n, ok(j) = get_param(ph.Inport(j), 'Line') > 0; end
        fprintf('  %s_F%d inputs(%d): %s\n', pfx{1}, c, n, mat2str(ok));
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
t = L_E.time; v = L_E.signals(1).values;
fprintf('EXT length: t0 %.4f  min %.4f  max %.4f m\n', v(1), min(v), max(v));
v = L_F.signals(1).values;
fprintf('FLX length: t0 %.4f  min %.4f  max %.4f m\n', v(1), min(v), max(v));
fprintf('EXT force max %.1f N | FLX force max %.1f N\n', max(abs(F_E.signals(1).values)), max(abs(F_F.signals(1).values)));
fprintf('V_E range %.1f..%.1f mV (Vrest -52)\n', min(VE.signals(1).values), max(VE.signals(1).values));
save(fullfile(here, '..', 'results', 'rig_bpa_sns_20260926.mat'), 'out');
close_system(mdl, 0);
fprintf('=== repair_mux DONE ===\n');
end
