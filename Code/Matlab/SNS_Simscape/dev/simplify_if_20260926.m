function simplify_if_20260926()
% the IF f-input is scalar (force along the line between frames):
% delete the vector machinery, wire F magnitude straight in
here = fileparts(mfilename('fullpath'));
sns = fileparts(here);
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));

for pfx = {'EXT', 'FLX'}
    % remove sps input from mux first
    sph = get_param([mdl '/' pfx{1} '_sps'], 'PortHandles');
    ln = get_param(sph.Inport(1), 'Line');
    if ln > 0, delete_line(ln); end
    % sps input <- BPA force magnitude (scalar)
    add_line(mdl, [pfx{1} '_BPA/1'], [pfx{1} '_sps/1'], 'autorouting', 'on');
    % delete now-unused blocks: Mux, F1..F3, invD, one, converters c1..c3
    junk = {[pfx{1} '_Mux'], [pfx{1} '_F1'], [pfx{1} '_F2'], [pfx{1} '_F3'], ...
        [pfx{1} '_invD'], [pfx{1} '_one'], ...
        [pfx{1} '_c1'], [pfx{1} '_c2'], [pfx{1} '_c3']};
    for j = 1:numel(junk)
        try
            delete_block([mdl '/' junk{j}]);
        catch
        end
    end
    % also turn off SenseX/Y/Z on the sensor (keep only Dist)
    set_param([mdl '/' pfx{1} '_TS'], 'SenseX', 'off', 'SenseY', 'off', 'SenseZ', 'off', 'SenseDist', 'on');
    fprintf('%s simplified\n', pfx{1});
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
fprintf('=== simplify_if DONE ===\n');
end
