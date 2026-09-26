function probe_taps_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_leg_rig_ba003_imported';
load_system(fullfile(sns, [mdl '.slx']));
blks = find_system(mdl, 'SearchDepth', 1, 'Regexp', 'on', 'Name', '_t[123]w?$', ...
    'LookUnderMasks', 'all', 'Type', 'Block');
fprintf('leftover taps: %d\n', numel(blks));
for k = 1:numel(blks), fprintf('  %s\n', blks{k}); end
% sensor port state
for pfx = {'EXT', 'FLX'}
    ph = get_param([mdl '/' pfx{1} '_TS'], 'PortHandles');
    allp = [ph.LConn ph.RConn];
    fprintf('%s_TS ports=%d lines=', numel(allp));
    for j = 1:numel(allp)
        ln = -1;
        try, ln = get_param(allp(j), 'Line'); catch, end
        fprintf('%d ', ln > 0);
    end
    fprintf('\n');
    fprintf('  SenseX=%s SenseY=%s SenseZ=%s SenseDist=%s\n', ...
        get_param([mdl '/' pfx{1} '_TS'], 'SenseX'), ...
        get_param([mdl '/' pfx{1} '_TS'], 'SenseY'), ...
        get_param([mdl '/' pfx{1} '_TS'], 'SenseZ'), ...
        get_param([mdl '/' pfx{1} '_TS'], 'SenseDist'));
end
close_system(mdl, 0);
end
