function probe_blockset()
% List blocks in mjLib and the plant block's dialog parameters.
load_system('mjLib');
blks = find_system('mjLib', 'SearchDepth', 2, 'LookUnderMasks', 'all', 'Type', 'Block');
for i = 1:numel(blks)
    fprintf('%s  [%s]\n', blks{i}, get_param(blks{i}, 'BlockType'));
end
% Dialog parameters of the first non-subsystem block (the plant)
for i = 1:numel(blks)
    bt = get_param(blks{i}, 'BlockType');
    if ~strcmp(bt, 'SubSystem') && ~strcmp(bt, 'BlockReference')
        dp = get_param(blks{i}, 'DialogParameters');
        fprintf('\nDialog params of %s:\n', blks{i});
        fn = fieldnames(dp);
        for k = 1:numel(fn)
            try
                v = get_param(blks{i}, fn{k});
                if ischar(v) || isstring(v) || isnumeric(v)
                    fprintf('  %-24s = %s\n', fn{k}, char(string(v)));
                else
                    fprintf('  %-24s (mask)\n', fn{k});
                end
            catch
                fprintf('  %-24s (no value)\n', fn{k});
            end
        end
        break
    end
end
bdclose('mjLib');
end
