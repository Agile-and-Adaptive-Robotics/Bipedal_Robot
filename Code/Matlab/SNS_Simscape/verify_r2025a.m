%VERIFY_R2025A - compare original vs exported-copy block inventories (run from this folder)
pairs = {'SNS_Library', 'SNS_Library_R2025a'; 'KneeReflexDemo', 'KneeReflexDemo_R2025a'};
for k = 1:size(pairs, 1)
    src = pairs{k, 1}; dst = pairs{k, 2};
    load_system(src);
    load_system(dst);
    ha = find_system(src, 'LookUnderMasks', 'all', 'Type', 'Block');
    hb = find_system(dst, 'LookUnderMasks', 'all', 'Type', 'Block');
    ta = sort(cellfun(@(h) get_param(h, 'BlockType'), ha, 'UniformOutput', false));
    tb = sort(cellfun(@(h) get_param(h, 'BlockType'), hb, 'UniformOutput', false));
    fprintf('%s: orig %d blocks, exported %d blocks\n', src, numel(ha), numel(hb));
    if isequal(ta, tb)
        fprintf('  BLOCKTYPES_IDENTICAL\n');
    else
        fprintf('  BLOCKTYPES_DIFFER:\n');
        u = unique([ta; tb]);
        for i = 1:numel(u)
            na = sum(strcmp(ta, u{i})); nb = sum(strcmp(tb, u{i}));
            if na ~= nb
                fprintf('    %-22s orig %d  exported %d\n', u{i}, na, nb);
            end
        end
    end
    if strcmp(src, 'KneeReflexDemo')
        links = find_system(src, 'LookUnderMasks', 'none', 'Type', 'Block');
        nLink = 0; libNames = {};
        for h = 1:numel(links)
            if ~strcmp(get_param(links{h}, 'LinkStatus'), 'none')
                nLink = nLink + 1;
                libNames{end+1} = strtok(get_param(links{h}, 'ReferenceBlock'), '/');
            end
        end
        fprintf('  Library-linked top-level blocks: %d', nLink);
        if nLink > 0, fprintf(' (from: %s)', strjoin(unique(libNames), ', ')); end
        fprintf('\n');
    end
    close_system(src, 0);
    close_system(dst, 0);
end
fprintf('VERIFY_DONE\n');
