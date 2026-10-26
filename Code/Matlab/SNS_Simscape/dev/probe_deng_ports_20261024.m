function probe_deng_ports_20261024()
% Read-only: enumerate port structure of the Deng trio to pin the
% port-width mismatch root cause (2026-10-24).
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root); addpath(fullfile(root, 'demos'));

for nm = {'SNS_Deng_RG', 'SNS_Deng_PF', 'SNS_Deng_CPGDemo'}
    m = nm{1};
    try
        load_system(m);
        fprintf('--- %s ---\n', m);
        ports = find_system(m, 'SearchDepth', 1, 'BlockType', 'Outport');
        for k = 1:numel(ports)
            fprintf('  outport: %s | Port=%s\n', ports{k}, ...
                    get_param(ports{k}, 'Port'));
        end
        ins = find_system(m, 'SearchDepth', 1, 'BlockType', 'Inport');
        for k = 1:numel(ins)
            fprintf('  inport:  %s | Port=%s\n', ins{k}, ...
                    get_param(ins{k}, 'Port'));
        end
        if strcmp(m, 'SNS_Deng_CPGDemo')
            % list Model blocks + their connections
            mb = find_system(m, 'SearchDepth', 1, 'LookUnderMasks', 'all', ...
                             'ReferenceBlock', '');
            blks = find_system(m, 'SearchDepth', 1, 'Type', 'Block');
            for k = 1:numel(blks)
                rt = get_param(blks{k}, 'ReferenceBlock');
                if ~isempty(rt)
                    fprintf('  model-block: %s -> %s\n', blks{k}, rt);
                end
            end
            l = find_system(m, 'SearchDepth', 1, 'FindAll', 'on', ...
                            'Type', 'line');
            fprintf('  %d top-level lines\n', numel(l));
        end
        close_system(m, 0);
    catch e
        fprintf('--- %s FAILED: %s ---\n', m, e.message);
    end
end
fprintf('probe done\n');
end
