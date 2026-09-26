function probe_blocks_20260926()
load_system('sm_lib');
for lib = {'sm_lib/Frames and Transforms/Rigid Transform', 'sm_lib/Body/Elements/Brick', ...
           'sm_lib/Body/Elements/Cylinder', 'sm_lib/Body/Elements/Point Mass'}
    try
        load_system(lib{1});
        dp = fieldnames(get_param(lib{1}, 'DialogParameters'));
        fprintf('%s:\n', lib{1});
        for k = 1:numel(dp)
            v = '';
            try, v = get_param(lib{1}, dp{k}); catch, end
            if ischar(v) && ~isempty(v), v = [' = ' v]; else, v = ''; end
            fprintf('   %s%s\n', dp{k}, v);
        end
    catch ME
        fprintf('%s: LOAD FAILED %s\n', lib{1}, ME.message);
    end
end
end
