function probe_lib2_20260926()
load_system('sm_lib');
l4 = find_system('sm_lib', 'SearchDepth', 5, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
want = {'cyl', 'mass'};
for k = 1:numel(l4)
    nm = get_param(l4{k}, 'Name');
    if contains(lower(nm), want)
        fprintf('FOUND: %s (mask: %s)\n', l4{k}, get_param(l4{k}, 'MaskType'));
    end
end
end
