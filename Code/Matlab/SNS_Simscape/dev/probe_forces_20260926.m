function probe_forces_20260926()
load_system('sm_lib');
l3 = find_system('sm_lib/Forces and Torques', 'SearchDepth', 2, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'Block');
for k = 1:numel(l3)
    fprintf('FT: %s\n', l3{k});
end
try
    dp = fieldnames(get_param('sm_lib/Forces and Torques/External Force', 'DialogParameters'));
    fprintf('ExtForce params: %s\n', strjoin(dp', ', '));
catch ME
    fprintf('ExtForce probe: %s\n', ME.message);
end
l4 = find_system('sm_lib/Frames and Transforms', 'SearchDepth', 1, 'Type', 'Block');
for k = 1:numel(l4), fprintf('FRAMES: %s\n', l4{k}); end
end
