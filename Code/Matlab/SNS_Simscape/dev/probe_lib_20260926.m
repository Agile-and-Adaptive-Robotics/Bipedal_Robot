function probe_lib_20260926()
load_system('sm_lib');
top = find_system('sm_lib', 'SearchDepth', 1, 'Type', 'Block');
fprintf('sm_lib top level: %s\n', strjoin(top, ' | '));
l2 = find_system('sm_lib', 'SearchDepth', 2, 'LookUnderMasks', 'all', 'Type', 'Block');
want = {'Brick', 'Cylinder', 'Point Mass', 'Sphere', 'File Solid', 'Extruded Solid', 'Revolved Solid'};
for k = 1:numel(l2)
    nm = get_param(l2{k}, 'Name');
    if any(strcmpi(nm, want))
        fprintf('FOUND: %s (mask: %s)\n', l2{k}, get_param(l2{k}, 'MaskType'));
    end
end
end
