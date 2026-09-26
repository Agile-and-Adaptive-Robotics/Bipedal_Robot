function fix_dims2_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));
set_param([sub '/spine_col'], 'CylinderRadius', '0.06', 'CylinderLength', '0.45', 'Density', '1200');
% verify all four solids
for b = {'foot_R','foot_L','spine_col','spine_lump'}
    try
        d = get_param([sub '/' b{1}], 'Density');
        fprintf('%s density=%s\n', b{1}, d);
    catch, end
    try
        fprintf('  dims=%s\n', get_param([sub '/' b{1}], 'BrickDimensions'));
    catch, end
    try
        fprintf('  cyl r=%s len=%s\n', get_param([sub '/' b{1}], 'CylinderRadius'), get_param([sub '/' b{1}], 'CylinderLength'));
    catch, end
end
try
    set_param(mdl, 'StopTime', '0.01');
    set_param(mdl, 'SimulationCommand', 'update');
    fprintf('UPDATE OK\n');
catch ME
    fprintf('UPDATE FAILED: %s\n', ME.message);
end
save_system(mdl);
close_system(mdl, 0);
end
