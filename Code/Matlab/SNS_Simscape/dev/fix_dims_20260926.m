function fix_dims_20260926()
sns = fileparts(fileparts(mfilename('fullpath')));
mdl = 'mdl_humanoid_lower_ah001_imported';
sub = [mdl '/x09_BA_001_1'];
load_system(fullfile(sns, [mdl '.slx']));

% discover cylinder param names
dp = fieldnames(get_param([sub '/spine_col'], 'DialogParameters'));
cyl = dp(contains(dp, 'Cyl', 'IgnoreCase', true));
fprintf('cyl params: %s\n', strjoin(cyl', ', '));

for b = {'foot_R', 'foot_L'}
    set_param([sub '/' b{1}], 'BrickDimensions', '[0.24 0.035 0.09]', 'Density', '1200');
end
set_param([sub '/spine_lump'], 'BrickDimensions', '[0.1 0.1 0.1]', 'Density', '61090.6438814214');
try
    set_param([sub '/spine_col'], 'CylindricalRadius', '0.06', 'CylindricalLength', '0.45', 'Density', '1200');
catch ME
    fprintf('cyl set failed: %s\n', ME.message);
    fn = fieldnames(get_param([sub '/spine_col'], 'DialogParameters'));
    for k = 1:numel(fn)
        if contains(fn{k}, 'Rad') || contains(fn{k}, 'Len')
            fprintf('  candidate: %s\n', fn{k});
        end
    end
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
