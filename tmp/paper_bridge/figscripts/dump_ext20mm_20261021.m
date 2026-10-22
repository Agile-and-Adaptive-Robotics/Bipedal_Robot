% Dump the new 20 mm extensor test entries + per-run DAQ stats
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
base = 'D:/GitHub/Bipedal_Robot/Testing_Data/2026_06_Festo';

% xlsx: sheet inventory + dump non-empty sheets beyond the old ones
f = fullfile(base, 'Results_table_20mm.xlsx');
fprintf('SHEETS: %s\n', strjoin(sheetnames(f), ' | '));
for s = 1:numel(sheetnames(f))
    nm = sheetnames(f);
    if ~contains(nm{s}, 'Ext') && ~contains(nm{s}, 'ext'), continue; end
    T = readtable(f, 'Sheet', nm{s}, 'VariableNamingRule', 'preserve');
    fprintf('--- SHEET %s (%d rows) ---\n', nm{s}, height(T));
    disp(T);
end

% per-run Stats from the mats
d = dir(fullfile(base, 'Ext_20mm', 'ExtTest_2_*.mat'));
for i = 1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name));
    if isfield(S, 'Stats')
        fprintf('%s Stats:\n', d(i).name);
        disp(S.Stats);
    else
        fprintf('%s fields: %s\n', d(i).name, strjoin(fieldnames(S), ', '));
    end
end
