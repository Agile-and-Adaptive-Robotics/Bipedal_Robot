% Focused dump: ExtTest20mm_2 sheet as label:value rows + per-run Stats
base = 'D:/GitHub/Bipedal_Robot/Testing_Data/2026_06_Festo';
f = fullfile(base, 'Results_table_20mm.xlsx');
R = readcell(f, 'Sheet', 'ExtTest20mm_2');
fprintf('--- ExtTest20mm_2 (%d x %d) ---\n', size(R, 1), size(R, 2));
for i = 1:size(R, 1)
    vals = '';
    for j = 3:min(size(R, 2), 26)
        v = R{i, j};
        if isnumeric(v) && ~isnan(v)
            vals = [vals sprintf('%9.4g', v)]; %#ok<AGROW>
        end
    end
    lab = R{i, 2}; if ~ischar(lab) && ~isstring(lab), lab = ''; end
    fprintf('%-28s |%s\n', char(lab), vals);
end
d = dir(fullfile(base, 'Ext_20mm', 'ExtTest_2_*.mat'));
for i = 1:numel(d)
    S = load(fullfile(d(i).folder, d(i).name));
    fn = fieldnames(S);
    fprintf('%s fields: %s\n', d(i).name, strjoin(fn, ','));
    if isfield(S, 'Stats')
        disp(S.Stats);
    end
end
% old-tab convention check: which torque row fed the paper's first three points
R1 = readcell(f, 'Sheet', 'ExtTest20mm_1');
for i = 1:size(R1, 1)
    lab = R1{i, 2};
    if (ischar(lab) || isstring(lab)) && contains(char(lab), 'orque')
        fprintf('TAB1 ROW [%s]:', char(lab));
        for j = 3:6
            v = R1{i, j};
            if isnumeric(v) && ~isnan(v), fprintf(' %g', v); end
        end
        fprintf('\n');
    end
end
