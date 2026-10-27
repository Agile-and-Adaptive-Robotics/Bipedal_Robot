% probe_mats_v3_20261027.m -- walk every v3 candidate mat and print the
% exact meetCells comparison (recomputed vs recorded) BEFORE asserting;
% also every surveyMeta field's class/size for the first failure.
function probe_mats_v3_20261027()
R = 'D:\GitHub\Bipedal_Robot';
capDir = fullfile(R, 'tmp', 'paper_bridge', 'capability_table_20261022');
resDir = fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'BiPulley_surveys');
Tpm = readtable(fullfile(capDir, 'BiPulley_survey2_perMuscle_20261023.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
T3 = readtable(fullfile(capDir, 'BiPulley_survey3_orig_4x20_20261027.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
nbad = 0;
for k = 1:height(Tpm)
    site = strtrim(Tpm.site{k});
    mats = {};
    tags = {};
    ix = find(strcmp(T3.site, site), 1);
    if ~isempty(ix)
        mats{end + 1} = strtrim(T3.matFile{ix}); %#ok<AGROW>
        tags{end + 1} = 'O4x20'; %#ok<AGROW>
    else
        d = dir(fullfile(resDir, ['BiPulley_survey2_', site, '_nB1_G1_*.mat']));
        [~, im] = max([d.datenum]);
        mats{end + 1} = fullfile(d(im).folder, d(im).name); %#ok<AGROW>
        tags{end + 1} = 'Orung1'; %#ok<AGROW>
    end
    mats{end + 1} = strtrim(Tpm.matFile{k});
    tags{end + 1} = 'UPD';
    for mi = 1:numel(mats)
        try
            S = load(mats{mi}, 'surveyMeta', 'tauHuman', 'dom', 'g');
            tauH = S.tauHuman;
            dom = logical(S.dom);
            tauR = S.g.tauB;
            taut = ~logical(S.g.slack) & ~logical(S.g.infeas) & isfinite(tauR);
            meet = taut & isfinite(tauH) & ...
                (sign(tauR) == sign(tauH) | abs(tauH) <= 1e-9) & ...
                abs(tauR) >= abs(tauH) - 1e-12;
            meetCells = sum(meet(:) & dom(:));
            rec = S.surveyMeta.meetCells;
            fprintf('%-12s %-6s recomputed %3d  recorded %s  class %s %s  eq=%d\n', ...
                site, tags{mi}, meetCells, mat2str(rec), class(rec), ...
                mat2str(size(rec)), meetCells == double(rec));
            if ~(isscalar(rec) && isfinite(rec) && isnumeric(rec) && ...
                    meetCells == double(rec))
                nbad = nbad + 1;
                fn = fieldnames(S.surveyMeta);
                for q = 1:numel(fn)
                    v = S.surveyMeta.(fn{q});
                    fprintf('    field %-18s class %-10s size %s\n', fn{q}, ...
                        class(v), mat2str(size(v)));
                end
            end
        catch e
            nbad = nbad + 1;
            fprintf('%-12s %-6s LOAD/RUN FAIL: %s\n', site, tags{mi}, e.message);
        end
    end
end
fprintf('PROBE DONE, problem rows: %d\n', nbad);
end
