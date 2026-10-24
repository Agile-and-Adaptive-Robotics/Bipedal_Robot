% probe_ext20_20261023.m -- read-only inspection for the extensor-20mm
% display-item edit: cam mat fields, ExtTest20mm_2 xlsx tab, ExtTest_2 mats.
root = 'D:/GitHub/Bipedal_Robot';

%% 1. cam mat of record fields
C = load(fullfile(root, 'Code/Matlab/Mesh_Optimization/Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
fprintf('CAM MAT top-level:\n'); disp(fieldnames(C).');
if isfield(C, 'ctxCamBest')
    ctx = C.ctxCamBest;
    fprintf('ctxCamBest fields:\n'); disp(fieldnames(ctx).');
    fprintf('phiD: n=%d range [%.2f %.2f] deg\n', numel(ctx.phiD), min(ctx.phiD), max(ctx.phiD));
    if isfield(ctx, 'humanAngleD')
        fprintf('humanAngleD: n=%d range [%.2f %.2f]\n', numel(ctx.humanAngleD), min(ctx.humanAngleD), max(ctx.humanAngleD));
        fprintf('humanTorque range [%.2f %.2f] N*m\n', min(ctx.humanTorque), max(ctx.humanTorque));
    end
    for f = {'Xi0','Xi1','Xi2','Xi3','Dia','BPAcount','targetPressure','N','pos'}
        if isfield(ctx, f{1})
            fprintf('ctx.%s = %s\n', f{1}, mat2str(ctx.(f{1})));
        end
    end
    if isfield(ctx, 'geo')
        g = ctx.geo; fn = fieldnames(g);
        fprintf('geo fields: %s\n', strjoin(fn.', ', '));
        for f = {'seedP6XY','seedP7XY','seedP8XY','tibiaWrapCircumcircle'}
            if isfield(g, f{1})
                fprintf('geo.%s = %s\n', f{1}, mat2str(g.(f{1})));
            end
        end
    end
end
if isfield(C, 'predBest')
    fprintf('predBest fields:\n'); disp(fieldnames(C.predBest).');
    if isfield(C.predBest, 'TorqueZ')
        tz = C.predBest.TorqueZ(:);
        fprintf('predBest.TorqueZ: n=%d finite=%d range [%.2f %.2f] N*m\n', ...
            numel(tz), sum(isfinite(tz)), min(tz), max(tz));
    end
end
if isfield(C, 'xBest'), fprintf('xBest = %s\n', mat2str(C.xBest, 8)); end
if isfield(C, 'XiUsed'), fprintf('XiUsed = %s\n', mat2str(C.XiUsed, 8)); end

%% 2. ExtTest20mm_2 tab of Results_table_20mm.xlsx
fprintf('\n=== ExtTest20mm_2 tab (first 25 rows x 26 cols) ===\n');
T = readcell(fullfile(root, 'Testing_Data/2026_06_Festo/Results_table_20mm.xlsx'), ...
    'Sheet', 'ExtTest20mm_2');
for r = 1:min(size(T,1), 25)
    rowstr = '';
    for c = 1:min(size(T,2), 26)
        v = T{r,c};
        if isnumeric(v) && isfinite(v)
            rowstr = [rowstr sprintf('%9.3g', v)]; %#ok<AGROW>
        elseif ischar(v)
            rowstr = [rowstr sprintf(' %8.8s', v)]; %#ok<AGROW>
        else
            rowstr = [rowstr '         ']; %#ok<AGROW>
        end
    end
    fprintf('r%02d:%s\n', r, rowstr);
end

%% 3. ExtTest_2 mats: last rows (held values)
fprintf('\n=== ExtTest_2 mats (last row + guard cols) ===\n');
d = dir(fullfile(root, 'Testing_Data/2026_06_Festo/Ext_20mm', 'ExtTest_2_*.mat'));
for k = 1:numel(d)
    S = load(fullfile(d(k).folder, d(k).name));
    if isfield(S, 'data')
        M = S.data;
        fprintf('%s: %d rows; last row: t=%.0f ms F=%.1f N P=%.1f kPa fill=%d exh=%d | K=%.1f LC=%.1f tz=%.3f tgt=%.3f\n', ...
            d(k).name, size(M,1), M(end,1), M(end,2), M(end,3), M(end,4), M(end,5), ...
            M(end,6), M(end,7), M(end,8), M(end,9));
    end
end
disp('DONE_probe_ext20');
