% qc_sweeps_20261020.m -- read-only QC on the sweep outputs + inventory of
% the four route-variant mats for the Table-5 comparison figure.
root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
gl = {'1/1', '2/1', '3/1'}; pl = [300 450 620];

fprintf('--- flag patterns on the margin slice ---\n');
sn = {'soleus_r', 'tib_ant_r', 'bifemsh_r', 'med_gas_r', 'rect_fem_r'};
for s = 1:numel(sn)
    D = load(fullfile(resDir, sprintf('PulleySurfaces_%s_20261020.mat', sn{s})));
    fprintf('%s: theta1 slice [%.1f %.1f] deg\n', sn{s}, ...
        D.theta1Slice_deg(1), D.theta1Slice_deg(end));
    for g = 1:3
        for pp = 1:3
            ninf = sum(D.infeasSlice(:, g, pp));
            nsl = sum(D.slackSlice(:, g, pp));
            nzero = sum(abs(D.tauSlice(:, g, pp)) < 1e-9);
            npos = sum(abs(D.tauSlice(:, g, pp)) >= 1e-9 & ~isnan(D.tauSlice(:, g, pp)));
            fprintf('  G=%s p=%d: taut %2d  slack %2d  infeas %2d  zeroTau %2d\n', ...
                gl{g}, pl(pp), npos, nsl, ninf, nzero);
        end
    end
    % sign check vs requirement
    fprintf('  sign robot tau (taut cells): %+d | sign tau_req: %+d\n', ...
        sign(sum(sign(D.tauSlice(abs(D.tauSlice) > 1e-9)))), ...
        sign(sum(sign(D.tauReqSlice))));
end

fprintf('\n--- med_gas dense grid alive map (G=1/1, p=620), rows=knee [%.0f..%.0f], cols=ankle [%.0f..%.0f] ---\n', ...
    D2minmax(resDir));
fprintf('\n--- four route-variant mats (Table 5 figure inputs) ---\n');
vv = {'Vas_Pam_20mm_Result_20260925.mat', 'Vas_Pam_20mm_Result_cam_20261014_0416.mat', ...
    'Bifemsh_20mm_Result_20260921_2037.mat', 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'};
for v = 1:numel(vv)
    V = load(fullfile(resDir, vv{v}));
    fprintf('\n== %s ==\n', vv{v});
    fn = fieldnames(V);
    tgt = '';
    if isfield(V, 'predBest'); tgt = 'predBest';
    elseif isfield(V, 'best'); tgt = 'best'; end
    if ~isempty(tgt)
        fprintf('%s fields: %s\n', tgt, strjoin(fieldnames(V.(tgt))', ' '));
        P = V.(tgt);
        for c = {'TorqueZ', 'TorqueInsZ', 'Torque', 'phiD', 'phi'}
            if isfield(P, c{1})
                x = P.(c{1});
                if isnumeric(x)
                    fprintf('  %s: [%s] n=%d\n', c{1}, num2str([min(x(:)) max(x(:))], '%.3f '), numel(x));
                end
            end
        end
    end
    if isfield(V, 'ctx')
        c2 = fieldnames(V.ctx);
        for c = {'humanAngleD', 'humanTorque', 'humanTorqueAbs', 'phiD', 'Pbins'}
            if isfield(V.ctx, c{1})
                x = V.ctx.(c{1});
                if isnumeric(x)
                    fprintf('  ctx.%s: [%s] n=%d\n', c{1}, num2str([min(x(:)) max(x(:))], '%.3f '), numel(x));
                end
            end
        end
    end
    if isfield(V, 'ctxCamBest'); fprintf('  has ctxCamBest: %s\n', strjoin(fieldnames(V.ctxCamBest)', ' ')); end
    if isfield(V, 'predBest') && isfield(V.predBest, 'pressure'); fprintf('  predBest.pressure? no\n'); end
end
fprintf('\nQC DONE\n');

function s = D2minmax(resDir) %#ok<DEFNU>
s = '';
D = load(fullfile(resDir, 'PulleySurfaces_med_gas_r_20261020.mat'));
alive = abs(D.tauDense(:, :, 1, 1, 3)) > 1e-9;
infe = D.infeasDense(:, :, 1, 3);
fprintf('     ankle:'); fprintf('%7.1f', D.theta2Dense_deg); fprintf('\n');
for i = 1:size(alive, 1)
    fprintf('  knee %6.1f:', D.theta1Dense_deg(i));
    for j = 1:size(alive, 2)
        if infe(i, j); fprintf('      I');
        elseif alive(i, j); fprintf('%6.1f', D.tauDense(i, j, 1, 1, 3));
        else; fprintf('      .'); end
    end
    fprintf('\n');
end
end
