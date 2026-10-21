% gain2_summary_20261020.m -- read the gain-2 re-optimization outputs and
% print the quotable comparison vs the fixed-gain [2 1] sweep.
root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
sites = {'med_gas_r', 'soleus_r', 'tib_ant_r', 'rect_fem_r', 'bifemsh_r'};
for s = 1:numel(sites)
    f = dir(fullfile(resDir, sprintf('BiPulley_opensim_%s_gain2_*.mat', sites{s})));
    if isempty(f)
        fprintf('%s: NO gain2 mat\n', sites{s});
        continue;
    end
    S = load(fullfile(resDir, f(1).name), 'spec', 'finalObj', 'st', 'xBest', 'exitflagS', 'exitflagP');
    aK = S.finalObj.ActiveCrossing;
    tau = S.finalObj.Torque_ins(:, :, :, aK);
    tn = squeeze(vecnorm(tau, 2, 2));
    axk = S.spec.kin.axes(aK, :) / norm(S.spec.kin.axes(aK, :));
    tv = squeeze(tau);
    tax = reshape(sum(reshape(tv, [], 3) .* axk, 2), size(tn));
    fprintf(['%s (%s):\n  xBest rest=%.4f tendon=%.4f | exit S=%g P=%g | margin=%+.4f' ...
        ' worstC=%+.3g\n  peak|tau_ax|=%.2f N*m  infeas=%d slack=%d | Tackle=%s gain=%s mode=%s\n'], ...
        sites{s}, f(1).name, S.xBest(1), S.xBest(2), S.exitflagS, S.exitflagP, ...
        S.st.margin, S.st.worstC, max(abs(tax(:))), ...
        sum(S.finalObj.PulleyInfeasible(:)), sum(S.finalObj.PulleySlack(:)), ...
        mat2str(S.finalObj.TackleLineParts), mat2str(S.finalObj.PulleyGain), ...
        strjoin(cellfun(@char, S.finalObj.RoutingMode, 'UniformOutput', false), '/'));
end
fprintf('CSV rows:\n');
csv = fileread(fullfile(resDir, 'BiPulley_batch_summary_gain2_20261020.csv'));
fprintf('%s\n', csv);
