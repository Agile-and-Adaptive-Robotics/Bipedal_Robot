% inspect_for_sweeps_20261020.m -- read-only inventory for the section-5
% pulley sweeps (EDIT_SPEC_20261020 §5). Scratch; writes nothing.
root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
sites = {'med_gas_r', 'soleus_r', 'tib_ant_r', 'rect_fem_r', 'bifemsh_r'};

for s = 1:numel(sites)
    f = dir(fullfile(resDir, sprintf('BiPulley_opensim_%s_*.mat', sites{s})));
    S = load(fullfile(resDir, f(1).name));
    fprintf('\n=== %s (%s) ===\n', sites{s}, f(1).name);
    fprintf('vars: %s\n', strjoin(fieldnames(S)', ' '));
    sp = S.spec;
    fprintf('frames: %s | cross=[%d %d] | pts=%d | fmax=%g N (%s)\n', ...
        strjoin(sp.frames, ','), sp.cross(1), sp.cross(2), size(sp.points, 1), ...
        sp.fmax, sp.fmaxSource);
    fprintf('kin.types: %s %s | ranges: [%g %g; %g %g] | kneeSlot=%d\n', ...
        sp.kin.types{1}, sp.kin.types{2}, sp.kin.thetaRangesDeg(1, :), ...
        sp.kin.thetaRangesDeg(2, :), sp.kin.kneeSlot);
    fprintf('kin.axes: [%.3f %.3f %.3f] [%.3f %.3f %.3f] | pivots: [%.3f %.3f %.3f] [%.3f %.3f %.3f]\n', ...
        sp.kin.axes(1, :), sp.kin.axes(2, :), sp.kin.pivots(1, :), sp.kin.pivots(2, :));
    fprintf('pulley: nPulleyBPA=%s tackle=%s mode=%s\n', mat2str(sp.pulley.nPulleyBPA), ...
        mat2str(sp.pulley.tackleLineParts), joinMode(sp.pulley.routingMode));
    if isfield(sp.kin, 'pointFolds'); fprintf('pointFolds{2}: %s\n', mat2str(sp.kin.pointFolds{2})); end
    fprintf('notes: %s\n', sp.notes);
    fprintf('T size: %s | theta1 deg: %s | theta2 deg: %s\n', mat2str(size(S.T)), ...
        mat2str(round(S.theta1 * 180 / pi, 2)), mat2str(round(S.theta2 * 180 / pi, 2)));
    fprintf('xBest: %s | finalObj RestingL=%.4f Kmax=%.4f TendonL=%.4f Pres=%g\n', ...
        mat2str(S.xBest(:)'), S.finalObj.RestingL, S.finalObj.Kmax, ...
        S.finalObj.TendonL, S.finalObj.Pressure);
    fprintf('finalObj RoutingMode: %s | ActiveCrossing=%d\n', ...
        joinMode(S.finalObj.RoutingMode), S.finalObj.ActiveCrossing);
    % Recorded torque stats per crossing (norm about hinge axis approximation:
    % report vector norm and dominant component)
    for k = 1:2
        tq = S.finalObj.Torque_ins(:, :, :, k);
        tn = squeeze(vecnorm(tq, 2, 2));   % N1 x N2
        fprintf('crossing %d: |tau| min %.3f max %.3f  nNaN %d  nZero %d\n', k, ...
            min(tn(:)), max(tn(:)), sum(isnan(tn(:))), sum(tn(:) == 0));
    end
    % base Xi constants saved?
    fprintf('baseXi0=%g baseXi1=%g baseXi2=%g baseWraps=%g basePressure=%g\n', ...
        S.baseXi0, S.baseXi1, S.baseXi2, S.baseWraps, S.basePressure);
end

% Mono pulley bifemsh mat
M = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'));
fprintf('\n=== mono bifemsh pulley ===\n');
fprintf('vars: %s\n', strjoin(fieldnames(M)', ' '));
fprintf('ctx humanAngleD: [%g %g] n=%d | humanTorque: [%.3f %.3f]\n', ...
    M.ctx.humanAngleD(1), M.ctx.humanAngleD(end), numel(M.ctx.humanAngleD), ...
    M.ctx.humanTorque(1), M.ctx.humanTorque(end));
if isfield(M, 'predBest')
    fprintf('predBest fields: %s\n', strjoin(fieldnames(M.predBest)', ' '));
    if isfield(M.predBest, 'TorqueInsZ')
        fprintf('TorqueInsZ: min %.3f max %.3f n=%d | phiD range [%g %g]\n', ...
            min(M.predBest.TorqueInsZ), max(M.predBest.TorqueInsZ), ...
            numel(M.predBest.TorqueInsZ), M.ctx.phiD(1), M.ctx.phiD(end));
    end
end

% Route-variant mats for the four-variant figure (Table 5)
cands = {'Vas_Pam_20mm_Result_20260925.mat', 'Vas_Pam_20mm_Result_cam_20261014_0416.mat', ...
    'Bifemsh_20mm_Result_20260921_2037.mat', 'Bifemsh_20mm_Result_20260920_1443.mat'};
for c = 1:numel(cands)
    pth = fullfile(resDir, cands{c});
    if ~isfile(pth); fprintf('\n--- %s MISSING\n', cands{c}); continue; end
    V = load(pth);
    fprintf('\n=== variant %s ===\n', cands{c});
    fprintf('vars: %s\n', strjoin(fieldnames(V)', ' '));
    if isfield(V, 'best')
        fprintf('best fields: %s\n', strjoin(fieldnames(V.best)', ' '));
        if isfield(V.best, 'TorqueZ'); fprintf('best.TorqueZ: [%.3f %.3f] n=%d\n', min(V.best.TorqueZ), max(V.best.TorqueZ), numel(V.best.TorqueZ)); end
        if isfield(V.best, 'TorqueInsZ'); fprintf('best.TorqueInsZ: [%.3f %.3f] n=%d\n', min(V.best.TorqueInsZ), max(V.best.TorqueInsZ), numel(V.best.TorqueInsZ)); end
    end
    if isfield(V, 'ctx') && isfield(V.ctx, 'humanTorqueAbs')
        fprintf('ctx.humanTorqueAbs: [%.3f %.3f] | req = [%.3f %.3f]?\n', ...
            min(V.ctx.humanTorqueAbs), max(V.ctx.humanTorqueAbs), ...
            min(1.05 * V.ctx.humanTorqueAbs), max(1.05 * V.ctx.humanTorqueAbs));
    end
end
fprintf('\nDONE\n');

function out = joinMode(v)
if ischar(v)
    out = v;
elseif iscell(v)
    out = strjoin(cellfun(@char, v, 'UniformOutput', false), '/');
else
    out = mat2str(v);
end
end
