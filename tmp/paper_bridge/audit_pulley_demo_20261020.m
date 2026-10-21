% audit_pulley_demo_20261020.m -- READ-ONLY part 2 of the pulley audit.
% Corrected per-crossing torque stats, budget echo, and an IN-MEMORY
% (nothing saved) fixed-gain / routing-mode sweep on the soleus mat to prove
% the tackle and bowden machinery executes and to time one construction.

root = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');

%% A. Corrected per-crossing torque + budgets
files = { ...
    'BiPulley_opensim_med_gas_r_20260926_1655.mat'; ...
    'BiPulley_opensim_soleus_r_20260926_1657.mat'; ...
    'BiPulley_opensim_tib_ant_r_20260926_1659.mat'; ...
    'BiPulley_opensim_rect_fem_r_20260926_1701.mat'; ...
    'BiPulley_opensim_bifemsh_r_20260926_1702.mat'};
fprintf('\n===== A. corrected per-crossing stats =====\n');
for k = 1:numel(files)
    S = load(fullfile(resDir, files{k}));
    fo = S.finalObj;
    aK = fo.ActiveCrossing;
    fprintf('--- %s (budgetS=%g, budgetP=%g, useParallel=%d) ---\n', ...
        S.spec.name, S.budgetS, S.budgetP, S.useParallel);
    fprintf('  theta1 (%s) = %s deg\n', sprintf('%.0f', S.N1), ...
        mat2str(S.theta1(:).'*180/pi, 4));
    fprintf('  theta2 (%s) = %s deg\n', sprintf('%.0f', S.N2), ...
        mat2str(S.theta2(:).'*180/pi, 4));
    for kk = 1:2
        tauK = squeeze(abs(fo.Torque_ins(:, :, :, kk)));   % N1 x N2
        fprintf(['  crossing %d: |Torque_ins| over 5x5 grid: min=%.4g ' ...
            'max=%.4g N*m ; NaN cells=%d ; infeasible=%d ; slack=%d\n'], ...
            kk, min(tauK(:)), max(tauK(:)), sum(isnan(tauK(:))), ...
            sum(fo.PulleyInfeasible(:)), sum(fo.PulleySlack(:)));
    end
    fprintf('  tauTarget = %g N*m (fmax %g N x 0.05 m)\n', S.tauTarget, S.spec.fmax);
end

%% B. In-memory sweep on soleus: gain and routing mode (NOTHING WRITTEN)
% NOTE (audit finding): scalar tackleLineParts broadcasts to BOTH crossings
% and trips BiPam_pulley:TwoTackles when >1 -- the CSV '1/1' notation is the
% PER-CROSSING gain pair [g1 g2]. Sweep tackle at crossing 1: [G 1].
fprintf('\n===== B. evaluation-only sweep demo (soleus mat) =====\n');
S = load(fullfile(resDir, 'BiPulley_opensim_soleus_r_20260926_1657.mat'));
sp = S.spec; fo = S.finalObj;
base = struct('name', sp.name, 'points', S.points, 'cross', sp.cross, ...
    'T', S.T, 'rest', fo.RestingL, 'kmax', fo.Kmax, 'tendon', fo.TendonL);
modes = {'moving_via', 'bowden'};
gains = [1, 2, 3];
for mi = 1:numel(modes)
    for gi = 1:numel(gains)
        gpair = [gains(gi), 1];
        tic;
        obj = BiPam_pulley(base.name, base.points, base.cross, 20, base.T, ...
            base.rest, base.kmax, base.tendon, 0.025, 620, ...
            0.005, 5e5, 1e4, 6, struct('nPulleyBPA', 1, ...
            'tackleLineParts', gpair, 'routingMode', modes{mi}));
        dt = toc;
        aK = obj.ActiveCrossing;
        tau = squeeze(abs(obj.Torque_ins(:, :, :, aK)));
        fprintf(['mode=%-10s gain=[%d %d]  PulleyActive=%s activeK=%d  ' ...
            'max|tau|=%7.3f N*m  infeas=%d slack=%d  build=%.2f s\n'], ...
            modes{mi}, gpair, mat2str(obj.PulleyActive), aK, max(tau(:)), ...
            sum(obj.PulleyInfeasible(:)), sum(obj.PulleySlack(:)), dt);
    end
end

%% C. same demo on bifemsh + med_gas at G=1,2 (timing only)
fprintf('\n===== C. timing on bifemsh_r + med_gas_r (evaluation-only) =====\n');
for f = {'BiPulley_opensim_bifemsh_r_20260926_1702.mat', ...
         'BiPulley_opensim_med_gas_r_20260926_1655.mat'}
    S2 = load(fullfile(resDir, f{1}));
    sp2 = S2.spec; fo2 = S2.finalObj;
    for G = [1, 2]
        tic;
        obj = BiPam_pulley(sp2.name, S2.points, sp2.cross, 20, S2.T, ...
            fo2.RestingL, fo2.Kmax, fo2.TendonL, 0.025, 620, ...
            0.005, 5e5, 1e4, 6, struct('nPulleyBPA', 1, ...
            'tackleLineParts', [G, 1], 'routingMode', 'moving_via'));
        dt = toc;
        aK = obj.ActiveCrossing;
        tau = squeeze(abs(obj.Torque_ins(:, :, :, aK)));
        fprintf('%s G=%d: max|tau|=%7.3f infeas=%d slack=%d build=%.2f s\n', ...
            sp2.name, G, max(tau(:)), sum(obj.PulleyInfeasible(:)), ...
            sum(obj.PulleySlack(:)), dt);
    end
end

%% D. mono bifemsh mat: angle vector + margin context
fprintf('\n===== D. mono mat extras =====\n');
M = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'));
fn = fieldnames(M);
thetish = fn(~cellfun(@isempty, regexpi(fn, 'theta|orient|angle', 'once')));
for j = 1:numel(thetish)
    v = M.(thetish{j});
    if isnumeric(v) && numel(v) >= 2
        fprintf('%s: size=%s range=[%.4g %.4g]\n', thetish{j}, ...
            mat2str(size(v)), min(v(:)), max(v(:)));
    end
end
fprintf('ctx.targetPressure=%g kPa ; predBest.pulleyGain=%s ; ', ...
    M.ctx.targetPressure, mat2str(M.predBest.pulleyGain));
fprintf('predBest.infeasibleCount=%g slackCount=%g\n', ...
    M.predBest.infeasibleCount, M.predBest.slackCount);
fprintf('fBest=%g (Opt_run_pulley objective at gain 1)\n', M.fBest);
clear obj S S2 M fo fo2 sp sp2 base
fprintf('\nDemo complete. Nothing written.\n');
