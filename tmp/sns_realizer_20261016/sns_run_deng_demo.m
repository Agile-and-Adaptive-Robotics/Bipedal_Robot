function sns_run_deng_demo()
% Run SNS_Deng_CPGDemo (one 10 nA / 1 ms pulse at t = 0.1 s -> 20 s) and
% verify continuous alternation + agreement with the numpy ODE reference
% (spinal/deng_cpg_ode.py -> results\deng_cpg_ref.mat, same equations,
% same dt = 0.1 ms Euler).

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                     % SNS_Simscape (SNS_Library)
mdl = 'SNS_Deng_CPGDemo';
load_system(fullfile(here, [mdl '.slx']));
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'sigs');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');

get = @(nm) out.sigs.get(nm).Values;
vE = get('V_RG_ext'); vF = get('V_RG_flx');
tE = vE.Time(:); yE = vE.Data(:); yF = vF.Data(:);

% cycle structure: peaks of ext after the kick
iPk = find(islocalmax(yE, 'MinProminence', 1));
tp = tE(iPk);
tp = tp(tp > 0.3);
per = mean(diff(tp));
fprintf(['DENG CPG DEMO: %d RG_ext bursts over 20 s from ONE pulse; ' ...
         'period %.3f s; V range [%.1f, %.1f] mV\n'], numel(tp), per, ...
        min(min(yE), min(yF)), max(max(yE), max(yF)));

% antiphase check: correlation of ext vs flx (base-MATLAB Pearson — the
% laptop license has no Statistics Toolbox, corr() is not defined there)
za = yE - mean(yE); zb = yF - mean(yF);
r = sum(za.*zb) / sqrt(sum(za.^2) * sum(zb.^2));
fprintf('RG ext-vs-flx correlation: %.3f (negative = alternation)\n', r);

% compare with the numpy reference
ref = load(fullfile(here, '..', 'results', 'deng_cpg_ref.mat'));
vnE = interp1(ref.t(:), ref.traces(:, 1), tE, 'linear', 'extrap');
vnF = interp1(ref.t(:), ref.traces(:, 2), tE, 'linear', 'extrap');
dE = max(abs(yE - vnE)); dF = max(abs(yF - vnF));
fprintf('vs numpy ODE reference: max|dV_ext| %.3f mV, max|dV_flx| %.3f mV\n', dE, dF);

% figure
fig = figure('Visible', 'off', 'Position', [80 80 1000 780]);
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(tl); hold on;
plot(tE, yE, 'LineWidth', 1.1); plot(tE, yF, 'LineWidth', 1.1);
grid on; ylabel('RG V (mV)'); legend('HC\_ext', 'HC\_flx');
title(sprintf('Deng 2019 two-layer CPG: one 10 nA/1 ms pulse at t=0.1 s (period %.2f s)', per));
nexttile(tl); hold on;
ph = get('V_PF_hip_e'); plot(ph.Time(:), ph.Data(:), 'LineWidth', 1.1);
ph = get('V_PF_hip_f'); plot(ph.Time(:), ph.Data(:), 'LineWidth', 1.1);
grid on; ylabel('PF hip V (mV)'); legend('HC\_ext', 'HC\_flx');
nexttile(tl); hold on;
ph = get('V_PF_ka_e'); plot(ph.Time(:), ph.Data(:), 'LineWidth', 1.1);
ph = get('V_PF_ka_f'); plot(ph.Time(:), ph.Data(:), 'LineWidth', 1.1);
grid on; ylabel('PF knee/ankle V (mV)'); legend('HC\_ext', 'HC\_flx');
nexttile(tl); hold on;
ph = get('act_ext'); plot(ph.Time(:), ph.Data(:), 'LineWidth', 1.1);
ph = get('act_flx'); plot(ph.Time(:), ph.Data(:), 'LineWidth', 1.1);
grid on; ylabel('hip activation'); xlabel('t (s)'); legend('ext', 'flx');
exportgraphics(fig, fullfile(here, 'sns_deng_demo.png'), 'Resolution', 130);
close_system(mdl, 0);
fprintf('saved %s\n', fullfile(here, 'sns_deng_demo.png'));
end
