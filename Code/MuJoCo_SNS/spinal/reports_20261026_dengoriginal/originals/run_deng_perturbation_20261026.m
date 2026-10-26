% Order B2: Deng 2023 perturbation protocol on the Simulink twin.
% 2 nA to the RG FLEXOR, t = 2.0-2.5 s. Measures burst periods
% before/during/after vs the unperturbed 1.938 s. Checks P1-P3 of
% PREREG_REPRO_ABLATION_20261026.md.
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root); addpath(fullfile(root, 'demos'));
mdl = 'SNS_Deng_CPGDemo';
load_system(fullfile(root, 'demos', [mdl '.slx']));
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'sigs', ...
          'StopTime', '10');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
close_system(mdl, 0);
ts = out.sigs.get('V_RG_ext').Values;
t = ts.Time(:); v = ts.Data(:);
ipk = find(islocalmax(v, 'MinProminence', 5));
tp = t(ipk);
wpre = tp(tp < 2.0);
wdur = tp(tp >= 2.0 & tp <= 2.6);
wpos = tp(tp > 2.6);
perPre = mean(diff(wpre));
perDur = mean(diff(wdur)) * (numel(wdur) >= 2);
perPos = mean(diff(wpos));
fprintf('DENG PERTURBATION (2 nA -> RG flexor, 2.0-2.5 s):\n');
fprintf('  bursts: %d pre / %d during+0.1 / %d post\n', ...
    numel(wpre), numel(wdur), numel(wpos));
fprintf('  period pre %.4f s | during %.4f s | post %.4f s (ref 1.938)\n', ...
    perPre, perDur, perPos);
p1 = numel(wdur) + numel(wpos) > 0;              % rhythm persists
p2 = abs(perDur - 1.938) / 1.938 > 0.10;         % >10% window deviation
p3 = abs(perPos - 1.938) / 1.938 < 0.05;         % post returns
fprintf('  P1 persist: %d | P2 window-shift >10%%: %d | P3 post-return: %d\n', ...
    p1, p2, p3);
save(fullfile(root, 'results', 'deng_perturbation_twin_20261026.mat'), ...
    't', 'v', 'tp', 'perPre', 'perDur', 'perPos');
fig = figure('Visible', 'off', 'Position', [60 60 1000 420]);
plot(t, v, 'b'); hold on; grid on;
xline(2.0, 'r--', '2 nA on'); xline(2.5, 'r--', 'off');
xlabel('t (s)'); ylabel('RG\_ext V (mV)');
title(sprintf('Deng 2023 perturbation on Simulink twin: period %.3f/%.3f/%.3f s (pre/dur/post)', ...
    perPre, perDur, perPos));
exportgraphics(fig, fullfile(root, 'results', 'pictures', ...
    'deng_perturbation_twin_20261026.png'), 'Resolution', 220);
close(fig);
disp('saved results\deng_perturbation_twin_20261026.mat + png');
