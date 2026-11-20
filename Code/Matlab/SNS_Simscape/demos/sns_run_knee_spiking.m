%% sns_run_knee_spiking.m — run KneeReflexDemo_Spiking, compare vs baseline
%
% Computes the SAME metrics as the committed KneeReflexDemo baseline (see
% reports_spiking_20261002/goal2_run_baselines.m, values recorded
% 2026-10-02 on this machine from the committed demo):
%   baseline: theta(0.05 s) 15.79 deg | settle mean 43.5 deg
%             (range 41.6-44.9, last 2 s) | A_ext 0.53 / A_flex 0.53
% plus the spiking-layer observables (IN firing rates, ripple).
% Saves goal2_knee_spiking.{mat,png} into the campaign report folder.

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);
addpath(fileparts(cdto));
rep = cdto;
while exist(fullfile(rep, 'Code', 'MuJoCo_SNS', 'spinal'), 'dir') == 0
    pn = fileparts(rep);
    if strcmp(pn, rep), error('repo root not found'); end
    rep = pn;
end
rep = fullfile(rep, 'Code', 'MuJoCo_SNS', 'spinal', 'reports_spiking_20261002');

load_system('SNS_Library');
load_system('KneeReflexDemo_Spiking');
set_param('KneeReflexDemo_Spiking', 'StopTime', '5');
out = sim('KneeReflexDemo_Spiking');

th  = out.log_th;    thd = out.log_thd;
Ae  = out.log_A_ext; Af  = out.log_A_flex;
Vme = out.log_V_MN_ext; Vmf = out.log_V_MN_flex;
Fe  = out.log_F_ext; Ff  = out.log_F_flex;
spkIa = out.log_IN_Ia_ext_spk; spkIb = out.log_IN_Ib_ext_spk;

t = th.Time; tf = t(end);
w = t > tf - 2;
S = struct();
S.rise   = interp1(t, th.Data, 0.05, 'linear')*180/pi;
S.settleMean = mean(th.Data(w))*180/pi;
S.settleMin  = min(th.Data(w))*180/pi;
S.settleMax  = max(th.Data(w))*180/pi;
S.final  = th.Data(end)*180/pi;
S.AeWin  = mean(Ae.Data(Ae.Time > tf-2));
S.AfWin  = mean(Af.Data(Af.Time > tf-2));
S.rateIa = count_edges(spkIa.Time, spkIa.Data)/tf;
S.rateIb = count_edges(spkIb.Time, spkIb.Data)/tf;
S.VmeMean = mean(Vme.Data(Vme.Time > tf-2));
S.VmfMean = mean(Vmf.Data(Vmf.Time > tf-2));
% antagonist alternation: dominant frequency of A_ext in the last 3 s
wl = t > tf-3;
[pxx, fx] = pwelch(Ae.Data(wl) - mean(Ae.Data(wl)), 2048, 1024, 2048, ...
    1/mean(diff(t(wl))));
[~, im] = max(pxx(fx <= 20));
S.altHz = fx(im);
fprintf(['SPIKING knee: theta(0.05s) %.2f deg | settle mean %.1f deg (range %.1f-%.1f) | ' ...
    'final %.2f deg | A_ext %.2f A_flex %.2f | IN rates Ia %.1f Hz Ib %.1f Hz | alt %.2f Hz\n'], ...
    S.rise, S.settleMean, S.settleMin, S.settleMax, S.final, S.AeWin, S.AfWin, ...
    S.rateIa, S.rateIb, S.altHz);
fprintf(['BASELINE knee: theta(0.05s) 15.79 deg | settle mean 43.5 deg (range 41.6-44.9) | ' ...
    'A_ext 0.53 A_flex 0.53 (2026-10-02 run of the committed demo)\n']);

fig = figure('Visible', 'off', 'Position', [100 100 900 1150]);
subplot(6,1,1); plot(t, th.Data*180/pi, 'LineWidth', 1.2); grid on;
ylabel('knee \theta (deg)');
title('SPIKING knee reflex — spiking interneurons + hybrid synapses -> non-spiking MNs -> BPAs');
subplot(6,1,2); plot(Vme.Time, Vme.Data, 'LineWidth', 1.2); hold on;
plot(Vmf.Time, Vmf.Data, 'LineWidth', 1.2); grid on;
ylabel('MN membrane (mV)'); legend('MN_{ext}', 'MN_{flex}'); ylim([-60 -30]);
subplot(6,1,3); plot(Ae.Time, Ae.Data, 'LineWidth', 1.2); hold on;
plot(Af.Time, Af.Data, 'LineWidth', 1.2); grid on;
ylabel('activation'); legend('extensor', 'flexor'); ylim([0 1.05]);
subplot(6,1,4); plot(Fe.Time, Fe.Data, 'LineWidth', 1.2); hold on;
plot(Ff.Time, Ff.Data, 'LineWidth', 1.2); grid on;
ylabel('BPA force (N)'); legend('extensor', 'flexor');
subplot(6,1,5); plot(spkIa.Time, spkIa.Data, 'LineWidth', 0.8); grid on;
ylabel('Ia_{ext} spikes'); ylim([-0.1 1.3]);
subplot(6,1,6); plot(spkIb.Time, spkIb.Data, 'LineWidth', 0.8); grid on;
ylabel('Ib_{ext} spikes'); xlabel('time (s)'); ylim([-0.1 1.3]);
exportgraphics(fig, fullfile(rep, 'goal2_knee_spiking.png'), 'Resolution', 150);
save(fullfile(rep, 'goal2_knee_spiking.mat'), 'out', 'S');
fprintf('saved %s + .png\n', fullfile(rep, 'goal2_knee_spiking.mat'));

function n = count_edges(tt, dd)
v = dd(:) >= 0.5;
n = sum(v & ~[false; v(1:end-1)]);
end
