%% sns_run_rgmn_spiking.m — run SNS_SpikingRG_MN, verify the hybrid motif
%
% PASS criteria (representative subnetwork of the spiking conversion):
%   1. BOTH half-centers burst (>= 3 bursts each in 12 s)
%   2. bursts are ANTIPHASE (rate-envelope correlation < -0.3)
%   3. the MN pool V is GRADED (modulates with the RG_E burst envelope,
%      ripple << modulation) — the spikes-to-analog conversion works
%   4. MN_F1 is antiphase with the MN_E pool

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
load_system('SNS_SpikingRG_MN');
set_param('SNS_SpikingRG_MN', 'StopTime', '12');
out = sim('SNS_SpikingRG_MN');

vE  = out.log_V_RG_E;  sE = out.log_spk_RG_E;
vF  = out.log_V_RG_F;  sF = out.log_spk_RG_F;
adp = out.log_V_Adp_E;
m1  = out.log_V_MN_E1; m2 = out.log_V_MN_E2; mF = out.log_V_MN_F1;
aE  = out.log_A_E;

% rate envelopes (1-s boxcar on the spike lines, sampled at 50 Hz)
fs = 50; tq = (0:1/fs:12)';
rE = boxrate(sE.Time, sE.Data, tq, 1.0);
rF = boxrate(sF.Time, sF.Data, tq, 1.0);
cc = corr(rE(:), rF(:));

% bursts: contiguous regions where the 1-s rate envelope > 20% of its max
bE = bursts(tq, rE);
bF = bursts(tq, rF);

% MN graded-ness: modulation depth vs ripple (detrended std in a window
% where the pool is driven: MN V above its 25th percentile)
mn = interp1(m1.Time, m1.Data, tq, 'linear');
drv = mn > prctile(mn, 25);
ripple = std(diff(mn(drv))) / max(std(mn), eps);
mod_depth = max(mn) - min(mn);
ccE = corr(rE(:), mn(:));
mFq = interp1(mF.Time, mF.Data, tq, 'linear');
ccF = corr(rF(:), mFq(:));

fprintf(['RG->MN subnetwork: RG_E %d bursts / RG_F %d bursts (12 s); envelope corr(E,F) %+.2f; ' ...
    'MN_E1 V range %.1f..%.1f mV (depth %.1f mV, corr with RG_E rate %+.2f); ' ...
    'MN_F1 corr with RG_F rate %+.2f; activation range %.2f..%.2f\n'], ...
    bE.n, bF.n, cc, min(mn), max(mn), mod_depth, ccE, ccF, min(aE.Data), max(aE.Data));
ok = bE.n >= 3 && bF.n >= 3 && cc < -0.3 && ccE > 0.5 && ccF > 0.5 && mod_depth > 5;
if ok, fprintf('RG->MN SUBNETWORK: PASS\n'); else, fprintf('RG->MN SUBNETWORK: CHECK (criteria above)\n'); end

fig = figure('Visible', 'off', 'Position', [100 100 900 1100]);
subplot(6,1,1); plot(tq, rE, 'LineWidth', 1.2); hold on; plot(tq, rF, 'LineWidth', 1.2); grid on;
ylabel('RG rate (Hz)'); legend('RG_E', 'RG_F');
title('Spiking RG half-centers -> hybrid synapses -> non-spiking MN pool');
subplot(6,1,2); plot(vE.Time, vE.Data, 'LineWidth', 0.6); grid on; ylabel('RG_E V (mV)');
subplot(6,1,3); plot(tq, mn, 'LineWidth', 1.2); grid on; ylabel('MN_{E1} V (mV)');
subplot(6,1,4); plot(tq, interp1(mF.Time, mF.Data, tq, 'linear'), 'LineWidth', 1.2); grid on;
ylabel('MN_{F1} V (mV)');
subplot(6,1,5); plot(adp.Time, adp.Data, 'LineWidth', 1.2); grid on; ylabel('Adp_E V (mV)');
subplot(6,1,6); plot(aE.Time, aE.Data, 'LineWidth', 1.2); grid on;
ylabel('activation'); xlabel('time (s)');
exportgraphics(fig, fullfile(rep, 'goal2_rgmn_spiking.png'), 'Resolution', 150);
save(fullfile(rep, 'goal2_rgmn_spiking.mat'), 'out', 'rE', 'rF', 'tq', 'mn', 'bE', 'bF', 'ok');
fprintf('saved %s + .png\n', fullfile(rep, 'goal2_rgmn_spiking.mat'));

function r = boxrate(tt, dd, tq, win)
% spikes/s in a trailing boxcar of width win, sampled on tq
v = dd(:) >= 0.5;
ts = tt(v);
r = zeros(size(tq));
for k = 1:numel(tq)
    r(k) = sum(ts > tq(k) - win & ts <= tq(k)) / win;
end
end

function b = bursts(tq, r)
thr = 0.2*max(r);
on = r > thr;
d = diff([false; on(:); false]);
i0 = find(d == 1); i1 = find(d == -1) - 1;
b = struct('n', numel(i0), 'start', tq(i0), 'stop', tq(i1));
end

function p = prctile(x, q)
p = quantile(x, q/100);
end
