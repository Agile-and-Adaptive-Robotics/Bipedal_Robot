%% sns_run_beer_spiking.m — run BeerCupReflexDemo_Spiking (reflex ON vs OFF)
%
% Same protocol as sns_run_beer_demo.m (kReflex 1 / 0). Baseline numbers
% from the committed demo (2026-10-02 run, this machine):
%   max sag after 2 s: ON 2.00 deg / OFF 7.85 deg | A_bi 0.392->0.417,
%   A_tri 0.400->0.316 (ON). Saves goal2_beer_spiking.{mat,png}.

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
load_system('BeerCupReflexDemo_Spiking');

runs = struct('label', {'reflex ON', 'reflex OFF'}, 'k', {1, 0});
for r = 1:2
    assignin('base', 'kReflex', runs(r).k);
    out = sim('BeerCupReflexDemo_Spiking');
    runs(r).th   = out.log_th;   runs(r).mCup = out.log_mCup;
    runs(r).F    = out.log_F_bi; runs(r).A    = out.log_A_bi;
    runs(r).Ft   = out.log_F_tri; runs(r).At  = out.log_A_tri;
    runs(r).spk  = out.log_IN_Ia_bi_spk;
    w = runs(r).th.Time > 2;
    runs(r).sagAfter2 = max(abs(runs(r).th.Data(w)))*180/pi;
    runs(r).final = runs(r).th.Data(end)*180/pi;
    runs(r).rateIa = count_edges(runs(r).spk.Time, runs(r).spk.Data)/runs(r).th.Time(end);
    fprintf(['%s (spiking): theta final %+.2f deg, max|theta| after 2 s %.2f deg; ' ...
        'biceps F %.1f N, A_bi %.3f -> %.3f; triceps A_tri %.3f -> %.3f; Ia_bi rate %.1f Hz (m_cup %.2f kg)\n'], ...
        runs(r).label, runs(r).final, runs(r).sagAfter2, runs(r).F.Data(end), ...
        runs(r).A.Data(1), runs(r).A.Data(end), runs(r).At.Data(1), runs(r).At.Data(end), ...
        runs(r).rateIa, runs(r).mCup.Data(end));
end
assignin('base', 'kReflex', 1);
fprintf(['BASELINE: max sag after 2 s ON 2.00 / OFF 7.85 deg | A_bi 0.392->0.417, ' ...
    'A_tri 0.400->0.316 (ON) (2026-10-02 run of the committed demo)\n']);

thOn = runs(1).th; thOff = runs(2).th;
fig = figure('Visible', 'off', 'Position', [100 100 900 1150]);
subplot(5,1,1); plot(thOff.Time, thOff.Data*180/pi, 'LineWidth', 1.4); hold on;
plot(thOn.Time, thOn.Data*180/pi, 'LineWidth', 1.6); grid on;
ylabel('cup sag \theta (deg)'); legend('reflex OFF', 'reflex ON', 'Location', 'southeast');
title('SPIKING beer-cup reflex: spiking Ia interneurons + hybrid synapses hold the cup level');
subplot(5,1,2); plot(runs(1).mCup.Time, runs(1).mCup.Data, 'LineWidth', 1.4); grid on;
ylabel('beer mass (kg)');
subplot(5,1,3); plot(runs(1).A.Time, runs(1).A.Data, 'LineWidth', 1.6); hold on;
plot(runs(1).At.Time, runs(1).At.Data, 'LineWidth', 1.6, 'LineStyle', '--'); ...
plot(runs(2).A.Time, runs(2).A.Data, 'LineWidth', 1.0, 'Color', [0.6 0.6 0.6]); grid on;
ylabel('activation'); legend('biceps (ON)', 'triceps (ON)', 'biceps (OFF)');
subplot(5,1,4); plot(runs(1).F.Time, runs(1).F.Data, 'LineWidth', 1.4); hold on;
plot(runs(1).Ft.Time, runs(1).Ft.Data, 'LineWidth', 1.4, 'LineStyle', '--'); grid on;
ylabel('BPA force (N)'); legend('biceps (ON)', 'triceps (ON)');
subplot(5,1,5); plot(runs(1).spk.Time, runs(1).spk.Data, 'LineWidth', 0.6); grid on;
ylabel('Ia_{bi} spikes'); xlabel('time (s)'); ylim([-0.1 1.3]);
exportgraphics(fig, fullfile(rep, 'goal2_beer_spiking.png'), 'Resolution', 150);
save(fullfile(rep, 'goal2_beer_spiking.mat'), 'runs');
fprintf('saved %s + .png\n', fullfile(rep, 'goal2_beer_spiking.mat'));

function n = count_edges(tt, dd)
v = dd(:) >= 0.5;
n = sum(v & ~[false; v(1:end-1)]);
end
