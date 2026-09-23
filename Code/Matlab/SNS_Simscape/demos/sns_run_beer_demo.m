%% sns_run_beer_demo.m — run BeerCupReflexDemo twice (reflex ON vs OFF), save plots
%
% The story: beer pours into the cup (0 -> 0.5 kg over 8 s). The model starts
% in equilibrium holding the EMPTY cup (activation initialized, descending
% drive sized accordingly — no startup settling transient). With the reflex
% pathway enabled (Ia stretch reflex, Ib autogenic inhibition, and reciprocal
% Ia inhibition onto the triceps MN) the biceps adds force and the triceps is
% progressively RELEASED, so the cup stays near level; with kReflex = 0 the
% arm sags under the growing load. Triceps activation DECREASES during the
% pour (antagonist inhibition) — plotted in panel 3.
cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder

load_system('SNS_Library');
load_system('BeerCupReflexDemo');

runs = struct('label', {'reflex ON', 'reflex OFF'}, 'k', {1, 0});
for r = 1:2
    assignin('base', 'kReflex', runs(r).k);   % all six reflex gmax read kReflex
    out = sim('BeerCupReflexDemo');
    runs(r).th = out.log_th;   runs(r).mCup = out.log_mCup;
    runs(r).F  = out.log_F_bi; runs(r).A = out.log_A_bi;
    runs(r).Ft = out.log_F_tri; runs(r).At = out.log_A_tri;
    w = runs(r).th.Time > 2;   % post-startup window for honest stats
    fprintf(['%s: theta final = %+.2f deg, max|theta| overall %.2f deg / after 2 s %.2f deg; ' ...
        'biceps F = %.1f N, A_bi %.3f -> %.3f; triceps A_tri %.3f -> %.3f (m_cup %.2f kg)\n'], ...
        runs(r).label, runs(r).th.Data(end)*180/pi, ...
        max(abs(runs(r).th.Data))*180/pi, max(abs(runs(r).th.Data(w)))*180/pi, ...
        runs(r).F.Data(end), runs(r).A.Data(1), runs(r).A.Data(end), ...
        runs(r).At.Data(1), runs(r).At.Data(end), runs(r).mCup.Data(end));
end
assignin('base', 'kReflex', 1);   % restore

thOn = runs(1).th; thOff = runs(2).th;
fig = figure('Visible', 'off', 'Position', [100 100 900 1150]);
subplot(5,1,1); plot(thOff.Time, thOff.Data*180/pi, 'LineWidth', 1.4); hold on;
plot(thOn.Time, thOn.Data*180/pi, 'LineWidth', 1.6); grid on;
ylabel('cup sag \theta (deg)'); legend('reflex OFF', 'reflex ON', 'Location', 'southeast');
title('Beer-cup demo: Ia/Ib reflexes + antagonist inhibition hold the cup level while beer pours in (starts in equilibrium)');
subplot(5,1,2); plot(runs(1).mCup.Time, runs(1).mCup.Data, 'LineWidth', 1.4); grid on;
ylabel('beer mass (kg)');
subplot(5,1,3); plot(runs(1).A.Time, runs(1).A.Data, 'LineWidth', 1.6); hold on;
plot(runs(1).At.Time, runs(1).At.Data, 'LineWidth', 1.6, 'LineStyle', '--'); ...
plot(runs(2).A.Time, runs(2).A.Data, 'LineWidth', 1.0, 'Color', [0.6 0.6 0.6]); grid on;
ylabel('activation'); legend('biceps (ON)', 'triceps (ON)', 'biceps (OFF)');
title('biceps activation RISES with the pour while triceps activation FALLS (reciprocal Ia inhibition)');
subplot(5,1,4); plot(runs(1).F.Time, runs(1).F.Data, 'LineWidth', 1.4); hold on;
plot(runs(1).Ft.Time, runs(1).Ft.Data, 'LineWidth', 1.4, 'LineStyle', '--'); grid on;
ylabel('BPA force (N)'); legend('biceps (ON)', 'triceps (ON)');
subplot(5,1,5); plot(runs(2).F.Time, runs(2).F.Data, 'LineWidth', 1.4, 'LineStyle', '--'); grid on;
ylabel('biceps force, OFF (N)'); xlabel('time (s)');
exportgraphics(fig, fullfile('..','results','pictures','sns_beer_results.png'), 'Resolution', 150);
save(fullfile('..','results','sns_beer_results.mat'), 'runs');
sagOn  = max(abs(thOn.Data(thOn.Time > 2)))*180/pi;
sagOff = max(abs(thOff.Data(thOff.Time > 2)))*180/pi;
fprintf('max sag after 2 s: reflex ON %.2f deg, reflex OFF %.2f deg\n', sagOn, sagOff);
fprintf('triceps activation (ON): %.3f -> %.3f over the pour\n', runs(1).At.Data(1), runs(1).At.Data(end));
fprintf('BEER DEMO RUN OK\n');
