%% sns_run_beer_demo.m — run BeerCupReflexDemo twice (reflex ON vs OFF), save plots
%
% The story: beer pours into the cup (0 -> 0.5 kg over 8 s). With the Ia/Ib
% reflex pathway enabled the biceps reflexively adds force and the cup stays
% near level; with kReflex = 0 the arm sags under the growing load.
cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder

load_system('SNS_Library');
load_system('BeerCupReflexDemo');

runs = struct('label', {'reflex ON', 'reflex OFF'}, 'k', {1, 0});
for r = 1:2
    set_param('BeerCupReflexDemo/syn_Ia_exc', 'gmax', sprintf('%.4g*%d', 0.005, runs(r).k));
    set_param('BeerCupReflexDemo/syn_Ib_inh', 'gmax', sprintf('%.4g*%d', 0.0015, runs(r).k));
    out = sim('BeerCupReflexDemo');
    runs(r).th = out.log_th;  runs(r).mCup = out.log_mCup;
    runs(r).F  = out.log_F_bi; runs(r).A = out.log_A_bi;
    fprintf('%s: theta final = %+.2f deg, biceps F = %.1f N (m_cup %.2f kg)\n', ...
        runs(r).label, runs(r).th.Data(end)*180/pi, runs(r).F.Data(end), runs(r).mCup.Data(end));
end
set_param('BeerCupReflexDemo/syn_Ia_exc', 'gmax', sprintf('%.4g*%d', 0.005, 1));   % restore
set_param('BeerCupReflexDemo/syn_Ib_inh', 'gmax', sprintf('%.4g*%d', 0.0015, 1));

thOn = runs(1).th; thOff = runs(2).th;
fig = figure('Visible', 'off', 'Position', [100 100 900 950]);
subplot(4,1,1); plot(thOff.Time, thOff.Data*180/pi, 'LineWidth', 1.4); hold on;
plot(thOn.Time, thOn.Data*180/pi, 'LineWidth', 1.6); grid on;
ylabel('cup sag \theta (deg)'); legend('reflex OFF', 'reflex ON', 'Location', 'southeast');
title('Beer-cup demo: Ia stretch reflex + Ib autogenic inhibition hold the cup level while beer pours in');
subplot(4,1,2); plot(runs(1).mCup.Time, runs(1).mCup.Data, 'LineWidth', 1.4); grid on;
ylabel('beer mass (kg)');
subplot(4,1,3); plot(runs(1).F.Time, runs(1).F.Data, 'LineWidth', 1.4); hold on;
plot(runs(2).F.Time, runs(2).F.Data, 'LineWidth', 1.4, 'LineStyle', '--'); grid on;
ylabel('biceps BPA force (N)'); legend('reflex ON', 'reflex OFF');
subplot(4,1,4); plot(runs(1).A.Time, runs(1).A.Data, 'LineWidth', 1.4); hold on;
plot(runs(2).A.Time, runs(2).A.Data, 'LineWidth', 1.4, 'LineStyle', '--'); grid on;
ylabel('activation'); xlabel('time (s)'); legend('reflex ON', 'reflex OFF');
exportgraphics(fig, fullfile('..','results','pictures','sns_beer_results.png'), 'Resolution', 150);
save(fullfile('..','results','sns_beer_results.mat'), 'runs');
sagOn  = max(abs(thOn.Data))*180/pi;
sagOff = max(abs(thOff.Data))*180/pi;
fprintf('max sag: reflex ON %.2f deg, reflex OFF %.2f deg\n', sagOn, sagOff);
fprintf('BEER DEMO RUN OK\n');
