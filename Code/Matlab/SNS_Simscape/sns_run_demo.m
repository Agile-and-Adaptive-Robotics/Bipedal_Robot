%% sns_run_demo.m — run KneeReflexDemo and save diagnostic plots
% defines parameters, simulates, saves PNGs + a .mat of logged signals.
cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);

% model params (same as PreLoadFcn)
I_knee = 0.06; b_knee = 0.5; K_knee = 0.5; th0 = 0.26; ROM = 1.5708;
r_arm = 0.035; Tload = 0.5; desc_ext = 4.2; desc_flex = 2.5;
Fmax_ext = 500; Fmax_flex = 450; epsScale = 0.15;

load_system('SNS_Library');
load_system('KneeReflexDemo');
set_param('KneeReflexDemo', 'StopTime', '5');
out = sim('KneeReflexDemo');

th   = out.log_th;    thd  = out.log_thd;
Ae   = out.log_A_ext; Af   = out.log_A_flex;
Vme  = out.log_V_MN_ext; Vmf = out.log_V_MN_flex;
Fe   = out.log_F_ext; Ff  = out.log_F_flex;

fig = figure('Visible', 'off', 'Position', [100 100 900 1100]);
subplot(5,1,1); plot(th.Time, th.Data*180/pi, 'LineWidth', 1.2); grid on;
ylabel('knee \theta (deg)'); title('SNS knee reflex demo — non-spiking RC neurons + E/I synapses + BPA antagonist pair');
subplot(5,1,2); plot(Vme.Time, Vme.Data, 'LineWidth', 1.2); hold on; plot(Vmf.Time, Vmf.Data, 'LineWidth', 1.2); grid on;
ylabel('MN membrane (mV)'); legend('MN_{ext}', 'MN_{flex}'); ylim([-60 -30]);
subplot(5,1,3); plot(Ae.Time, Ae.Data, 'LineWidth', 1.2); hold on; plot(Af.Time, Af.Data, 'LineWidth', 1.2); grid on;
ylabel('activation'); legend('extensor', 'flexor'); ylim([0 1.05]);
subplot(5,1,4); plot(Fe.Time, Fe.Data, 'LineWidth', 1.2); hold on; plot(Ff.Time, Ff.Data, 'LineWidth', 1.2); grid on;
ylabel('BPA force (N)'); legend('extensor', 'flexor');
subplot(5,1,5); plot(thd.Time, thd.Data, 'LineWidth', 1.2); grid on;
ylabel('vel (rad/s)'); xlabel('time (s)');
exportgraphics(fig, 'sns_demo_results.png', 'Resolution', 150);
save('sns_demo_results.mat', 'out');
fprintf('sim OK: theta final = %.2f deg, A_ext = %.2f, A_flex = %.2f\n', th.Data(end)*180/pi, Ae.Data(end), Af.Data(end));
