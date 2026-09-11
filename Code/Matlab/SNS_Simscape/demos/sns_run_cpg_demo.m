%% sns_run_cpg_demo.m — run BPACPGLegDemo, verify alternating CPG, save plots
cdto = fileparts(mfilename('fullpath'));
cd(cdto);
addpath(cdto);
addpath(fileparts(cdto));   % SNS_Library lives in the parent folder

load_system('SNS_Library');
load_system('BPACPGLegDemo');
set_param('BPACPGLegDemo', 'StopTime', '10');
out = sim('BPACPGLegDemo');

th  = out.log_th;   Ae = out.log_A_ext; Af = out.log_A_flex;
Vex = out.log_V_RG_ext; Vfl = out.log_V_RG_flex;
Fe  = out.log_F_ext;    Ff  = out.log_F_flex;
t = th.Time;

% oscillation check: sign changes of the RG voltage difference
dV = Vex.Data - Vfl.Data;
s = sign(dV); s(s == 0) = 1;
switches = sum(abs(diff(s)) > 0);
period = 0;
if switches >= 2
    % mean half-period from crossing times
    idx = find(abs(diff(s)) > 0);
    tc = t(idx);
    period = 2*mean(diff(tc));
end
fprintf('CPG: %d half-cycles in %g s', switches, t(end));
if period > 0
    fprintf(' (period ~ %.3f s, freq ~ %.2f Hz)\n', period, 1/period);
else
    fprintf('\n');
end

fig = figure('Visible', 'off', 'Position', [100 100 900 1000]);
subplot(5,1,1); plot(t, th.Data*180/pi, 'LineWidth', 1.4); grid on;
ylabel('knee \theta (deg)'); title('SNS CPG leg demo — half-center RG \rightarrow antagonist BPAs on 1-DOF knee');
subplot(5,1,2); plot(t, Vex.Data, 'LineWidth', 1.4); hold on; plot(t, Vfl.Data, 'LineWidth', 1.4); grid on;
ylabel('RG membrane (mV)'); legend('RG_{ext}', 'RG_{flex}');
subplot(5,1,3); plot(t, Ae.Data, 'LineWidth', 1.4); hold on; plot(t, Af.Data, 'LineWidth', 1.4); grid on;
ylabel('activation'); legend('extensor', 'flexor'); ylim([-0.05 1.05]);
subplot(5,1,4); plot(t, Fe.Data, 'LineWidth', 1.4); hold on; plot(t, Ff.Data, 'LineWidth', 1.4); grid on;
ylabel('BPA force (N)'); legend('extensor', 'flexor');
subplot(5,1,5); plot(t, th.Data*180/pi, 'LineWidth', 1.4); grid on;
ylabel('\theta zoom (deg)'); xlabel('time (s)');
exportgraphics(fig, fullfile('..','results','pictures','sns_cpg_results.png'), 'Resolution', 150);
save(fullfile('..','results','sns_cpg_results.mat'), 'out');
fprintf('CPG run OK: theta range %.1f..%.1f deg\n', min(th.Data)*180/pi, max(th.Data)*180/pi);
