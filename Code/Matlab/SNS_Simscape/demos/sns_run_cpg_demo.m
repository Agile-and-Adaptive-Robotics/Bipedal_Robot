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

% oscillation check: HYSTERESIS state switches of the RG voltage difference
% (plain zero-crossing counts boundary chatter as extra half-cycles)
dV = Vex.Data - Vfl.Data;
hyst = 0.1*max(dV) - 0.1*min(dV);          % 10% of the dV swing
upLvl = max(dV) - hyst;  dnLvl = min(dV) + hyst;
state = 1*(dV(1) >= 0);  tsw = [];
for i = 2:numel(dV)
    if state == 0 && dV(i) > upLvl, state = 1; tsw(end+1) = t(i); %#ok<SAGROW>
    elseif state == 1 && dV(i) < dnLvl, state = 0; tsw(end+1) = t(i); %#ok<SAGROW>
    end
end
switches = numel(tsw);
period = 2*mean(diff(tsw));
fprintf('CPG: %d hysteresis switches in %g s', switches, t(end));
if switches >= 2
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
