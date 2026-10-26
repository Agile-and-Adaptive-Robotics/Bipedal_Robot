% Air-stepping smoke run: simulate SNS_AirStep 40 s and verify curriculum
% stages 1-2 behavior (both legs stepping, antiphase across legs).
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root); addpath(fullfile(root, 'demos'));
mdl = 'SNS_AirStep';
load_system(fullfile(root, 'demos', [mdl '.slx']));
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'sigs', ...
          'StopTime', '40');
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
close_system(mdl, 0);
sigs = out.sigs;
jn = {'theta_L_hip', 'theta_L_knee', 'theta_R_hip', 'theta_R_knee'};
vals = cell(1, 4); tms = cell(1, 4);
for k = 1:4
    ts = sigs.get(jn{k}).Values;
    vals{k} = ts.Data(:); tms{k} = ts.Time(:);
end
names = {'L hip', 'L knee', 'R hip', 'R knee'};
for k = 1:4
    v = rad2deg(vals{k}); t = tms{k};
    ipk = find(islocalmax(v, 'MinProminence', 5));
    ipk = ipk(t(ipk) > 2);
    if numel(ipk) >= 2, per = sprintf('%.2f', mean(diff(t(ipk)))); else, per = 'n/a'; end
    fprintf('%-8s: range [%7.2f, %7.2f] deg | %3d cycles | period %s s\n', ...
        names{k}, min(v), max(v), numel(ipk), per);
end
% antiphase check: L hip vs R hip correlation
n = min(numel(vals{1}), numel(vals{3}));
a = vals{1}(1:n) - mean(vals{1}(1:n));
b_ = vals{3}(1:n) - mean(vals{3}(1:n));
r = sum(a .* b_) / sqrt(sum(a .^ 2) * sum(b_ .^ 2));
fprintf('L-vs-R hip correlation: %.3f (negative = antiphase)\n', r);
ok = r < -0.3;
if ok
    disp('AIR STEP SMOKE: PASS - bilateral antiphase stepping under the shared CPG');
else
    disp('AIR STEP SMOKE: CHECK - inspect traces');
end
