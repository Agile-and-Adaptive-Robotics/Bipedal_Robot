function goal4_export_knee_traces()
% GOAL 4 (dissertation figures) — export knee-angle + MN-voltage traces of
% the committed KneeReflexDemo (non-spiking baseline) and its spiking twin
% KneeReflexDemo_Spiking so the dissertation figure can be drawn to the
% project figure standards in Python. Read-only with respect to the models:
% load_system + sim, StopTime only (same as goal2_run_baselines.m).
%
% Output: goal4_knee_traces.mat next to this file with
%   base.t/th  baseline theta (rad), base.Vme/Vmf  MN voltages (mV)
%   spk.*      same channels for the spiking twin
% Metrics printed for cross-check against goal2_knee_spiking.mat:
%   baseline rise 15.79 deg, settle 43.5 (41.6-44.9); twin rise 15.79,
%   settle 44.4 (41.9-46.7).

rep = fileparts(mfilename('fullpath'));
d = rep;
while exist(fullfile(d, 'Code', 'Matlab', 'SNS_Simscape', 'SNS_Library.slx'), 'file') == 0
    dn = fileparts(d);
    if strcmp(dn, d), error('repo root not found'); end
    d = dn;
end
sns = fullfile(d, 'Code', 'Matlab', 'SNS_Simscape');
addpath(sns); addpath(fullfile(sns, 'demos'));
cd(fullfile(sns, 'demos'));
load_system('SNS_Library');

base = grab('KneeReflexDemo');
spk  = grab('KneeReflexDemo_Spiking');

save(fullfile(rep, 'goal4_knee_traces.mat'), 'base', 'spk');
fprintf('saved goal4_knee_traces.mat\n');

for nm = {'base', 'spk'}
    S = eval(nm{1});
    w = S.t > S.t(end) - 2;
    fprintf('%-4s rise %.2f deg | settle mean %.1f (%.1f-%.1f) | final %.2f\n', ...
        nm{1}, interp1(S.t, S.th, 0.05, 'linear')*180/pi, ...
        mean(S.th(w))*180/pi, min(S.th(w))*180/pi, max(S.th(w))*180/pi, ...
        S.th(end)*180/pi);
end

    function S = grab(model)
        load_system(model);
        set_param(model, 'StopTime', '5');
        out = sim(model);
        S.t   = out.log_th.Time;
        S.th  = out.log_th.Data(:)';            % rad
        try
            S.Vme = out.log_V_MN_ext.Data(:)';  % mV
            S.Vmf = out.log_V_MN_flex.Data(:)'; % mV
        catch
            S.Vme = []; S.Vmf = [];             % baseline demo may not log MN V
        end
    end
end
