function goal2_dbg_rgmn()
% Debug SNS_SpikingRG_MN: why do the RG half-centers never spike?
rep = fileparts(mfilename('fullpath'));
d = rep;
while exist(fullfile(d, 'Code', 'Matlab', 'SNS_Simscape', 'SNS_Library.slx'), 'file') == 0
    pn = fileparts(d);
    if strcmp(pn, d), error('repo root not found'); end
    d = pn;
end
cd(fullfile(d, 'Code', 'Matlab', 'SNS_Simscape', 'demos'));
addpath(fileparts(pwd));   % SNS_Library lives in the parent
mdl = 'SNS_SpikingRG_MN';
load_system('SNS_Library');
load_system(mdl);
set_param(mdl, 'StopTime', '0.5');
out = sim(mdl);
vE = out.log_V_RG_E; sE = out.log_spk_RG_E; m1 = out.log_V_MN_E1;
fprintf('V_RG_E: min %.2f max %.2f (n=%d)\n', min(vE.Data), max(vE.Data), numel(vE.Data));
fprintf('spk_RG_E: max %.1f, count(>=0.5) %d\n', max(sE.Data), sum(sE.Data >= 0.5));
fprintf('V_MN_E1: min %.2f max %.2f\n', min(m1.Data), max(m1.Data));
w = vE.Time <= 0.1;
fprintf('first 100 ms of V_RG_E (every ~10th sample):\n');
idx = 1:10:numel(vE.Data);
fprintf('%.3f ', vE.Data(idx(w)));
fprintf('\n');
fprintf('times: '); fprintf('%.3f ', vE.Time(idx(w))); fprintf('\n');
end
