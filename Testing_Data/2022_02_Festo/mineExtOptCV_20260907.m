%Mine the optimized extensor CV (20260819): which 3-holdout combos generalize?
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
S = load('minimizeExtPin10_results_20260819_2transforms_Z2.mat');
cv = S.results_cv;
if isfield(S,'validLabels'), lab = string(S.validLabels); end
lab = ["40cm","40cm-tendon","42cm","42cm-tendon","43cm","43cm-tendon","46cm","47cm","48cm"]; %foldIdx uses full 9-test numbering
a0 = minimizeExtX3(0, Inf, Inf, 0);
fprintf('===== EXTENSOR optimized CV (10 folds, numHold=3, Xi1/Xi2 locked to flexor) =====\n');
fprintf('%-28s %8s %8s %8s %8s %9s %7s\n','holdout','Xi0(m)','Xi3','trRMSE','valRMSE','dist','passPool');
rows = cell(numel(cv),1);
for k = 1:numel(cv)
    d = cv{k}.distance_all(:); [~,b] = min(d);
    xa = cv{k}.optParams_all(b,:);
    Xi0 = xa(1)/100; Xi3 = xa(4);
    tr = cv{k}.trainScores_all(b,:); va = cv{k}.validation_all(b,:);
    ho = cv{k}.foldIdx(b,:);
    F = minimizeExtX3(Xi0, 10^xa(2), 10^xa(3), Xi3);
    pass = all(F(ho,1:2) <= a0(ho,1:2), 'all');
    hn = strjoin(lab(ho)','+');
    fprintf('%-28s %8.4f %8.3f %8.3f %8.3f %9.3f %7d\n', hn, Xi0, Xi3, tr(1), va(1), d(b), pass);
    rows{k} = struct('holdout',hn,'dist',d(b),'Xi3',Xi3,'pass',pass);
end
dists = cellfun(@(r) r.dist, rows);
[~,ord] = sort(dists);
fprintf('\nExtensor holdout combos ranked by generalization (best first):\n');
for i = 1:numel(ord)
    fprintf('  %2d. %-28s dist=%.3f  Xi3=%.2f  pass=%d\n', i, rows{ord(i)}.holdout, rows{ord(i)}.dist, rows{ord(i)}.Xi3, rows{ord(i)}.pass);
end
x3 = cellfun(@(r) r.Xi3, rows);
fprintf('\nXi3 at best candidate per fold: min %.2f  median %.2f  max %.2f\n', min(x3), median(x3), max(x3));
