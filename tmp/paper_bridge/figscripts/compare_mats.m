%compare_mats.m — what changed between the archived ext-10mm-52cm mat and
%the fresh rerun? Constants, routing, Xi, and curve finiteness.
root = 'D:/GitHub/Bipedal_Robot';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
festo = fullfile(root, 'Testing_Data/2022_02_Festo');
A = load(fullfile(festo, 'KneeExt_10mm_52cm_stale_20261019.mat'));
B = load(fullfile(festo, 'KneeExt_10mm_52cm.mat'));
for f = {'rest','kmax','tendon','fitting','pres','Dia','CrossPoint'}
    fa = A.(f{1}); fb = B.(f{1});
    if isnumeric(fa), eq = isequal(fa, fb); else, eq = strcmp(fa, fb); end
    fprintf('CONST %-10s old=%s  new=%s  same=%d\n', f{1}, mat2str(fa, 6), mat2str(fb, 6), eq);
end
fprintf('LOCATION same: %d (old size %s, new %s)\n', isequal(A.Location, B.Location), ...
    mat2str(size(A.Location)), mat2str(size(B.Location)));
fprintf('OLD TorqueR finite: %d/%d\n', sum(isfinite(A.TorqueR(:,3))), numel(A.TorqueR(:,3)));
if isprop(A.Vas_Pam_adj, 'Torque_p')
    fprintf('OLD adj Torque_p finite: %d/%d\n', sum(isfinite(A.Vas_Pam_adj.Torque_p(:,3))), numel(A.Vas_Pam_adj.Torque_p(:,3)));
    fprintf('OLD adj Xi: %.6g %.6g %.6g %.6g\n', A.Vas_Pam_adj.Xi0, A.Vas_Pam_adj.Xi1, A.Vas_Pam_adj.Xi2, A.Vas_Pam_adj.Xi3);
else
    fn = fieldnames(A.Vas_Pam_adj);
    fprintf('OLD adj fields: %s\n', strjoin(fn', ', '));
end
fprintf('OLD Vas_Pam.MuscleLength(1:5): %s\n', mat2str(A.Vas_Pam.MuscleLength(1:5)', 5));
fprintf('NEW Vas_Pam.MuscleLength(1:5): %s\n', mat2str(B.Vas_Pam.MuscleLength(1:5)', 5));
fprintf('OLD TorqueZ: %s\n', mat2str(A.TorqueZ', 4));
fprintf('OLD contraction finite: %d/%d ; max=%.4f min=%.4f\n', ...
    sum(isfinite(A.Vas_Pam.Contraction)), numel(A.Vas_Pam.Contraction), ...
    max(A.Vas_Pam.Contraction([], 'omitnan'), [], 'omitnan'), min(A.Vas_Pam.Contraction, [], 'omitnan'));
fprintf('NEW contraction finite: %d/%d ; max=%.4f min=%.4f\n', ...
    sum(isfinite(B.Vas_Pam.Contraction)), numel(B.Vas_Pam.Contraction), ...
    max(B.Vas_Pam.Contraction, [], 'omitnan'), min(B.Vas_Pam.Contraction, [], 'omitnan'));
disp('DONE_compare');
