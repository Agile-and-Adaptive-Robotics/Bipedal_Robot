% Dump per-test resting + tendon lengths for the GoF table columns
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
fprintf('--- FLEXOR pinned tests ---\n');
fx = {'Plot_KneeFlxPin10mm_48cm.mat', 'Plot_KneeFlxPin10mm_46cm.mat', ...
      'Plot_KneeFlxPin10mm_47cm.mat', 'Plot_KneeFlxPin10mm_40cm.mat', ...
      'Plot_KneeFlxPin10mm_42cm.mat'};
for i = 1:numel(fx)
    S = load(fx{i}, 'Bifemsh_Pam');
    fprintf('%s: rest = %.4f m, tendon = %.4f m (%.1f mm)\n', fx{i}, ...
        S.Bifemsh_Pam.RestingL, S.Bifemsh_Pam.TendonL, ...
        1000 * S.Bifemsh_Pam.TendonL);
end
fprintf('--- EXTENSOR pinned tests (ExtPinBPASet ke struct array) ---\n');
E = load('ExtPinBPASet.mat', 'ke');
for i = 1:numel(E.ke)
    fprintf('ke(%d): rest = %.4f m, tendon = %.4f m (%.1f mm)\n', i, ...
        E.ke(i).rest, E.ke(i).ten, 1000 * E.ke(i).ten);
end
