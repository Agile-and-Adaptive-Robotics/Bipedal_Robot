% Probe Vas_Pam result mats for xBest/XiUsed and compare against builder Xi
R = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(R,'Code','Matlab')));
addpath(fullfile(R,'Code','Matlab','Mesh_Optimization'));
addpath(fullfile(R,'Testing_Data','2022_02_Festo'), '-end');

ctx = buildKneeExtContext20mm();
xiBuilt = [ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.Xi3];
fprintf('builder Xi: [%.6g %.6g %.6g %.6g]\n', xiBuilt);

d = dir(fullfile(R,'Code/Matlab/Mesh_Optimization/Results','Vas_Pam_20mm_Result*.mat'));
for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name);
    info = whos('-file', f);
    vars = strjoin({info.name}, ',');
    hasX = any(strcmp({info.name},'xBest'));
    hasXi = any(strcmp({info.name},'XiUsed'));
    fprintf('\n%s\n  vars: %s\n', d(k).name, vars);
    if hasX
        S = load(f,'xBest');
        fprintf('  xBest = [%.4f %.4f %.4f | %.4f %.4f %.4f | rest %.4f tendon %.4f]\n', ...
            S.xBest(1),S.xBest(2),S.xBest(3),S.xBest(4),S.xBest(5),S.xBest(6),S.xBest(7),S.xBest(8));
    end
    if hasXi
        S2 = load(f,'XiUsed');
        fprintf('  XiUsed = [%.6g %.6g %.6g %.6g]  matchBuilder=%d\n', S2.XiUsed, isequal(S2.XiUsed, xiBuilt));
    end
end
