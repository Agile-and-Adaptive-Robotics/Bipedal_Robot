%Xi3 landscape under the NEW extensor bracket point (rib midpoint).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
POOL = [1 2 5 6 7 8];
ALL9 = 1:9;
%Baseline a0: load the stored one (saves an Inf-stiffness evaluation that
%spams "Matrix is singular" warnings from minimizeExtX3's rigid case).
a0file = 'minimizeExt10mmX3_results_20260907_2trans.mat';
if exist(a0file, 'file')
    a0 = load(a0file, 'a0'); a0 = a0.a0;
else
    warnstate = warning('off', 'MATLAB:singularMatrix');
    a0 = minimizeExtX3(0, Inf, Inf, 0);
    warning(warnstate);
end
g = [1.508e4, 1.218e4];                             %flexor Xi1/Xi2 (locked in the CVs)
xi3s = [0, 0.05, 0.10, 0.159, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00];
xi0s = [-0.002, -0.0064, -0.012];
fprintf('New-Pbr pool(6) mean RMSE, rows=Xi0, cols=Xi3   (baseline pool mean = %.3f)\n', mean(a0(POOL,1)));
fprintf('%9s', 'Xi0\Xi3'); fprintf('%7.3f', xi3s); fprintf('\n');
best = struct('m', Inf, 'xi0', 0, 'xi3', 0);
for i0 = 1:numel(xi0s)
    fprintf('%9.4f', xi0s(i0));
    for i3 = 1:numel(xi3s)
        f = minimizeExtX3(xi0s(i0), g(1), g(2), xi3s(i3));
        m = mean(f(POOL,1));
        fprintf('%7.3f', m);
        if m < best.m, best.m = m; best.xi0 = xi0s(i0); best.xi3 = xi3s(i3); end
    end
    fprintf('\n');
end
fprintf('Best: pool RMSE %.3f at Xi0=%.4f, Xi3=%.3f\n\n', best.m, best.xi0, best.xi3);

fprintf('New-Pbr pool(6) mean FVU normalized by baseline, rows=Xi0, cols=Xi3 (1.00 = baseline)\n');
fprintf('%9s', 'Xi0\Xi3'); fprintf('%7.3f', xi3s); fprintf('\n');
for i0 = 1:numel(xi0s)
    fprintf('%9.4f', xi0s(i0));
    for i3 = 1:numel(xi3s)
        f = minimizeExtX3(xi0s(i0), g(1), g(2), xi3s(i3));
        m = mean(f(POOL,2) ./ a0(POOL,2));
        fprintf('%7.3f', m);
    end
    fprintf('\n');
end
fprintf('\n');
fprintf('%9s', 'Xi0\Xi3'); fprintf('%7.3f', xi3s); fprintf('\n');
for i0 = 1:numel(xi0s)
    fprintf('%9.4f', xi0s(i0));
    for i3 = 1:numel(xi3s)
        f = minimizeExtX3(xi0s(i0), g(1), g(2), xi3s(i3));
        m = mean(f(ALL9,1));
        fprintf('%7.3f', m);
    end
    fprintf('\n');
end
fprintf('(second table = mean over all 9 tests)\n');
