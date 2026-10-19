%dump_FlxPinDiag.m — score every candidate "Unoptimized" series + Hunt on the
%evaluator's own sampling basis (griddedInterpolant linear onto Aexp) to find
%which curve the figures should carry and why the preopt table disagreed.
root = 'D:/GitHub/Bipedal_Robot';
festo = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festo); warning('off','all'); setenv('FLX_T5_YMM','');

[f0, b0] = minimizeFlxPin(0, Inf, Inf, 1:5);      % zero-correction
g = [3.937370e-3, 3.997873e+04, 1.473375e+04];    % row 77 ledger
[f1, b1] = minimizeFlxPin(g(1), g(2), g(3), 1:5); % adopted
load('ForceStrainForFit.mat','z'); z(isnan(z)) = 0;
Xp = linspace(0,620,20); Xp = Xp(2:20); Yp = linspace(0,1,30);

gofat = @(Ak, Mv, Aexp, Mexp) local_gof(Ak, Mv, Aexp, Mexp);  %#ok<NASGU>

for i = [1 2 4 5]
    % Hunt series on zero-correction kinematics (identical call to original)
    Yq = b0(i).strain ./ ((b0(i).rest - b0(i).Kmax) / b0(i).rest);
    Xq = b0(i).P .* ones(size(Yq));
    Vq = interp2(Xp, Yp', z, Xq, Yq, 'spline');
    Fold = Vq .* b0(i).unitD;
    Fq = hypot(Fold(:,1), Fold(:,2));
    Mo = -b0(i).mA .* Fq;
    Mo(Yq > 1) = 0;
    Mo(b0(i).strain < 0) = NaN;

    cand = {'Hunt(Mold)', 'stale .M', 'M_p zero-corr', 'M_p adopted'};
    ser = {Mo, b0(i).M, b0(i).M_p(:,3), b1(i).M_p(:,3)};
    fprintf('--- test %d (f0 = %.4f/%.4f/%.4f, f1 = %.4f/%.4f/%.4f) ---\n', ...
        i, f0(i,:), f1(i,:));
    for c = 1:numel(ser)
        Mv = ser{c}(:);
        [tRMSE, tFVU, tMax] = local_gof(b0(i).Ak, Mv, b0(i).Aexp, b0(i).Mexp);
        fprintf('  %-14s RMSE %7.4f  FVU %7.4f  Max %7.4f\n', cand{c}, tRMSE, tFVU, tMax);
    end
    % hybrid: paired points (A_h, M_h) vs measured
    if isequal(b0(i).A_h(:), b0(i).Aexp(:))
        rh = b0(i).M_h(:) - b0(i).Mexp(:);
    else
        rh = interp1(b0(i).A_h(:), b0(i).M_h(:), b0(i).Aexp(:), 'linear', 'extrap') - b0(i).Mexp(:);
    end
    vh = isfinite(rh);
    fprintf('  %-14s RMSE %7.4f  FVU %7.4f  Max %7.4f  (paired A_h==Aexp: %d)\n', ...
        'Hybrid M_h', sqrt(mean(rh(vh).^2)), ...
        sum(rh(vh).^2)/sum((b0(i).Mexp(vh)-mean(b0(i).Mexp(vh))).^2), ...
        max(abs(rh(vh))), isequal(b0(i).A_h(:), b0(i).Aexp(:)));
    fprintf('  |stale.M - M_p0| max %.4f   |stale.M - Mold| max %.4f   numel(A_h vs Aexp): %d/%d\n', ...
        max(abs(b0(i).M - b0(i).M_p(:,3))), max(abs(b0(i).M - Mo)), ...
        numel(b0(i).A_h), numel(b0(i).Aexp));
end
disp('DONE_diag');

function [RMSE, fvu, mx] = local_gof(Ak, Mv, Aexp, Mexp)
    [As, ix] = sort(Ak(:));
    G = griddedInterpolant(As, Mv(ix));
    Mp = G(Aexp(:));
    r = Mp - Mexp(:);
    v = ~isnan(r) & isfinite(Mexp(:));
    RMSE = sqrt(mean(r(v).^2));
    fvu = sum(r(v).^2) / sum((Mexp(v) - mean(Mexp(v))).^2);
    mx = max(abs(r(v)));
end
