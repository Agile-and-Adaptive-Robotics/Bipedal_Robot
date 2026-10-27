%dbg_Hunt_20261026.m -- instrument the Hunt audit's domain/GoF mechanics
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');

Z = load(fullfile(root, 'Code/Matlab/Functions/ForceStrainForFit.mat'));
z = Z.z;
fprintf('z size [%d %d], vars: %s\n', size(z,1), size(z,2), strjoin(fieldnames(Z).', ', '));
fprintf('z NaN per col: '); fprintf('%3d', sum(isnan(z),1)); fprintf('\n');
fprintf('z(1,:) = '); fprintf('%9.3f', z(1,:)); fprintf('\n');
fprintf('z(30,:) = '); fprintf('%9.3f', z(30,:)); fprintf('\n');

[f0, bpa0] = minimizeFlxPin(0, Inf, Inf, 1:5);
ordr = [1 5 2 4];
names = {'48.5','41.0','45.7','40.6t'};
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    KMAX = (b.rest - b.Kmax)/b.rest;
    fprintf('\ntest %s: numel Ak=%d Aexp=%d A_h=%d; M=%s Mexp=%s M_h=%s\n', ...
        names{k}, numel(b.Ak), numel(b.Aexp), numel(b.A_h), ...
        mat2str(size(b.M)), mat2str(size(b.Mexp)), mat2str(size(b.M_h)));
    fprintf('  Ak range [%.2f %.2f]; Aexp range [%.2f %.2f]; A_h range [%.2f %.2f]\n', ...
        min(b.Ak), max(b.Ak), min(b.Aexp), max(b.Aexp), min(b.A_h), max(b.A_h));
    fprintf('  rest=%.4f Kmax=%.4f KMAX=%.4f P uniq=%s strain range [%.4f %.4f]\n', ...
        b.rest, b.Kmax, KMAX, mat2str(unique(round(b.P,2))), min(b.strain), max(b.strain));
    fprintf('  strain>=0 count over Ak: %d of %d\n', sum(b.strain >= 0), numel(b.strain));
    fprintf('  Ak where strain crosses 0: ');
    idx = find(diff(sign(b.strain)) ~= 0);
    fprintf('%s\n', mat2str(b.Ak(idx).'));
end

% reproduce lookup series for test 1 and print domain + residuals
i = ordr(1); b = bpa0(1);
zz = z; zz(isnan(zz)) = 0;
Xp = linspace(0, 620, 20); Xp = Xp(2:20);
Yp = linspace(0, 1, 30);
Yq = b.strain ./ ((b.rest - b.Kmax) / b.rest);
Xq = b.P .* ones(size(Yq));
Vq = interp2(Xp, Yp', zz, Xq, Yq, 'spline');
Fold = Vq .* b.unitD;
Fq = hypot(Fold(:,1), Fold(:,2));
Mo = -b.mA .* Fq;
Mo(Yq > 1) = 0;
Mo(b.strain < 0) = NaN;
[Ak_s, ia] = sort(b.Ak);
sv = Mo(ia); valid = isfinite(sv); svv = sv(valid);
fprintf('\ntest48.5 lookup: valid grid pts %d of %d; valid Ak range [%.2f %.2f]\n', ...
    sum(valid), numel(valid), min(Ak_s(valid)), max(Ak_s(valid)));
G = griddedInterpolant(Ak_s(valid), svv);
sel = b.Aexp >= min(Ak_s(valid)) & b.Aexp <= max(Ak_s(valid));
fprintf('sel count: %d of %d; sel Aexp: %s\n', sum(sel), numel(sel), mat2str(b.Aexp(sel).'));
yv = b.Mexp(sel); pv = G(b.Aexp(sel));
fprintf('yv size %s; pv size %s\n', mat2str(size(yv)), mat2str(size(pv)));
fprintf('yv: '); fprintf('%8.3f', yv); fprintf('\n');
fprintf('pv: '); fprintf('%8.3f', pv); fprintf('\n');
[r, fv, mx] = Go_OfF(yv, pv);
fprintf('Go_OfF: %.4f %.4f %.4f\n', r, fv, mx);
yresid = yv - pv;
fprintf('resid: '); fprintf('%8.3f', yresid); fprintf('\n');
n = numel(yv);
SSresid = sum(yresid.^2,'omitnan');
SStot = sum((yv-mean(yv,1)).^2);
fprintf('manual: n=%d SSresid=%.4f SStot=%.4f RMSE=%.4f FVU=%.4f Max=%.4f\n', ...
    n, SSresid, SStot, sqrt(SSresid/n), SSresid*(n-1)/(SStot*n), max(abs(yresid)));
disp('DONE_dbg');
