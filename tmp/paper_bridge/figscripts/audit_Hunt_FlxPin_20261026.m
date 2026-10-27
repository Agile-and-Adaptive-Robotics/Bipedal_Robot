%audit_Hunt_FlxPin_20261026.m -- P1-b numerical audit of the Hunt-equation
% series in fig:FlxPin / tab:FlxPinGof / tab:FlxPinPreGof.
% READ-ONLY over Ben's sources (HuntEq.m, ForceStrainForFit.mat,
% minimizeFlxPin.m). Questions, answered numerically:
%   Q1 what does ForceStrainForFit.mat's z actually contain (units, grid,
%      NaN pattern), and does it equal a fresh HuntEq table (isokinetic S=0,
%      epsilon grid of LookupTable_generate.m)?
%   Q2 does the published generator's series (gen_FlxPin_allmodels_20261023
%      lines 55-72: z with NaN->0, Xp=linspace(0,620,20)(2:20) pressure axis,
%      Yq=strain/KMAX, spline interp2, Mo=-mA*F) reproduce the published
%      Hunt GoF column (RMSE 5.779/6.633/5.315/3.945)?
%   Q3 what does DIRECTLY calling HuntEq give (equation of record), with
%      epsilon_max = 0.171 (the equation's own basis) vs the per-actuator
%      KMAX; how do the series and GoF differ from the published column?
%   Q4 per-model-domain vs common-domain GoF, and the valid-point mask.
% Outputs: printed audit + Hunt_series_audit_20261026.csv.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');

%% Q1: the lookup table
Z = load(fullfile(root, 'Code/Matlab/Functions/ForceStrainForFit.mat'));
fprintf('ForceStrainForFit.mat variables: %s\n', strjoin(fieldnames(Z).', ', '));
z = Z.z;
fprintf('z size: %d x %d; finite %d; NaN %d; min %.4f; max %.4f\n', ...
    size(z,1), size(z,2), sum(isfinite(z), 'all'), sum(isnan(z), 'all'), ...
    min(z, [], 'omitnan'), max(z, [], 'omitnan'));
% NaN pattern: per column, first NaN row
firstNaN = sum(isnan(z), 1);
fprintf('NaN count per column (1..19): '); fprintf('%3d', firstNaN); fprintf('\n');
zval = z; zval(isnan(zval)) = 0;
fprintf('row 1 (eps=0) across cols:  '); fprintf('%8.2f', zval(1,:)); fprintf('\n');
fprintf('row 30 (eps=emax) across:   '); fprintf('%8.2f', zval(30,:)); fprintf('\n');
fprintf('col 19 (P=620) down rows 1..30 step 3: ');
fprintf('%8.2f', zval(1:3:30,19)); fprintf('\n');

% fresh HuntEq table on the generator grid (evalc swallows fzero chatter)
eps_g = linspace(0, 0.171, 30)';
P_g = linspace(0, 620, 19);
evalc('B = HuntEq(eps_g, 0.171, P_g, 0);');
Af = B(:, 2:end);
fprintf('\nfresh HuntEq table: size %dx%d; min %.2f; max %.2f\n', ...
    size(Af,1), size(Af,2), min(Af,[],'omitnan'), max(Af,[],'omitnan'));
nz = isnan(z);
diffRaw = abs(z - Af);
fprintf('z vs fresh HuntEq (S=0): max abs diff over common finite cells = %.4f\n', ...
    max(diffRaw(~nz), [], 'omitnan'));
FmaxGen = max(Af(1,:));
fprintf('generator Fmax = max(A(row eps=0)) = %.2f N\n', FmaxGen);
fprintf('z/FreshNorm max diff = %.4f (z == fresh/Fmax?)\n', ...
    max(abs(z(~nz) - Af(~nz)/max(Af(1,:))), [], 'omitnan'));

%% baseline evaluator (zero-correction) for the four cross-validation tests
[f0, bpa0] = minimizeFlxPin(0, Inf, Inf, 1:5);
ordr = [1 5 2 4];      % family indices in panel order 48.5, 41.0, 45.7, 40.6t
names = {'48.5', '41.0', '45.7', '40.6t'};
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    KMAX = (b.rest - b.Kmax)/b.rest;
    fprintf('test %s: rest %.4f m Kmax %.4f -> KMAX %.4f; P unique: %s; strain range [%.4f %.4f]\n', ...
        names{k}, b.rest, b.Kmax, KMAX, mat2str(unique(round(b.P,2))), ...
        min(b.strain), max(b.strain));
end

%% Q2: reproduce the published lookup series (verbatim generator code)
load(fullfile(root, 'Code/Matlab/Functions/ForceStrainForFit.mat'), 'z');
z(isnan(z)) = 0;
Xp = linspace(0, 620, 20); Xp = Xp(2:20);
Yp = linspace(0, 1, 30);
Mlookup = cell(4,1);
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    Yq = b.strain ./ ((b.rest - b.Kmax) / b.rest);
    Xq = b.P .* ones(size(Yq));
    Vq = interp2(Xp, Yp', z, Xq, Yq, 'spline');
    Fold = Vq .* b.unitD;
    Fq = hypot(Fold(:,1), Fold(:,2));
    Mo = -b.mA .* Fq;
    Mo(Yq > 1) = 0;
    Mo(b.strain < 0) = NaN;
    Mlookup{k} = Mo;
end

%% Q3: direct HuntEq series
Mdir171 = cell(4,1); MdirKMAX = cell(4,1);
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    KMAX = (b.rest - b.Kmax)/b.rest;
    eps = b.strain(:);
    % equation basis epsilon_max = 0.171
    evalc('B1 = HuntEq(eps, 0.171, b.P, 0);');
    F1 = B1(:,2);
    Mo1 = -b.mA .* F1;
    Mo1(b.strain > 0.171) = 0; Mo1(b.strain < 0) = NaN;
    Mdir171{k} = Mo1;
    % per-actuator basis epsilon_max = KMAX
    B2 = HuntEq(eps, KMAX, b.P, 0);
    F2 = B2(:,2);
    Mo2 = -b.mA .* F2;
    Mo2(b.strain > KMAX) = 0; Mo2(b.strain < 0) = NaN;
    MdirKMAX{k} = Mo2;
end

%% GoF with the evaluator's own SSE machinery (prediction at Ak -> Aexp)
published = [5.779 1.1761 9.798; 6.633 3.2313 11.409; 5.315 0.9408 9.991; 3.945 0.6750 7.527];
fprintf('\n=== Hunt-series GoF, per-model domain (RMSE / FVU / Max) ===\n');
gofLookup = zeros(4,3); gof171 = zeros(4,3); gofK = zeros(4,3);
commonMaskAll = cell(4,1);
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    fprintf('%s: published [%.3f %.4f %.3f]\n', names{k}, published(k,1), published(k,2), published(k,3));
    gofLookup(k,:) = gofMasked(b, Mlookup{k}, 'lookup');
    gof171(k,:)    = gofMasked(b, Mdir171{k}, 'dir .171');
    gofK(k,:)      = gofMasked(b, MdirKMAX{k}, 'dir KMAX');
end

%% Q4: common-domain scores (Hunt valid AND char valid AND hybrid defined)
fprintf('\n=== common-domain GoF (Hunt epsmax=0.171 vs char-alone vs hybrid) ===\n');
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    % char-alone prediction b.M has NaN where strain < -0.03
    charValid = isfinite(b.M);
    huntValid = isfinite(Mdir171{k});
    % the hybrid lives at its own angles A_h (same count as Aexp here)
    [Ak_s, ia] = sort(b.Ak);
    [Ah_s, ih] = sort(b.A_h);
    hybValid = isfinite(b.M_h);
    sH = Mdir171{k}(ia); sH = sH(huntValid);
    sC = b.M(ia);        sC = sC(charValid);
    sY = b.M_h(ih);      sY = sY(hybValid);
    aexp = b.Aexp(:);
    % unique grids for the interpolants
    [xH, ~, jH] = unique(Ak_s(huntValid), 'stable'); vH = accumarray(jH, sH, [], @mean);
    [xC, ~, jC] = unique(Ak_s(charValid), 'stable'); vC = accumarray(jC, sC, [], @mean);
    [xY, ~, jY] = unique(Ah_s(hybValid), 'stable');  vY = accumarray(jY, sY, [], @mean);
    Gi = griddedInterpolant(xH, vH);
    Gc = griddedInterpolant(xC, vC);
    Gh = griddedInterpolant(xY, vY);
    inH = aexp >= min(xH) & aexp <= max(xH);
    inC = aexp >= min(xC) & aexp <= max(xC);
    inHb = aexp >= min(xY) & aexp <= max(xY);
    common = inH & inC & inHb;
    commonMaskAll{k} = common;
    yc = b.Mexp(common); yc = yc(:);
    [rH, fH, mH] = gofManual(yc, Gi(aexp(common)));
    [rC, fC, mC] = gofManual(yc, Gc(aexp(common)));
    [rB, fB, mB] = gofManual(yc, Gh(aexp(common)));
    fprintf('%s: n_common %d of %d\n', names{k}, sum(common), numel(aexp));
    fprintf('   Hunt   %7.4f %8.4f %7.4f\n', rH, fH, mH);
    fprintf('   char   %7.4f %8.4f %7.4f\n', rC, fC, mC);
    fprintf('   hybrid %7.4f %8.4f %7.4f\n', rB, fB, mB);
    % per-model (own-domain) char for reference
    gofMasked(b, b.M, 'charOwn');
end

%% CSV: valid mask + series values
fid = fopen(fullfile(qcdir, 'Hunt_series_audit_20261026.csv'), 'w');
fprintf(fid, ['test,idx,Aexp_deg,measured_Nm,P_kPa,', ...
    'lookup_published_atAexp,direct_epsmax0171_atAexp,direct_epsmaxKMAX_atAexp,', ...
    'char_alone_atAexp,hybrid_atAexp,hunt_valid_atAexp,common_valid\n']);
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    KMAX = (b.rest - b.Kmax)/b.rest;
    [Ak_s, ia] = sort(b.Ak);
    [Ah_s, ih] = sort(b.A_h);
    charValid = isfinite(b.M); huntValid = isfinite(Mdir171{k});
    hybValid = isfinite(b.M_h);
    sH = Mdir171{k}(ia); sH = sH(huntValid);
    sC = b.M(ia);        sC = sC(charValid);
    sY = b.M_h(ih);      sY = sY(hybValid);
    [xH, ~, jH] = unique(Ak_s(huntValid), 'stable'); vHu = accumarray(jH, sH, [], @mean);
    [xC, ~, jC] = unique(Ak_s(charValid), 'stable'); vCu = accumarray(jC, sC, [], @mean);
    [xY, ~, jY] = unique(Ah_s(hybValid), 'stable');  vYu = accumarray(jY, sY, [], @mean);
    Gi = griddedInterpolant(xH, vHu);
    Gc = griddedInterpolant(xC, vCu);
    Gh = griddedInterpolant(xY, vYu);
    aexp = b.Aexp(:);
    inH = aexp >= min(xH) & aexp <= max(xH);
    inC = aexp >= min(xC) & aexp <= max(xC);
    inHb = aexp >= min(xY) & aexp <= max(xY);
    common = inH & inC & inHb;
    for p = 1:numel(aexp)
        ae = aexp(p);
        vL = interp1(Ak_s, Mlookup{k}(ia), ae, 'linear', NaN);
        v1 = interp1(Ak_s, Mdir171{k}(ia), ae, 'linear', NaN);
        vK = interp1(Ak_s, MdirKMAX{k}(ia), ae, 'linear', NaN);
        vC = interp1(Ak_s, b.M(ia), ae, 'linear', NaN);
        vY = interp1(Ah_s, b.M_h(ih), ae, 'linear', NaN);
        fprintf(fid, '%s,%d,%.3f,%.4f,%.1f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d\n', ...
            names{k}, p, ae, b.Mexp(p), b.P, vL, v1, vK, vC, vY, ...
            double(inH(p)), double(common(p)));
    end
end
fclose(fid);
fprintf('WROTE Hunt_series_audit_20261026.csv\n');
disp('DONE_audit_Hunt');

%% local: manual GoF triplet (Go_OfF formulas, verified in dbg_Hunt)
function [r, fv, mx] = gofManual(y, v)
y = y(:); v = v(:);
resid = y - v;
n = numel(y);
SSresid = sum(resid.^2);
SStot = sum((y - mean(y)).^2);
r = sqrt(SSresid/n);
fv = SSresid*(n-1)/(SStot*n);
mx = max(abs(resid));
end

%% local: GoF of a masked series on its own valid domain (SSE-style).
% Metrics replicate Go_OfF.m exactly (verified numerically in dbg_Hunt):
% RMSE = sqrt(SSresid/n); FVU = SSresid*(n-1)/(SStot*n); Max = max|resid|.
function out = gofMasked(b, series, label)
[Ak_s, ia] = sort(b.Ak);
sv = series(ia);
valid = isfinite(sv);
Ak_v = Ak_s(valid); sv = sv(valid);
[Ak_u, ~, iu] = unique(Ak_v, 'stable');   % griddedInterpolant needs unique x
sv_u = accumarray(iu, sv, [], @mean);
G = griddedInterpolant(Ak_u, sv_u);
ae = b.Aexp(:);
sel = ae >= min(Ak_u) & ae <= max(Ak_u);
y = b.Mexp(sel); y = y(:);
v = G(ae(sel)); v = v(:);
resid = y - v;
n = numel(y);
SSresid = sum(resid.^2);
SStot = sum((y - mean(y)).^2);
r = sqrt(SSresid/n);
fv = SSresid*(n-1)/(SStot*n);
mx = max(abs(resid));
out = [r fv mx];
if nargin >= 3
    fprintf('   %-9s [%8.4f %8.4f %8.4f]  n=%d\n', label, r, fv, mx, n);
end
end
