%Pick sweep v2 — fixed loop/indexing, dedup rule (Ben), extra GoF measures.
%Dedup rule (Ben): a candidate is redundant if its Xi0 is within 1 mm of a kept
%candidate AND Xi1 AND Xi2 are both within +/-1e4 of it AND its distance is not a
%big jump (|dDist| < DISTJUMP). Otherwise it is kept as a distinct representative.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

TOLX0 = 0.001;   %m   (1 mm)
TOLX12 = 1e4;    %N/m (Xi1 or Xi2 separation)
DISTJUMP = 0.05; %big change in the distance calc keeps a candidate anyway

S = load('minimizeFlxPin10_results_20260907_2brkt_2trans.mat', 'filtered_results', 'xCols', 'labels');
F = S.filtered_results; xC = S.xCols;
trainCols = 7:9; valCols = 10:12; dC = 13;   %[ind, hold(2), x(3), train(3), val(3), dist]
N = size(F,1);
fprintf('Front: %d candidates. Xi0 [%.4f, %.4f] m | Xi1 [%.3g, %.3g] | Xi2 [%.3g, %.3g]\n', ...
    N, min(F(:,xC(1))), max(F(:,xC(1))), min(F(:,xC(2))), max(F(:,xC(2))), min(F(:,xC(3))), max(F(:,xC(3))));

%% Dedup per Ben's rule
kept = zeros(1,0);
for p = 1:N
    redundant = false;
    for q = kept
        if abs(F(p,xC(1))-F(q,xC(1))) < TOLX0 && ...
           abs(F(p,xC(2))-F(q,xC(2))) < TOLX12 && ...
           abs(F(p,xC(3))-F(q,xC(3))) < TOLX12 && ...
           abs(F(p,dC)-F(q,dC)) < DISTJUMP
            redundant = true; break;
        end
    end
    if redundant == false
        kept(end+1) = p; %#ok<AGROW>
    end
end
K = numel(kept);
fprintf('Dedup: %d -> %d distinct representatives (TOLX0=1mm, TOLX12=1e4, DISTJUMP=%.2f)\n\n', N, K, DISTJUMP);

%% Table of representatives with standard GoF (stored columns)
fprintf('%-5s %8s %10s %10s %8s %8s %8s %8s\n', 'pick', 'Xi0(m)', 'Xi1', 'Xi2', 'dist', 'trnRMSE', 'valRMSE', 'valFVU');
for q = 1:K
    p = kept(q);
    fprintf('%-5d %8.4f %10.3g %10.3g %8.3f %8.3f %8.3f %8.3f\n', ...
        p, F(p,xC(1)), F(p,xC(2)), F(p,xC(3)), F(p,dC), ...
        mean(F(p,trainCols)), F(p,valCols(1)), F(p,valCols(2)));
end

%% Extra GoF per representative: R^2, MAE, mean |Lm_p - Lm_h| (all 5 tests)
nEval = min(K, 15);
res = nan(nEval, 4);   %[R2, MAE, dLm_mm, meanRMSE]
fprintf('\nExtra GoF for the first %d representatives (full 5-test evaluation):\n', nEval);
fprintf('%-5s %8s %10s %10s %7s %7s %8s %8s\n', 'pick', 'Xi0(m)', 'Xi1', 'Xi2', 'R^2', 'MAE', 'dLm(mm)', 'meanRMSE');
for q = 1:nEval
    p = kept(q);
    [f, bpa] = minimizeFlxPin2brk(F(p,xC(1)), F(p,xC(2)), F(p,xC(3)), [], true, '2trans');
    r2v = zeros(5,1); maev = zeros(5,1); dlm = zeros(5,1);
    for j = 1:5
        [aks, isrt] = sort(bpa(j).Ak);
        Mp = interp1(aks, bpa(j).M_p(isrt,3), bpa(j).Aexp);
        e = Mp - bpa(j).Mexp;
        r2v(j) = 1 - sum(e.^2)/sum((bpa(j).Mexp - mean(bpa(j).Mexp)).^2);
        maev(j) = mean(abs(e));
        Lm_p = bpa(j).Lmt_p - 2*bpa(j).fitn - bpa(j).ten;
        Lm_h = bpa(j).Lmt - 2*bpa(j).fitn - bpa(j).ten;
        Lmh = interp1(bpa(j).Ak, Lm_h, bpa(j).Aexp);
        dlm(j) = mean(abs(Lm_p(isrt)' - Lmh))*1000;  %mm
    end
    res(q,:) = [mean(r2v), mean(maev), mean(dlm), mean(f,1,'omitnan')];
    fprintf('%-5d %8.4f %10.3g %10.3g %7.3f %7.3f %8.2f %8.3f\n', ...
        p, F(p,xC(1)), F(p,xC(2)), F(p,xC(3)), res(q,1), res(q,2), res(q,3), res(q,4));
end
fprintf('(R^2/MAE from torque vs measured at Aexp; dLm = mean |Lm_p - Lm_h| in mm)\n');

%% Nearest candidate to the historical 1trans-era pick
ref = [0.0075, 5.41e4, 2.67e4];
dref = vecnorm(F(:,xC) - ref, 2, 2) ./ vecnorm(ref, 2, 2);
[rmin, imin] = min(dref);
fprintf('\nNearest front candidate to (7.5mm, 5.41e4, 2.67e4): rank %d, Xi0=%.4f Xi1=%.3g Xi2=%.3g (rel. dist %.3f)\n', ...
    imin, F(imin,xC(1)), F(imin,xC(2)), F(imin,xC(3)), rmin);
