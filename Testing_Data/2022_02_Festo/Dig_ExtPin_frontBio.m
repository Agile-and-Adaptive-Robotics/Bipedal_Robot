%Dig_ExtPin_frontBio.m — scan the captured noT3newXi extensor front through the
%biomimetic knee extensor (minimizeExt, 52cm). Dedup per Ben's rule; score = mean
%RMSE ratio + FVU ratio vs baseline (325kPa-style degeneracy: none here). Xi3 shown
%per candidate so the Xi3 > 0.02 view is direct.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

MATFILE = 'minimizeExt10mmX3_results_20260910_noT3.mat';
S = load(MATFILE, 'filtered_results', 'xCols', 'sol_actual', 'f');
F = S.filtered_results; xC = S.xCols; N = size(F,1);
be0 = minimizeExt(0, Inf, Inf, 0, 1);
fprintf('Extensor front: %d rows | pick=1: Xi0=%.4f Xi1=%.3g Xi2=%.3g Xi3=%.3f\n', ...
    N, S.sol_actual(1), S.sol_actual(2), S.sol_actual(3), S.sol_actual(4));
fprintf('Bio-ext baseline: RMSE %.3f FVU %.3f MaxR %.3f\n', be0(1), be0(2), be0(3));

% dedup
dub = zeros(0,4); dubRow = zeros(0,1);
for i = 1:N
    c = F(i,xC);
    dup = false;
    for j = 1:size(dub,1)
        k = dub(j,:);
        if abs(c(1)-k(1)) < 1e-3 && abs(c(2)-k(2)) < 0.05*max(abs(k(2)),eps) && ...
           abs(c(3)-k(3)) < 0.05*max(abs(k(3)),eps)
            dup = true; break;
        end
    end
    if dup == false
        dub(end+1,:) = c; dubRow(end+1) = i; %#ok<AGROW>
    end
end
nd = size(dub,1);
fprintf('Dedup: %d -> %d unique candidates\n', N, nd);

met = nan(nd,3); score = Inf(nd,1);
for i = 1:nd
    try
        b = minimizeExt(dub(i,1), dub(i,2), dub(i,3), dub(i,4), 1);
        met(i,:) = b(1:3);
        score(i) = mean([met(i,1)/be0(1), met(i,2)/be0(2)]);
    catch
        score(i) = Inf;
    end
end
[~, ord] = sort(score, 'ascend');

pickScore = NaN; pickRank = NaN;
for q = 1:nd
    if dubRow(ord(q)) == 1, pickRank = q; pickScore = score(ord(q)); break; end
end

fprintf('===== Bio-extensor results for the extensor front (best 12) =====\n');
fprintf('%-4s %-5s %8s %9s %9s %6s | %5s %5s %5s | %6s\n', ...
    'rank','front','Xi0(m)','Xi1','Xi2','Xi3','R','FVU','MaxR','score');
shown = 0; nAbove = 0;
for q = 1:nd
    i = ord(q);
    if ~isfinite(score(i)), continue; end
    shown = shown + 1;
    if dub(i,4) > 0.02, nAbove = nAbove + 1; end
    fprintf('%-4d %-5d %8.4f %9.3g %9.3g %6.3f | %5.2f %5.2f %5.2f | %6.3f\n', ...
        q, dubRow(i), dub(i,1), dub(i,2), dub(i,3), dub(i,4), ...
        met(i,1), met(i,2), met(i,3), score(i));
    if shown == 12, break; end
end
fprintf('Of the top %d shown, %d have Xi3 > 0.02\n', shown, nAbove);
save('Dig_out/Dig_ExtPin_frontBio_results.mat', 'dub', 'dubRow', 'met', 'score');
fprintf('Saved to Dig_out\\Dig_ExtPin_frontBio_results.mat\n');
