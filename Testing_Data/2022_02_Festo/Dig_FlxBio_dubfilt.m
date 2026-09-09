%Dig_FlxBio_dubfilt.m — pool the 8 corrected-encoder CV fronts, dedup near-repeats,
%then evaluate every unique candidate on the biomimetic knee flexor (minimizeFlx:
%straight calculation, no optimization). Ranks by combined 10mm+20mm improvement.
%Dedup rule: candidate is a near-repeat of a kept one if |dXi0|<1mm and Xi1, Xi2
%are each within 5% relative (1.234e4 vs 1.2e4 counts as the same, per Ben).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

tags = {'full','noT3','noT5','noT3noT5'};
arms = {'2trans','1trans'};
pool = zeros(0,3);
for a = 1:2
    for t = 1:4
        fn = sprintf('minimizeFlxPin10_results_20260908_2brkt_%s_%s.mat', arms{a}, tags{t});
        Sv = load(fn, 'filtered_results', 'xCols');
        pool = [pool; Sv.filtered_results(:, Sv.xCols)]; %#ok<AGROW>
    end
end

%Dedup -> dub_filt_results
dub_filt_results = zeros(0,3);
for i = 1:size(pool,1)
    c = pool(i,:);
    dup = false;
    for j = 1:size(dub_filt_results,1)
        k = dub_filt_results(j,:);
        if abs(c(1)-k(1)) < 1e-3 && abs(c(2)-k(2)) < 0.05*max(abs(k(2)),eps) && ...
           abs(c(3)-k(3)) < 0.05*max(abs(k(3)),eps)
            dup = true; break;
        end
    end
    if dup == false
        dub_filt_results(end+1,:) = c; %#ok<AGROW>
    end
end
nDup = size(pool,1) - size(dub_filt_results,1);
fprintf('Pooled %d filtered candidates -> %d unique in dub_filt_results (%d near-repeats removed)\n', ...
    size(pool,1), size(dub_filt_results,1), nDup);

%% Baselines on the biomimetic flexor
[b10, b20a, b20b] = minimizeFlx(0, Inf, Inf);
bR = [b10(1), b20a(1), b20b(1)]; bF = [b10(2), b20a(2), b20b(2)];
fprintf('Biomimetic baselines (RMSE): 10mm %.3f | 20mm-a %.3f | 20mm-b %.3f\n\n', bR(1), bR(2), bR(3));

%% Evaluate all unique candidates (serial: minimizeFlx loads data via relative paths,
% which parallel workers cannot see)
n = size(dub_filt_results,1);
score = Inf(n,1); met = nan(n,9);   %[rmse fv maxR] x [10mm, 620kPa, 325kPa]
nErr = 0; firstErr = '';
for i = 1:n
    try
        [u, va, vb] = minimizeFlx(dub_filt_results(i,1), dub_filt_results(i,2), dub_filt_results(i,3));
        met(i,:) = [u(1:3), va(1:3), vb(1:3)];
        R = [u(1), va(1), vb(1)] ./ bR;
        V = [u(2), va(2), vb(2)] ./ bF;
        %325kPa baseline FVU is 0 (degenerate) -> exclude that one ratio from the score
        score(i) = mean([R, V(1:2)]);
    catch MEi
        nErr = nErr + 1;
        if isempty(firstErr), firstErr = MEi.message; end
    end
    if mod(i, 50) == 0, fprintf('  %d/%d evaluated\n', i, n); end
end
fprintf('Evaluation done: %d ok, %d errors. %s\n', n - nErr, nErr, firstErr);
[~, ord] = sort(score, 'ascend');

fprintf('Top 8 unique results on the biomimetic knee (score = mean of RMSE+FVU ratios vs baseline; 1.00 = baseline):\n');
fprintf('%-4s %8s %9s %9s | %33s | %33s | %s\n', 'rank','Xi0(m)','Xi1','Xi2', ...
    '10mm: RMSE  FVU  MaxResid', '20a: RMSE  FVU  MaxResid', 'score  prior');
shown = 0;
for q = 1:n
    i = ord(q);
    if ~isfinite(score(i)), continue; end
    shown = shown + 1;
    prior = (dub_filt_results(i,1) >= 0.005 && dub_filt_results(i,1) <= 0.010 && ...
             dub_filt_results(i,2) > 10*dub_filt_results(i,3));
    fprintf('%-4d %8.4f %9.3g %9.3g | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %6.3f %6d\n', ...
        q, dub_filt_results(i,1), dub_filt_results(i,2), dub_filt_results(i,3), ...
        met(i,1), met(i,2), met(i,3), met(i,4), met(i,5), met(i,6), ...
        met(i,7), met(i,8), met(i,9), score(i), prior);
    if shown == 8, break; end
end
fprintf('(BenPrior = 1 if Xi0 in [5,10] mm and Xi1 > 10*Xi2, per Ben intuition)\n');

save('Dig_out/Dig_FlxBio_dubfilt_results.mat', 'dub_filt_results', 'score', 'met', 'ord', 'bR', 'bF');
fprintf('Saved dub_filt_results + rankings to Dig_out\\Dig_FlxBio_dubfilt_results.mat\n');
