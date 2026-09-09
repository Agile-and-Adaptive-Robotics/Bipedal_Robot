%Dig_FlxPin_frontScan.m — scan the deduped filtered fronts of the noT3 CVs
%(1trans and 2trans) through the biomimetic flexor calc. Top 8 per front, all 3 GoF
%per case; score = mean of RMSE ratios (10/620/325) + FVU ratios (10/620 only,
%325 baseline FVU degenerate). BenPrior = Xi0 in [5,10] mm and Xi1 > 10*Xi2.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

mats  = {'minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat', ...
         'minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat'};
names = {'1trans_noT3', '2trans_noT3'};

[b10, b620, b325] = minimizeFlx(0, Inf, Inf);
bR = [b10(1), b620(1), b325(1)]; bF = [b10(2), b620(2), b325(2)];

for m = 1:2
    S = load(mats{m}, 'filtered_results', 'xCols', 'k1', 'k2', 'k3', 'f');
    F = S.filtered_results; xC = S.xCols; N = size(F,1);

    dub = zeros(0,3); dubRow = zeros(0,1);
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

    met = nan(nd,9); score = Inf(nd,1);
    for i = 1:nd
        [u, v620, v325] = minimizeFlx(dub(i,1), dub(i,2), dub(i,3));
        met(i,:) = [u(1:3), v620(1:3), v325(1:3)];
        R = met(i,[1 4 7]) ./ bR; V = met(i,[2 5 8]) ./ bF;
        score(i) = mean([R, V(1:2)]);
    end
    [~, ord] = sort(score, 'ascend');

    pickScore = NaN; pickRank = NaN;
    for q = 1:nd
        if dubRow(ord(q)) == 1, pickRank = q; pickScore = score(ord(q)); break; end
    end

    fprintf('===== %s: %d filtered -> %d unique | pick=1 biomimetic score %.3f (rank %d of %d) =====\n', ...
        names{m}, N, nd, pickScore, pickRank, nd);
    fprintf('%-4s %-5s %8s %9s %9s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %6s %5s\n', ...
        'rank','front','Xi0(m)','Xi1','Xi2','R10','F10','M10','R620','F620','M620','R325','F325','M325','score','prior');
    shown = 0;
    for q = 1:nd
        i = ord(q);
        if ~isfinite(score(i)), continue; end
        shown = shown + 1;
        prior = (dub(i,1) >= 0.005 && dub(i,1) <= 0.010 && dub(i,2) > 10*dub(i,3));
        fprintf('%-4d %-5d %8.4f %9.3g %9.3g | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %6.3f %5d\n', ...
            q, dubRow(i), dub(i,1), dub(i,2), dub(i,3), ...
            met(i,1), met(i,2), met(i,3), met(i,4), met(i,5), met(i,6), ...
            met(i,7), met(i,8), met(i,9), score(i), prior);
        if shown == 8, break; end
    end
    fprintf('\n');
end
fprintf('Biomimetic baselines (RMSE): 10mm %.3f | 620 %.3f | 325 %.3f\n', bR);
