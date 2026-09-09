%Dig_FlxBio_refine.m — refinement grid around the hand-tune winner on the
%biomimetic flexor (straight calc): Xi0 8..12 mm, Xi1 5e4..5e5, Xi2 7e3..1.4e4.
%Score = mean of RMSE ratios (10mm/620/325) + FVU ratios (10mm/620 only —
%325kPa baseline FVU is degenerate). All 3 GoF printed.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

Dg = load('Dig_out/Dig_FlxBio_dubfilt_results.mat', 'bR', 'bF');
bR = Dg.bR; bF = Dg.bF;

xi0 = [0.008, 0.009, 0.010, 0.011, 0.012];
xi1 = [5e4, 1e5, 2e5, 5e5];
xi2 = [7e3, 8.5e3, 1e4, 1.2e4, 1.4e4];
nG = numel(xi0)*numel(xi1)*numel(xi2);
grid = zeros(nG, 3); gIdx = 0;
for a = 1:numel(xi0)
    for b = 1:numel(xi1)
        for c = 1:numel(xi2)
            gIdx = gIdx + 1;
            grid(gIdx,:) = [xi0(a), xi1(b), xi2(c)];
        end
    end
end

met = nan(nG, 9); score = Inf(nG,1); nErr = 0;
for i = 1:nG
    try
        [u, v620, v325] = minimizeFlx(grid(i,1), grid(i,2), grid(i,3));
        met(i,:) = [u(1:3), v620(1:3), v325(1:3)];
        R = met(i,[1 4 7]) ./ bR; V = met(i,[2 5 8]) ./ bF;
        score(i) = mean([R, V(1:2)]);
    catch
        nErr = nErr + 1;
    end
end
[~, ord] = sort(score, 'ascend');
fprintf('Refinement grid: %d points, %d errors. Best first (score = mean RMSE ratios + 10/620 FVU ratios):\n', nG, nErr);
fprintf('%-5s %8s %9s %9s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %6s\n', ...
    'rank','Xi0(m)','Xi1','Xi2','R10','F10','M10','R620','F620','M620','R325','F325','M325','score');
for q = 1:min(12, nG)
    i = ord(q);
    fprintf('%-5d %8.4f %9.3g %9.3g | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %6.3f\n', ...
        q, grid(i,1), grid(i,2), grid(i,3), met(i,1), met(i,2), met(i,3), ...
        met(i,4), met(i,5), met(i,6), met(i,7), met(i,8), met(i,9), score(i));
end
save('Dig_out/Dig_FlxBio_refine_results.mat', 'grid', 'met', 'score', 'xi0', 'xi1', 'xi2');
fprintf('Saved to Dig_out\\Dig_FlxBio_refine_results.mat\n');
