%Dig_FlxBio_refine2.m — extended flexor refinement grid: the first refinement hit
%grid edges at Xi0=12mm and Xi1=5e5, so this pushes Xi0 to 16mm and Xi1 to 1e6.
%Same score as before: mean RMSE ratios (10/620/325) + FVU ratios (10/620).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

Dg = load('Dig_out/Dig_FlxBio_dubfilt_results.mat', 'bR', 'bF');
bR = Dg.bR; bF = Dg.bF;

xi0 = [0.010, 0.011, 0.012, 0.013, 0.014, 0.016];
xi1 = [2e5, 5e5, 1e6];
xi2 = [8.5e3, 1e4, 1.2e4];
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
fprintf('Extended flexor refinement: %d points, %d errors. Best first:\n', nG, nErr);
fprintf('%-5s %8s %9s %9s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %6s\n', ...
    'rank','Xi0(m)','Xi1','Xi2','R10','F10','M10','R620','F620','M620','R325','F325','M325','score');
for q = 1:min(12, nG)
    i = ord(q);
    fprintf('%-5d %8.4f %9.3g %9.3g | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %6.3f\n', ...
        q, grid(i,1), grid(i,2), grid(i,3), met(i,1), met(i,2), met(i,3), ...
        met(i,4), met(i,5), met(i,6), met(i,7), met(i,8), met(i,9), score(i));
end
save('Dig_out/Dig_FlxBio_refine2_results.mat', 'grid', 'met', 'score', 'xi0', 'xi1', 'xi2');
fprintf('Saved to Dig_out\\Dig_FlxBio_refine2_results.mat\n');
