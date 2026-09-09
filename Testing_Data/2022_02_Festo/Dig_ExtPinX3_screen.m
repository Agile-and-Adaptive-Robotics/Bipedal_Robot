%Dig_ExtPinX3_screen.m — screen the best biomimetic-flexor Xi combinations through
%the pinned knee extensor (minimizeExtX3, POOL = allBPA [1 2 5 6 8], NOT all 9) and
%the biomimetic knee extensor (minimizeExt, 52cm). Xi3 grid per candidate; report
%best Xi3 per candidate. 325kPa-style degenerate FVU guard: none needed here.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
POOL = [1 2 5 6 8];
labels9 = ["40cm","40cm-tendon","42cm","42cm-tendon","43cm","43cm-tendon","46cm","47cm","48cm"];

%Candidates: name | Xi0 | Xi1 | Xi2  (provenance in name)
cand = { ...
 "2trans_noT3_front_r1", 0.0027, 4.3e4, 1.11e4; ...
 "2trans_noT3_front_r2", 0.0016, 4.3e4, 1.03e4; ...
 "2trans_noT3_front_r3", 0.0033, 4.45e4, 1.18e4; ...
 "2trans_noT3_front_r4", 0.0052, 4.49e4, 1.37e4; ...
 "2trans_noT3_front_r6", 0.0089, 5.62e4, 1.85e4; ...
 "1trans_noT3_pick1",    0.0060, 4.36e5, 2.0e4; ...
 "1trans_noT3_front_r3", 0.0064, 9.99e5, 2.0e4; ...
 "handgrid_best",        0.012,  5e5,   1.0e4; ...
 "handgrid_r2",          0.010,  5e5,   8.5e3; ...
 "legacyCV_reference",   0.0,    1.508e4, 1.218e4};

xi3s = [0.02, 0.05, 0.10, 0.20, 0.35];
a0e = minimizeExtX3(0, Inf, Inf, 0);                 %9x3 baseline (pinned ext)
be0 = minimizeExt(0, Inf, Inf, 0, 1);                %1x3 baseline (bio ext)

fprintf('Baselines: pinned pool RMSE %.3f | bio-ext RMSE %.3f FVU %.3f MaxR %.3f\n', ...
    mean(a0e(POOL,1)), be0(1), be0(2), be0(3));
fprintf('%-22s %8s %9s %9s | %5s | %8s %8s | %8s %8s %8s\n', ...
    'candidate','Xi0(m)','Xi1','Xi2','Xi3','poolRMSE','poolFVU','bioRMSE','bioFVU','bioMaxR');
nc = size(cand,1); out = nan(nc, 9);
for i = 1:nc
    xi0 = -cand{i,2};   %EXTENSOR SIGN: Xi0 = negative of the flexor value (per minimizeExt10mmX3.m -g(1); bounds [-2,0] cm)
    xi1 = cand{i,3}; xi2 = cand{i,4};
    bestP = Inf; bx3 = NaN; bmet = nan(1,3); bpool = nan(1,2); bs = nan(1,2);
    for q = 1:numel(xi3s)
        f = minimizeExtX3(xi0, xi1, xi2, xi3s(q));
        R = f(POOL,1) ./ a0e(POOL,1); V = f(POOL,2) ./ a0e(POOL,2);
        s = mean([R, V]);                       %1x2: [mean RMSE ratio, mean FVU ratio]
        if prod(s) < bestP
            bestP = prod(s); bs = s; bx3 = xi3s(q);
            bpool = [mean(f(POOL,1)), mean(f(POOL,2))];
            b = minimizeExt(xi0, xi1, xi2, xi3s(q), 1);
            bmet = b(1:3);
        end
    end
    fprintf('%-22s %8.4f %9.3g %9.3g | %5.2f | %8.3f %8.3f | %8.3f %8.3f %8.3f\n', ...
        cand{i,1}, xi0, xi1, xi2, bx3, bpool(1), bpool(2), bmet(1), bmet(2), bmet(3));
    out(i,1) = xi0; out(i,2) = xi1; out(i,3) = xi2; out(i,4) = bx3;
    out(i,5) = bs(1); out(i,6) = bs(2);          %pool mean RMSE / FVU ratios vs baseline
    out(i,7) = bmet(1); out(i,8) = bmet(2); out(i,9) = bmet(3);
end
fprintf('(poolRMSE/poolFVU = means over tests [1 2 5 6 8]; ratios vs baseline are in the saved mat)\n');
save('Dig_out/Dig_ExtPinX3_screen_results.mat', 'cand', 'out', 'xi3s', 'POOL');
fprintf('Saved to Dig_out\\Dig_ExtPinX3_screen_results.mat\n');
