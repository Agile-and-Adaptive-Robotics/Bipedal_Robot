%Direct pinned-data check: does the refinement corner (and the 2trans front-best)
%beat baseline on the noT3 pinned tests? Straight calc, no optimizer.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
tests = [1 2 4 5];
labels = ["48cm"; "46cm"; "40cm-tendon"; "41cm"];
a0 = minimizeFlxPin2brk(0, Inf, Inf, tests, true, '2trans');
fw = minimizeFlxPin2brk(0.012, 5e5, 1e4, tests, true, '2trans');      %refinement corner
fb = minimizeFlxPin2brk(0.0027, 4.3e4, 1.11e4, tests, true, '2trans'); %2trans front-best
fprintf('%-12s %14s %17s %17s\n', 'test', 'base RMSE FVU', 'winner RMSE FVU', 'frontbest RMSE FVU');
for j = 1:4
    r = tests(j);   %output rows are indexed by kf test number, not loop position
    fprintf('%-12s %8.2f %6.2f %10.2f %6.2f %12.2f %6.2f\n', labels(j), a0(r,1), a0(r,2), fw(r,1), fw(r,2), fb(r,1), fb(r,2));
end
tw = tests; fb2 = tests;
fprintf('winner beats baseline on all 4 (RMSE+FVU): %d\n', all(fw(tw,1:2) <= a0(tw,1:2), 'all'));
fprintf('front-best beats baseline on all 4 (RMSE+FVU): %d\n', all(fb(fb2,1:2) <= a0(fb2,1:2), 'all'));
