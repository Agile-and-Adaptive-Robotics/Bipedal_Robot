% Smoke test 1: evaluator spot-check for minimizeFlxPin2brk
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');

t0 = tic;
[f_inf, ~] = minimizeFlxPin2brk(0, Inf, Inf, [], true);
fprintf('2brk baseline (0,Inf,Inf): OK in %.1fs\n', toc(t0));
disp(array2table(f_inf, 'VariableNames', {'RMSE','FVU','MaxResidual'}));

t0 = tic;
[f_b1, ~] = minimizeFlxPin2brk(0.01, 2e4, 0.8e4, [], false);
fprintf('2brk sol (0.01, 2e4, 0.8e4), bracket2 OFF: OK in %.1fs\n', toc(t0));
disp(array2table(f_b1, 'VariableNames', {'RMSE','FVU','MaxResidual'}));

t0 = tic;
[f_b2, bpa2] = minimizeFlxPin2brk(0.01, 2e4, 0.8e4, [], true);
fprintf('2brk sol (0.01, 2e4, 0.8e4), bracket2 ON: OK in %.1fs\n', toc(t0));
disp(array2table(f_b2, 'VariableNames', {'RMSE','FVU','MaxResidual'}));

% Origin-bracket deflection sanity (should be nonzero with bracket2 on)
eA2max = cellfun(@(b) max(abs(b.eA2(:))), num2cell(bpa2));
fprintf('max |eA2| per BPA (m): %s\n', mat2str(eA2max, 3));

% Compare against the original single-bracket evaluator
[f_old, ~] = minimizeFlxPin(0.01, 2e4, 0.8e4);
disp(table(f_old(:,1), f_b1(:,1), f_b2(:,1), f_old(:,2), f_b1(:,2), f_b2(:,2), ...
    'VariableNames', {'RMSE_orig','RMSE_2brk_noB2','RMSE_2brk_B2','FVU_orig','FVU_2brk_noB2','FVU_2brk_B2'}));

assert(all(isfinite(f_inf(:))), 'nonfinite baseline metrics');
assert(all(isfinite(f_b1(:))), 'nonfinite metrics with bracket2 off');
assert(all(isfinite(f_b2(:))), 'nonfinite metrics with bracket2 on');
assert(any(eA2max > 0), 'bracket2 deflections all zero - bracket 2 not active?');
fprintf('SMOKE1 PASS\n');
