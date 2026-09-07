function [winner, stats] = abSolver()
%ABSOLVER One-fold A/B on the flexor 2brk problem: gamultiobj vs surrogateopt.
% gamultiobj: Pop 50 x Gen 100 (~5000 evals), 3-objective Pareto
% surrogateopt: 1500 evaluations of the baseline-normalized weighted sum
% Winner: surrogateopt if clearly faster and not much worse on holdout distance.
useB2 = true;
train = [2 3 4 5]; holdout = 1;
W = [0.5, 0.3, 0.2];

    function ff = abmin1(x)
        Xi0 = x(1)/100; Xi1 = 10^x(2); Xi2 = 10^x(3);
        try
            f_all = minimizeFlxPin2brk(Xi0, Xi1, Xi2, train, useB2);
            ff = mean(f_all(train,:) ./ a0(train,:), 1, 'omitnan');
        catch
            ff = [Inf, Inf, Inf];
        end
    end

    function fs = abmin1scalar(x)
        f3 = abmin1(x);
        if any(~isfinite(f3)), fs = Inf; else, fs = W * f3(:); end
    end

a0 = minimizeFlxPin2brk(0, Inf, Inf, [], useB2);
lb = [0*100, log10(5e3), log10(5e3)];
ub = [0.020*100, log10(5e7), log10(5e7)];

%gamultiobj
t1 = tic;
opts = optimoptions('gamultiobj', 'UseParallel', true, 'Display', 'off', ...
    'PopulationSize', 50, 'MaxGenerations', 100, ...
    'MutationFcn', {@mutationadaptfeasible}, 'CrossoverFraction', 0.8, ...
    'CrossoverFcn', {@crossoverscattered}, 'FunctionTolerance', 4e-3);
[xg, ~] = gamultiobj(@abmin1, 3, [], [], [], [], lb, ub, opts);
stats.gaTime = toc(t1);
valg = zeros(size(xg,1), 3);
parfor i = 1:size(xg,1)
    Xi0 = xg(i,1)/100; Xi1 = 10^xg(i,2); Xi2 = 10^xg(i,3);
    f_all = minimizeFlxPin2brk(Xi0, Xi1, Xi2, holdout, useB2);
    valg(i,:) = mean(f_all(holdout,:) ./ a0(holdout,:), 2);
end
stats.gaDist = min(vecnorm(valg, 2, 2));
fprintf('gamultiobj : %.0f s, best holdout distance %.4f (%d pareto pts)\n', ...
    stats.gaTime, stats.gaDist, size(xg,1));

%surrogateopt
t2 = tic;
surrOpts = optimoptions('surrogateopt', 'UseParallel', true, 'Display', 'off', ...
    'MaxFunctionEvaluations', 1500);
sol = surrogateopt(@abmin1scalar, lb, ub, surrOpts);
stats.surrTime = toc(t2);
%R2025a: first surrogateopt output IS the solution point (empty if all evals fail)
assert(~isempty(sol), 'surrogateopt returned empty - all evaluations failed');
xs = reshape(sol, 1, []);
f_all = minimizeFlxPin2brk(xs(1)/100, 10^xs(2), 10^xs(3), holdout, useB2);
stats.surrDist = norm(mean(f_all(holdout,:) ./ a0(holdout,:), 1));
fprintf('surrogateopt: %.0f s, holdout distance %.4f\n', stats.surrTime, stats.surrDist);

%Decision: take surrogateopt if it is clearly faster and not much worse
if stats.surrDist <= 1.5 * stats.gaDist && stats.surrTime < stats.gaTime
    winner = 'surrogateopt';
else
    winner = 'gamultiobj';
end
fprintf('A/B winner: %s\n', winner);
end
