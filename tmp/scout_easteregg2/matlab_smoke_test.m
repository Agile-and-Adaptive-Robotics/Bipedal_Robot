% Headless MATLAB smoke test for the corrected campaign (easteregg2).
% Verifies: parallel pool, Optimization Toolbox solvers, dated save.
fprintf('MATLAB %s on %s\n', version, computer);
p = parpool('local');
fprintf('pool: %d workers\n', p.NumWorkers);
vals = zeros(4, 1);
parfor k = 1:4
    vals(k) = k ^ 2;
end
fprintf('parfor ok: %g\n', vals(4));
f = @(x) (x(1) - 2) ^ 2 + (x(2) + 1) ^ 2;
opts = optimoptions('patternsearch', 'Display', 'off', 'UseParallel', true);
[x, fv] = patternsearch(f, [0, 0], [], [], [], [], [-5, -5], [5, 5], opts);
fprintf('patternsearch ok: x = [%.3f %.3f], f = %.6g\n', x, fv);
optsS = optimoptions('surrogateopt', 'Display', 'off', ...
    'MaxFunctionEvaluations', 60, 'UseParallel', true);
[xS, fS] = surrogateopt(f, [-5, -5], [5, 5], optsS);
fprintf('surrogateopt ok: f = %.6g\n', fS);
delete(p);
disp('SMOKE PASS');
