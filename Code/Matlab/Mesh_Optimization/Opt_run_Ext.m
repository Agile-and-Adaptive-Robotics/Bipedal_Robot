clear
clc
rehash

ctx = buildKneeExtContext20mm();

obj = @(x) objective_KneeExt20mm(x, ctx);

rng default

if isempty(gcp('nocreate'))
    parpool;
end

optsG = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 100);

[xG, fG] = surrogateopt(obj, ctx.lb, ctx.ub, optsG);

optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 300, ...
    'MeshTolerance', 1e-5, ...
    'StepTolerance', 1e-5);

[xBest, fBest] = patternsearch(obj, xG, [], [], [], [], ctx.lb, ctx.ub, [], optsP);

predBest = predictKneeExt20mm(xBest, ctx);

fprintf('\n========== OPTIMIZED DESIGN VALUES ==========\n')

fprintf('\nObjective value:\n')
fprintf('fBest = %.6g\n', fBest)

fprintf('\np1 original, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', ctx.x0(1:3))

fprintf('\np1 optimized, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', predBest.p1)

fprintf('\np2 original, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', ctx.x0(4:6))

fprintf('\np2 optimized, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', predBest.p2)

fprintf('\np1 change, m:\n')
fprintf('[%+.6f, %+.6f, %+.6f]\n', predBest.p1 - ctx.x0(1:3))

fprintf('\np2 change, m:\n')
fprintf('[%+.6f, %+.6f, %+.6f]\n', predBest.p2 - ctx.x0(4:6))

fprintf('\nBPA/tendon lengths:\n')
fprintf('rest   = %.6f m\n', predBest.rest)
fprintf('tendon = %.6f m\n', predBest.tendon)

fprintf('\nKMAX/kmax:\n')
fprintf('KMAX = %.6f\n', predBest.KMAX)
fprintf('kmax = %.6f m\n', predBest.kmax)

fprintf('\nConstraint check:\n')
fprintf('restLmt           = %.6f m\n', predBest.restLmt)
fprintf('extensionDistance = %.6f m\n', predBest.extensionDistance)
fprintf('cRestLength       = %.6f m\n', predBest.cRestLength)

fprintf('=============================================\n')

figure
plot(ctx.phiD, predBest.TorqueZ, ctx.humanAngleD, ctx.humanTorqueAbs)
legend('BPA torque','Human target')
xlabel('Knee angle, deg')
ylabel('Torque magnitude, N m')