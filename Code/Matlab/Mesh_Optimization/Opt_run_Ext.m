clear
clc
rehash

ctx = buildKneeExtContext20mm();

obj = @(x) objective_KneeExt20mm(x, ctx);

nonlcon = @(x) nonlconExt20mm(x, ctx);

objconstr = @(x) objconstrExt20mm(x, obj, nonlcon);

rng default

if isempty(gcp('nocreate'))
    parpool;
end

optsG = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 100);

[xG, fG] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsG);

optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 300, ...
    'MeshTolerance', 1e-5, ...
    'StepTolerance', 1e-5);

% Do not throw away a better initial design before patternsearch.
f0 = obj(ctx.x0);
if f0 <= fG
    xStart = ctx.x0;
else
    xStart = xG;
end

[xBest, fBest] = patternsearch(obj, xStart, [], [], [], [], ctx.lb, ctx.ub, nonlcon, optsP);

predBest = predictKneeExt20mm(xBest, ctx);

fprintf('\n========== OPTIMIZED DESIGN VALUES ==========\n')

fprintf('\nObjective value:\n')
fprintf('fBest = %.6g\n', fBest)

fprintf('\np1 original, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', ctx.x0(1:3))

fprintf('\np1 optimized, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', predBest.p1)

fprintf('\npEnd original, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', ctx.x0(4:6))

fprintf('\npEnd optimized, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', predBest.pEnd)

fprintf('\np1 change, m:\n')
fprintf('[%+.6f, %+.6f, %+.6f]\n', predBest.p1 - ctx.x0(1:3))

fprintf('\npEnd change, m:\n')
fprintf('[%+.6f, %+.6f, %+.6f]\n', predBest.pEnd - ctx.x0(4:6))

fprintf('\nBPA/tendon lengths:\n')
fprintf('rest   = %.6f m\n', predBest.rest)
fprintf('tendon = %.6f m\n', predBest.tendon)

fprintf('\nKMAX/kmax:\n')
fprintf('KMAX = %.6f\n', predBest.KMAX)
fprintf('kmax = %.6f m\n', predBest.kmax)

fprintf('\nConstraint check:\n')
fprintf('maxPathLength0 = %.6f m\n', predBest.maxPathLength0)
fprintf('restLmt        = %.6f m\n', predBest.restLmt)
fprintf('cRestLength    = %.6f m\n', predBest.cRestLength)

fprintf('=============================================\n')

figure
plot(ctx.phiD, predBest.TorqueZ, ctx.humanAngleD, ctx.humanTorqueAbs)
legend('BPA torque','Human target')
xlabel('Knee angle, deg')
ylabel('Torque magnitude, N m')
