clear
clc
rehash

ctx = buildKneeFlexorContext20mm();

geo = buildGeoExclusion();
idxP2 = 4:6;

obj = @(x) objective_KneeFlexor20mm(x, ctx);
objconstr = @(x) objconstrExclusion(x, obj, geo, idxP2);
nonlcon = @(x) nonlconExclusion(x, geo, idxP2);

rng default

if isempty(gcp('nocreate'))
    parpool;
end

optsG = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 800);

[xG, fG] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsG);

optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 300, ...
    'MeshTolerance', 1e-5, ...
    'StepTolerance', 1e-5);

[xBest, fBest] = patternsearch(obj, xG, [], [], [], [], ctx.lb, ctx.ub, nonlcon, optsP);

predBest = predictKneeFlexor20mm(xBest, ctx);

fprintf('\n========== OPTIMIZED DESIGN VALUES ==========\n')

fprintf('\nObjective value:\n')
fprintf('fBest = %.6g\n', fBest)

fprintf('\np1, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', predBest.p1)

fprintf('\np2, m:\n')
fprintf('[%.6f, %.6f, %.6f]\n', predBest.p2)

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