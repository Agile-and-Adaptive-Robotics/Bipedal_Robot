%% Run this first instead of the full optimizer as a sanity check.

clear
clc
clear classes
rehash

ctx = buildKneeFlexorContext20mm();

geo = buildGeoExclusion();
idxP2 = 4:6;

J0 = objective_KneeFlexor20mm(ctx.x0, ctx);
pred0 = predictKneeFlexor20mm(ctx.x0, ctx);


disp(J0)
disp(pred0.failReason)
disp(pred0.extensionDistance)
disp(pred0.restLmt)
disp(pred0.cRestLength)

fprintf('\nInitial design values:\n')
fprintf('p1 = [%.6f %.6f %.6f] m\n', pred0.p1)
fprintf('p2 = [%.6f %.6f %.6f] m\n', pred0.p2)
fprintf('rest   = %.6f m\n', pred0.rest)
fprintf('tendon = %.6f m\n', pred0.tendon)
fprintf('KMAX   = %.6f\n', pred0.KMAX)
fprintf('kmax   = %.6f m\n', pred0.kmax)
fprintf('restLmt = %.6f m\n', pred0.restLmt)

figure
plot(ctx.phiD, pred0.TorqueZ, ctx.humanAngleD, -ctx.humanTorqueAbs)
legend('BPA abs torque','Human abs target')
xlabel('Knee angle, deg')
ylabel('Torque magnitude, N m')

%% After the first section runs and the plot looks reasonable, run the following code
obj = @(x) objective_KneeFlexor20mm(x, ctx);
objconstr = @(x) objconstrExclusion(x, obj, geo, idxP2);
nonlcon = @(x) nonlconExclusion(x, geo, idxP2);

[c0, ~] = nonlcon(ctx.x0);
fprintf('exclusion c0 = %.6f\n', c0)

optsTest = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'MaxFunctionEvaluations', 30);

[xTest, fTest] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsTest);