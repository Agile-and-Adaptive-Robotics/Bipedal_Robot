%% Run this first instead of the full optimizer as a sanity check.

clear
clc
% clear classes
rehash

ctx = buildKneeFlexorContext20mm();

geo = buildGeoExclusion();
idxP2 = 4:6;

J0 = objective_KneeFlexor20mm(ctx.x0, ctx);
pred0 = predictKneeFlexor20mm(ctx.x0, ctx);


disp(J0)
disp(pred0.failReason)

if ~pred0.ok
    disp('predictKneeFlexor20mm failed before creating full pred fields.')
    disp('Available pred0 fields:')
    disp(fieldnames(pred0))
    return
end

disp(pred0.extensionDistance)
disp(pred0.restLmt)
disp(pred0.cRestLength)
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
legend('BPA torque','Human target')
xlabel('Knee angle, deg')
ylabel('Torque magnitude, N m')

%% After the first section runs and the plot looks reasonable, run the following code
obj = @(x) objective_KneeFlexor20mm(x, ctx);
objconstr = @(x) objconstrExclusion(x, obj, geo, ctx, idxP2);
nonlcon = @(x) nonlconExclusion(x, geo, ctx, idxP2);

[c0, ~, collision0] = nonlcon(ctx.x0);
fprintf('p2 exclusion c0       = %.6f m\n', c0(1))
fprintf('BPA-tibia collision c0= %.6f m\n', c0(2))
fprintf('BPA-femur collision c0= %.6f m\n', c0(3))
fprintf('series-length c0      = %.6f m\n', c0(4))
fprintf('collision angle       = %.6f deg\n', collision0.angleD)
fprintf('inflated BPA radius   = %.6f m\n', collision0.bpaRadius)
fprintf('tendon radius         = %.6f m\n', collision0.tendonRadius)
fprintf('tendon length         = %.6f m\n', collision0.tendon)
fprintf('tendon length checked = %.6f m\n', collision0.tendonLengthChecked)
fprintf('current Lm            = %.6f m\n', collision0.currentMuscleLength)
fprintf('Lm + two fittings     = %.6f m\n', collision0.bpaFittingsLengthChecked)
fprintf('minimum tibia clear.  = %.6f m\n', collision0.minClearanceTibia)
fprintf('minimum femur clear.  = %.6f m\n', collision0.minClearanceFemur)
fprintf('required clearance    = %.6f m\n', collision0.requiredClearance)
fprintf('worst tibia region    = %s\n', collision0.worstTibiaRegion)
fprintf('worst tibia component = %s (radius %.6f m)\n', ...
    collision0.worstTibiaComponent, collision0.worstTibiaRadius)
fprintf('worst tibia ctr., t1  = [%.6f %.6f %.6f] m\n', ...
    collision0.worstTibiaCenterT1)
fprintf('worst femur component = %s (radius %.6f m)\n', ...
    collision0.worstFemurComponent, collision0.worstFemurRadius)
fprintf('worst femur ctr.      = [%.6f %.6f %.6f] m\n', ...
    collision0.worstFemurCenterFemur)

optsTest = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'MaxFunctionEvaluations', 30);

[xTest, fTest] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsTest);
