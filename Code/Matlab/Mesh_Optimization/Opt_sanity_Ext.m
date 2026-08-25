%% Run this first instead of the full optimizer as a sanity check.

clearvars
clc
close all
rehash

ctx = buildKneeExtContext20mm();

fprintf('\n========== EXTENSOR CONTEXT SANITY ==========\n')
fprintf('Target: %s\n', ctx.targetName)
fprintf('BPA diameter: %d mm\n', ctx.Dia)
fprintf('BPA count: %d\n', ctx.BPAcount)
fprintf('CrossPoint: %d\n', ctx.CrossPoint)
fprintf('Xi0 = %.6f m\n', ctx.Xi0)
fprintf('Xi1 = %.6g N/m\n', ctx.Xi1)
fprintf('Xi2 = %.6g N/m\n', ctx.Xi2)
fprintf('Xi3 = %.6f\n', ctx.Xi3)
fprintf('Initial tendon = %.6f m\n', ctx.initialTendon0)
fprintf('Initial geometry-dependent tendon maximum = %.6f m\n', ...
    ctx.initialTendonMax)
fprintf('Initial max undeformed Lmt = %.6f m at %.2f deg\n', ...
    ctx.initialMaxLmt0, ctx.initialMaxLmtAngleD)
fprintf('Initial rest = max(Lmt0) - 2*fitting - tendon - Xi0 = %.6f m\n', ...
    ctx.initialRest0)
fprintf('=============================================\n')

%% Evaluate initial design
J0 = objective_KneeExt20mm(ctx.x0, ctx);
pred0 = predictKneeExt20mm(ctx.x0, ctx);

fprintf('\n========== INITIAL DESIGN ==========\n')
fprintf('J0 = %.6g\n', J0)
fprintf('failReason = %s\n', pred0.failReason)

if ~pred0.ok
    error('Initial prediction failed: %s', pred0.failReason)
end

fprintf('\np1 = [%.6f %.6f %.6f] m\n', pred0.p1)
fprintf('p8 = [%.6f %.6f %.6f] m\n', pred0.p8)
fprintf('rest   = %.6f m\n', pred0.rest)
fprintf('tendon = %.6f m\n', pred0.tendon)
fprintf('KMAX   = %.6f\n', pred0.KMAX)
fprintf('kmax   = %.6f m\n', pred0.kmax)

fprintf('\nPath-length check:\n')
fprintf('maxPathLength0 = %.6f m\n', pred0.maxPathLength0)
fprintf('maxPathAngleD0 = %.2f deg\n', pred0.maxPathAngleD0)
fprintf('restLmt        = %.6f m\n', pred0.restLmt)
fprintf('cRestLength    = %.6f m\n', pred0.cRestLength)

if isfield(pred0, 'maxPathLengthDeformed')
    fprintf('maxPathLength deformed diagnostic = %.6f m\n', ...
        pred0.maxPathLengthDeformed)
end

if isfield(pred0, 'maxPathAngleDDeformed')
    fprintf('maxPathAngleD deformed diagnostic = %.2f deg\n', ...
        pred0.maxPathAngleDDeformed)
end

fprintf('\nStrain check:\n')
fprintf('max(strain_f) = %.6f  [includes Xi3, force strain]\n', max(pred0.strain_f))
fprintf('min(strain_f) = %.6f\n', min(pred0.strain_f))
fprintf('max(strain_p) = %.6f  [excludes Xi3, measured-comparison strain]\n', max(pred0.strain_p))
fprintf('min(strain_p) = %.6f\n', min(pred0.strain_p))
fprintf('====================================\n')


%% Routing activation / discontinuity diagnostics
routeInfo = pred0.routeInfo;

fprintf('\n========== ROUTING CONTACT ACTIVATION ==========\n')
fprintf('Sweep increment = %.9f deg\n', routeInfo.sweepStepD)
fprintf('Sweep start     = %.9f deg\n', routeInfo.sweepD(1))
fprintf('Solver max      = %.9f deg\n', max(ctx.phiD))
fprintf('Solver min      = %.9f deg\n\n', min(ctx.phiD))
fprintf('Point   Frame    Added angle, deg       x, m          y, m          z, m\n')
fprintf('-----   -----    ----------------       ----------    ----------    ----------\n')

for j = 1:8
    if isnan(routeInfo.addedAngleD(j))
        fprintf('p%d      %-5s    NOT ACTIVE\n', j, char(routeInfo.pointFrame(j)))
    else
        q = routeInfo.addedPoint(j,:);
        fprintf('p%d      %-5s    %16.9f    % .8f   % .8f   % .8f\n', ...
            j, char(routeInfo.pointFrame(j)), routeInfo.addedAngleD(j), ...
            q(1), q(2), q(3))
    end
end

fprintf('\nActivation order: ')
for k = 1:numel(routeInfo.activationOrder)
    j = routeInfo.activationOrder(k);
    if ~isnan(routeInfo.addedAngleD(j))
        fprintf('p%d ', j)
    end
end
fprintf('\n')

fprintf('\n========== ACTIVE POINTS OVER SOLVER RANGE ==========\n')
prevActive = false(8,1);
for ii = ctx.N:-1:1
    a = routeInfo.active(:,ii);
    if ii == ctx.N || any(a ~= prevActive)
        fprintf('phi = %9.4f deg :', ctx.phiD(ii))
        fprintf(' p%d', find(a))
        fprintf('\n')
    end
    prevActive = a;
end

raw = routeInfo.raw;
loc = pred0.Location;

fprintf('\n========== MAX POINT-TO-POINT JUMPS ==========\n')
for j = 1:8
    Rj = squeeze(raw(j,:,:)).';
    Lj = squeeze(loc(j,:,:)).';

    dRaw = vecnorm(diff(Rj,1,1),2,2);
    dLoc = vecnorm(diff(Lj,1,1),2,2);

    [maxRaw, iRaw] = max(dRaw);
    [maxLoc, iLoc] = max(dLoc);

    fprintf(['p%d: raw = %8.3f mm, %8.3f -> %8.3f deg' ...
             ' | Location = %8.3f mm, %8.3f -> %8.3f deg\n'], ...
        j, ...
        1000*maxRaw, ctx.phiD(iRaw), ctx.phiD(iRaw+1), ...
        1000*maxLoc, ctx.phiD(iLoc), ctx.phiD(iLoc+1))
end

[dPathJump, iPathJump] = max(abs(diff(pred0.pathLength0)));
[dBendJump, iBendJump] = max(abs(diff(pred0.bendMeasure)));
[dLossJump, iLossJump] = max(abs(diff(pred0.delta_L)));

fprintf('\n========== LARGEST CURVE JUMPS ==========\n')
fprintf('undeformed path: %8.3f mm, %8.3f -> %8.3f deg\n', ...
    1000*dPathJump, ctx.phiD(iPathJump), ctx.phiD(iPathJump+1))
fprintf('bendMeasure:     %8.3f mm, %8.3f -> %8.3f deg\n', ...
    1000*dBendJump, ctx.phiD(iBendJump), ctx.phiD(iBendJump+1))
fprintf('X3 delta_L:      %8.3f mm, %8.3f -> %8.3f deg\n', ...
    1000*dLossJump, ctx.phiD(iLossJump), ctx.phiD(iLossJump+1))
fprintf('===============================================\n')

%% Plot torque against target
H = readmatrix('OpenSim_Vasti_Results.txt', ...
    'FileType', 'text', ...
    'NumHeaderLines', 7);

ctx.humanAngleD = H(:,2);
Hang = H(:,2); %Human knee angle data point

Name1 = "Vastus Intermedius"; 
Name2 = "Vastus Lateralis";
Name3 = "Vastus Medialis";

T1 = H(:,3);
T2 = H(:,4);
T3 = H(:,5);

humanAbs = interp1( ...
    ctx.humanAngleD, ...
    ctx.humanTorqueAbs, ...
    ctx.phiD, ...
    'pchip', ...
    'extrap');

figure('Name','Extensor torque sanity','Color','w')
hold on
plot(ctx.phiD, abs(pred0.TorqueZ), 'LineWidth', 2)
plot(Hang, T1, '--', 'LineWidth', 2)
plot(Hang, T2, '-.', 'LineWidth', 2)
plot(Hang, T3, ':', 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Torque magnitude, N m')
legend('20 mm BPA prediction', Name1, Name2, Name3, ...
    'Location', 'best')
title('Initial extensor torque sanity check')

%% Plot path lengths
figure('Name','Extensor path length sanity','Color','w')
hold on
plot(ctx.phiD, pred0.pathLength, 'LineWidth', 2)
plot(ctx.phiD, pred0.pathLength0, '--', 'LineWidth', 2)
yline(pred0.restLmt, ':', 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Length, m')
legend('Deformed path length', 'Undeformed path length', 'rest+tendon+2fitting+Xi0', ...
    'Location', 'best')
title('Extensor path-length sanity check')

%% Plot strain
figure('Name','Extensor strain sanity','Color','w')
hold on
plot(ctx.phiD, pred0.strain_f, 'LineWidth', 2)
plot(ctx.phiD, pred0.strain_p, '--', 'LineWidth', 2)
yline(ctx.KMAX, ':', 'LineWidth', 2)
yline(ctx.minStrain, ':', 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Strain')
legend('strain_f includes Xi3', 'strain_p excludes Xi3', 'KMAX', 'minStrain', ...
    'Location', 'best')
title('Extensor strain sanity check')

Fmag = vecnorm(pred0.bpa.F_p, 2, 2);

rEff = abs(pred0.TorqueZ) ./ Fmag;

figure
yyaxis left
plot(ctx.phiD, Fmag, 'LineWidth', 2)
ylabel('Total BPA force, N')

yyaxis right
plot(ctx.phiD, 1000*rEff, 'LineWidth', 2)
ylabel('Effective Z moment arm, mm')

xlabel('Knee angle, deg')
grid on
figure
plot(ctx.phiD, pred0.delta_L*1000, 'LineWidth', 2)
xlabel('Knee angle, deg')
ylabel('X3 length loss, mm')
grid on

%% Optional tiny optimizer smoke test
runSmallOptimizer = true;

if runSmallOptimizer
    obj = @(x) objective_KneeExt20mm(x, ctx);
    nonlcon = @(x) nonlconExt20mm(x, ctx);
    objconstr = @(x) objconstrExt20mm(x, obj, nonlcon);

    optsTest = optimoptions('surrogateopt', ...
        'Display', 'iter', ...
        'UseParallel', false, ...
        'MaxFunctionEvaluations', 30);

    [xTest, fTest] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsTest);

    predTest = predictKneeExt20mm(xTest, ctx);

    fprintf('\n========== SMALL OPTIMIZER TEST ==========\n')
    fprintf('fTest = %.6g\n', fTest)
    fprintf('p1 = [%.6f %.6f %.6f] m\n', predTest.p1)
    fprintf('p8 = [%.6f %.6f %.6f] m\n', predTest.p8)
    fprintf('rest   = %.6f m\n', predTest.rest)
    fprintf('tendon = %.6f m\n', predTest.tendon)
    fprintf('maxPathLength = %.6f m\n', predTest.maxPathLength)
    fprintf('restLmt       = %.6f m\n', predTest.restLmt)
    fprintf('cRestLength   = %.6f m\n', predTest.cRestLength)
    fprintf('==========================================\n')
end
