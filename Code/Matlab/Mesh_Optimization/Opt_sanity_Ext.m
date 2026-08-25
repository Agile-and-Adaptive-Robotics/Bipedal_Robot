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
fprintf('Inflated BPA clearance diameter = %.3f mm\n', 2*1000*ctx.geo.bpaRadius)
fprintf('Inflated BPA clearance radius = %.3f mm\n', 1000*ctx.geo.bpaRadius)
fprintf('Route rows: %d total = 5 femur + 4 t1\n', ctx.routeRows)
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
fprintf('pEnd = [%.6f %.6f %.6f] m  [route row p9]\n', pred0.pEnd)
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

[~, idxExtension] = max(ctx.phiD);
reserveAtExtension = pred0.KMAX - pred0.strain_f(idxExtension);

fprintf('Reserve contraction at +10 deg = %.6f strain fraction\n', reserveAtExtension)
fprintf('(No equality constraint forces strain_f(+10 deg) = KMAX.)\n')
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

for j = 1:9
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

fprintf('Expected rough activation order after +10 deg baseline: p3 (~-6), p4 (~-45), p6 (~-60), p5 (~-85)\n')

fprintf('\n========== ACTIVE POINTS OVER SOLVER RANGE ==========\n')
prevActive = false(9,1);
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
for j = 1:9
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
fprintf('Femur crossing-collision frames = %d / %d\n', ...
    nnz(routeInfo.femurCrossingCollision), ctx.N)
if any(routeInfo.femurCrossingCollision)
    fprintf('First crossing collision while flexing = %.4f deg\n', ...
        max(ctx.phiD(routeInfo.femurCrossingCollision)))
end
fprintf('Femur overflow frames after p5 = %d / %d\n', ...
    nnz(routeInfo.femurOverflow), ctx.N)
if any(routeInfo.femurOverflow)
    fprintf('First femur overflow while flexing = %.4f deg\n', ...
        max(ctx.phiD(routeInfo.femurOverflow)))
end
fprintf('===============================================\n')

%% Plot full geometry and p1:p9 route
plotAnglesRouteD = [10 -30 -90 -120];
thPlot = linspace(0,2*pi,200).';

figure('Name','9-point extensor route geometry','Color','w')
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for qPlot = 1:numel(plotAnglesRouteD)

    [~,ii] = min(abs(ctx.phiD - plotAnglesRouteD(qPlot)));

    P = routeInfo.raw(:,:,ii);

    % routeInfo.raw keeps p6:p9 in native t1. Convert them t1 -> ICR ->
    % femur for a common-frame geometry plot.
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), P(j,:));
        P(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
    end

    nexttile
    hold on

    % --- femur cylinder clearance
    C = ctx.geo.femurCylCenter;
    R = ctx.geo.femurCylClearRadius;
    plot(C(1)+R*cos(thPlot), C(2)+R*sin(thPlot), '-', 'LineWidth', 1.2)

    % --- femur wall / line clearance
    plot([ctx.geo.femurLineX ctx.geo.femurLineX], ctx.geo.femurLineY, '-', ...
        'LineWidth', 1.2)

    % --- corrected condyle clearance
    Q = ctx.geo.femurOffsetBoundary;
    plot([Q(:,1);Q(1,1)], [Q(:,2);Q(1,2)], '-', 'LineWidth', 1.5)

    % --- t1 lower circle, plotted in femur frame
    L = [ ...
        ctx.geo.tibiaLowerCenter(1)+ctx.geo.tibiaLowerClearRadius*cos(thPlot), ...
        ctx.geo.tibiaLowerCenter(2)+ctx.geo.tibiaLowerClearRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];
    Lf = zeros(size(L));
    for kk = 1:size(L,1)
        qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), L(kk,:));
        Lf(kk,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
    end
    plot(Lf(:,1), Lf(:,2), '-', 'LineWidth', 1.2)

    % --- t1 upper circle, plotted in femur frame
    U = [ ...
        ctx.geo.tibiaUpperCenter(1)+ctx.geo.tibiaUpperClearRadius*cos(thPlot), ...
        ctx.geo.tibiaUpperCenter(2)+ctx.geo.tibiaUpperClearRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];
    Uf = zeros(size(U));
    for kk = 1:size(U,1)
        qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), U(kk,:));
        Uf(kk,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
    end
    plot(Uf(:,1), Uf(:,2), '-', 'LineWidth', 1.2)

    % --- tibia 3-point arc circle, plotted in femur frame
    A = [ ...
        ctx.geo.tibiaArcCenter(1)+ctx.geo.tibiaArcRadius*cos(thPlot), ...
        ctx.geo.tibiaArcCenter(2)+ctx.geo.tibiaArcRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];
    Af = zeros(size(A));
    for kk = 1:size(A,1)
        qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), A(kk,:));
        Af(kk,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
    end
    plot(Af(:,1), Af(:,2), '--', 'LineWidth', 1.0)

    % --- route
    plot(P(:,1), P(:,2), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5)

    for j = 1:9
        if routeInfo.active(j,ii)
            text(P(j,1), P(j,2), sprintf(' p%d',j), ...
                'FontSize', 8, 'Interpreter', 'none')
        end
    end

    axis equal
    grid on
    xlabel('Femur-frame x, m')
    ylabel('Femur-frame y, m')
    title(sprintf('\\phi = %.2f deg',ctx.phiD(ii)))
    legend('Femur cylinder clr', 'Femur line clr', ...
        'Corrected condyle clr', ...
        'Tibia lower clr', 'Tibia upper clr', 'Tibia p6-p8 arc', ...
        'p1:p9 route', 'Location', 'best')
end

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
runSmallOptimizer = false;

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
    fprintf('pEnd = [%.6f %.6f %.6f] m\n', predTest.pEnd)
    fprintf('rest   = %.6f m\n', predTest.rest)
    fprintf('tendon = %.6f m\n', predTest.tendon)
    fprintf('maxPathLength = %.6f m\n', predTest.maxPathLength)
    fprintf('restLmt       = %.6f m\n', predTest.restLmt)
    fprintf('cRestLength   = %.6f m\n', predTest.cRestLength)
    fprintf('==========================================\n')
end
