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
    'MaxFunctionEvaluations', 1000);

[xG, fG] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsG);

optsP = optimoptions('patternsearch', ...
    'Display', 'iter', ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 30000, ...
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

restcheck = max(predBest.bpa.MuscleLength)-2*predBest.bpa.FittingLength-predBest.bpa.TendonL-predBest.bpa.Xi0
distanceChange = max(predBest.bpa.MuscleLength)-min(predBest.bpa.MuscleLength)
distancePossible = predBest.bpa.RestingL-predBest.bpa.Kmax

%% ================================================================
%  OPTIMIZED-DESIGN PLOTS
%  ================================================================

%% Interpolate selected human target onto solver knee-angle grid
humanBest = interp1( ...
    ctx.humanAngleD, ...
    ctx.humanTorqueAbs, ...
    ctx.phiD, ...
    'pchip', ...
    'extrap');

%% 1. Optimized torque against selected human target
figure('Name','Optimized extensor torque','Color','w')
hold on
plot(ctx.phiD, predBest.TorqueZ, 'LineWidth', 2)
plot(ctx.humanAngleD, ctx.humanTorqueAbs, 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Torque magnitude, N m')
legend('20 mm BPA prediction', char(ctx.targetName), ...
    'Location', 'best')
title('Optimized extensor torque')

%% 2. Torque margin relative to human target
% Positive = BPA exceeds target
% Negative = BPA falls below target
torqueMargin = (predBest.TorqueZ - humanBest)/humanBest;

figure('Name','Extensor torque margin','Color','w')
hold on
plot(ctx.phiD, torqueMargin, 'LineWidth', 2)
yline(0, '--', 'LineWidth', 1.5)
grid on
xlabel('Knee angle, deg')
ylabel('BPA torque - human target, N m')
title('Optimized torque margin')

%% 3. Path lengths
figure('Name','Optimized extensor path lengths','Color','w')
hold on
plot(ctx.phiD, predBest.pathLength, 'LineWidth', 2)
plot(ctx.phiD, predBest.pathLength0, '--', 'LineWidth', 2)
yline(predBest.restLmt, ':', 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Length, m')
legend('Deformed path length', ...
       'Undeformed path length', ...
       'rest+tendon+2fitting+Xi0', ...
       'Location', 'best')
title('Optimized extensor path length')

%% 4. BPA strain
figure('Name','Optimized extensor strain','Color','w')
hold on
plot(ctx.phiD, predBest.strain_f, 'LineWidth', 2)
plot(ctx.phiD, predBest.strain_p, '--', 'LineWidth', 2)
yline(ctx.KMAX, ':', 'LineWidth', 2)
yline(ctx.minStrain, ':', 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Strain')
legend('strain_f includes Xi3', ...
       'strain_p excludes Xi3', ...
       'KMAX', ...
       'minStrain', ...
       'Location', 'best')
title('Optimized extensor strain')

%% 5. BPA force and effective moment arm
Fmag = vecnorm(predBest.bpa.F_p, 2, 2);
rEff = predBest.TorqueZ ./ max(Fmag, eps);

figure('Name','Optimized force and moment arm','Color','w')

yyaxis left
plot(ctx.phiD, Fmag, 'LineWidth', 2)
ylabel('Total BPA force, N')

yyaxis right
plot(ctx.phiD, 1000*rEff, 'LineWidth', 2)
ylabel('Effective Z moment arm, mm')

xlabel('Knee angle, deg')
grid on
title('Optimized BPA force and effective moment arm')

%% 6. Xi3 bending / routing length loss
figure('Name','Optimized Xi3 length loss','Color','w')
plot(ctx.phiD, 1000*predBest.delta_L, 'LineWidth', 2)
grid on
xlabel('Knee angle, deg')
ylabel('Xi3 length loss, mm')
title('Optimized Xi3 length loss')

%% 7. Active routing contacts
routeInfo = predBest.routeInfo;

figure('Name','Optimized routing contact activation','Color','w')
imagesc(ctx.phiD, 1:ctx.routeRows, double(routeInfo.active))
axis xy
xlabel('Knee angle, deg')
ylabel('Route point')
yticks(1:ctx.routeRows)
yticklabels(compose('p%d',1:ctx.routeRows))
title('Active routing contacts')
colorbar

%% 8. Route smoothness / discontinuity diagnostic

phiMid = ctx.phiD(2:end-1);

d2Path = diff(predBest.pathLength0, 2);
d2Loss = diff(predBest.delta_L, 2);
d2Torque = diff(predBest.TorqueZ, 2);

figure('Name','Optimized route smoothness diagnostic','Color','w')

yyaxis left
hold on
plot(phiMid, 1000*d2Path, 'LineWidth', 1.8)
plot(phiMid, 1000*d2Loss, '--', 'LineWidth', 1.8)
ylabel('Second difference in length, mm')

yyaxis right
plot(phiMid, d2Torque, 'LineWidth', 1.8)
ylabel('Second difference in torque, N m')

xlabel('Knee angle, deg')
grid on
legend('\Delta^2 undeformed path', ...
       '\Delta^2 Xi3 loss', ...
       '\Delta^2 torque', ...
       'Location', 'best')
title('Routing / torque smoothness diagnostic')

%% 9. Wrap angles and effective bend radii
if isfield(routeInfo, 'wrap')

    figure('Name','Optimized wrap angles','Color','w')
    hold on

    hWrap = gobjects(numel(routeInfo.wrap.labels),1);

    for cWrap = 1:numel(routeInfo.wrap.labels)
        hWrap(cWrap) = plot( ...
            ctx.phiD, ...
            routeInfo.wrap.angleDeg(cWrap,:), ...
            'LineWidth', 1.6);
    end

    grid on
    xlabel('Knee angle, deg')
    ylabel('Wrap angle, deg')
    legend(hWrap, routeInfo.wrap.labels, 'Location', 'best')
    title('Optimized wrap angles by contact')

    if isfield(routeInfo.wrap, 'radius')

        figure('Name','Optimized bend radii','Color','w')
        hold on

        hRadius = gobjects(numel(routeInfo.wrap.labels),1);

        for cWrap = 1:numel(routeInfo.wrap.labels)
            hRadius(cWrap) = plot( ...
                ctx.phiD, ...
                1000*routeInfo.wrap.radius(cWrap,:), ...
                'LineWidth', 1.6);
        end

        grid on
        xlabel('Knee angle, deg')
        ylabel('Effective bend radius, mm')
        legend(hRadius, routeInfo.wrap.labels, ...
            'Location', 'best')
        title('Optimized effective bend radii')
    end
end

%% Plot full optimized geometry and p1:p9 route

routeInfo = predBest.routeInfo;

plotAnglesRouteD = [9 7 -6 -44 -58 -91 -120];

thPlot = linspace(0,2*pi,200).';

figure( ...
    'Name','Optimized 9-point extensor route geometry', ...
    'Color','w')

tiledlayout( ...
    'flow', ...
    'TileSpacing','compact', ...
    'Padding','compact')

for qPlot = 1:numel(plotAnglesRouteD)

    [~,ii] = min(abs(ctx.phiD - plotAnglesRouteD(qPlot)));

    %% Native route
    Praw = routeInfo.raw(:,:,ii);
    P = Praw;

    % p6:p9 are native t1 coordinates.
    % Convert t1 -> ICR -> femur so everything is drawn in one frame.
    for j = 6:9

        qICR = RowVecTrans( ...
            ctx.T_ICR_t1(:,:,ii), ...
            P(j,:));

        P(j,:) = RowVecTrans( ...
            ctx.T_Pam(:,:,ii), ...
            qICR);
    end

    nexttile
    hold on

    %% Femur cylinder clearance
    C = ctx.geo.femurCylCenter;
    R = ctx.geo.femurCylClearRadius;

    plot( ...
        C(1)+R*cos(thPlot), ...
        C(2)+R*sin(thPlot), ...
        '-', ...
        'LineWidth',1.2)

    %% Femur straight-wall clearance
    plot( ...
        [ctx.geo.femurLineX ctx.geo.femurLineX], ...
        ctx.geo.femurLineY, ...
        '-', ...
        'LineWidth',1.2)

    %% True normal-offset condyle clearance
    Q = ctx.geo.femurOffsetBoundary;

    plot( ...
        [Q(:,1);Q(1,1)], ...
        [Q(:,2);Q(1,2)], ...
        '-', ...
        'LineWidth',1.5)

    %% Lower tibia clearance circle -> femur frame
    L = [ ...
        ctx.geo.tibiaLowerCenter(1) + ...
            ctx.geo.tibiaLowerClearRadius*cos(thPlot), ...
        ctx.geo.tibiaLowerCenter(2) + ...
            ctx.geo.tibiaLowerClearRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];

    Lf = zeros(size(L));

    for kk = 1:size(L,1)

        qICR = RowVecTrans( ...
            ctx.T_ICR_t1(:,:,ii), ...
            L(kk,:));

        Lf(kk,:) = RowVecTrans( ...
            ctx.T_Pam(:,:,ii), ...
            qICR);
    end

    plot( ...
        Lf(:,1), ...
        Lf(:,2), ...
        '-', ...
        'LineWidth',1.2)

    %% Upper tibia clearance circle -> femur frame
    U = [ ...
        ctx.geo.tibiaUpperCenter(1) + ...
            ctx.geo.tibiaUpperClearRadius*cos(thPlot), ...
        ctx.geo.tibiaUpperCenter(2) + ...
            ctx.geo.tibiaUpperClearRadius*sin(thPlot), ...
        zeros(numel(thPlot),1)];

    Uf = zeros(size(U));

    for kk = 1:size(U,1)

        qICR = RowVecTrans( ...
            ctx.T_ICR_t1(:,:,ii), ...
            U(kk,:));

        Uf(kk,:) = RowVecTrans( ...
            ctx.T_Pam(:,:,ii), ...
            qICR);
    end

    plot( ...
        Uf(:,1), ...
        Uf(:,2), ...
        '-', ...
        'LineWidth',1.2)

    %% Exact optimized condyle tangent radii
    radiusLineX = nan(1,5);
    radiusLineY = nan(1,5);

    if isfield(routeInfo,'wrap')

        qStart = routeInfo.wrap.femurCondyleStart(ii,:);
        qEnd   = routeInfo.wrap.femurCondyleEnd(ii,:);

        if all(isfinite(qStart)) && all(isfinite(qEnd))

            fcc = routeInfo.wrap.femurCondyleCenter;

            radiusLineX = [ ...
                fcc(1), qStart(1), ...
                NaN, ...
                fcc(1), qEnd(1)];

            radiusLineY = [ ...
                fcc(2), qStart(2), ...
                NaN, ...
                fcc(2), qEnd(2)];
        end
    end

    plot( ...
        radiusLineX, ...
        radiusLineY, ...
        '--', ...
        'LineWidth',1.2)

    %% Optimized route
    plot( ...
        P(:,1), ...
        P(:,2), ...
        'o-', ...
        'LineWidth',1.8, ...
        'MarkerSize',5)

    %% Label active route points
    for j = 1:9

        if routeInfo.active(j,ii)

            text( ...
                P(j,1), ...
                P(j,2), ...
                sprintf(' p%d',j), ...
                'FontSize',8, ...
                'Interpreter','none')
        end
    end


%% 10. Highlight optimized design endpoints
    plot( ...
        P(1,1), ...
        P(1,2), ...
        's', ...
        'MarkerSize',8, ...
        'LineWidth',1.5)

    plot( ...
        P(9,1), ...
        P(9,2), ...
        'd', ...
        'MarkerSize',8, ...
        'LineWidth',1.5)

    axis equal
    grid on

    xlabel('Femur-frame x, m')
    ylabel('Femur-frame y, m')

    title(sprintf( ...
        'Optimized route, \\phi = %.2f deg', ...
        ctx.phiD(ii)))

    legend( ...
        'Femur cylinder clr', ...
        'Femur line clr', ...
        'Corrected condyle clr', ...
        'Tibia lower clr', ...
        'Tibia upper clr', ...
        'Condyle tangent radii', ...
        'p1:p9 optimized route', ...
        'Optimized p1', ...
        'Optimized pEnd', ...
        'Location','best')
end