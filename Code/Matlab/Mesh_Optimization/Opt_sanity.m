%% Run this before the full optimizer.

clear
clear functions
clc
rehash

set(groot, ...
    'defaultAxesFontName','Arial', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesFontWeight','bold', ...
    'defaultAxesLabelFontSizeMultiplier',1, ...
    'defaultAxesTitleFontSizeMultiplier',1.2, ...
    'defaultAxesLineWidth',2, ...
    'defaultAxesBox','off', ...
    'defaultAxesXMinorTick','on', ...
    'defaultAxesYMinorTick','on', ...
    'defaultAxesTickLength',[0.025 0.05], ...
    'defaultAxesXGrid','off', ...
    'defaultAxesYGrid','off', ...
    'defaultLineLineWidth',2, ...
    'defaultLegendFontName','Arial', ...
    'defaultLegendFontSize',8, ...
    'defaultLegendFontWeight','bold', ...
    'defaultLegendBox','off')

ctx = buildKneeFlexorContext20mm();
geo = buildGeoExclusion();
idxP2 = 4:6;

pred0 = predictKneeFlexor20mm(ctx.x0, ctx);
J0 = objective_KneeFlexor20mm(ctx.x0, ctx);

if ~pred0.ok
    error('x0 prediction failed: %s', pred0.failReason)
end

fprintf('\n========== FLEXOR X3 SANITY ==========\n')
fprintf('J0                         = %.9g\n', J0)
fprintf('Xi0                        = %.9g m\n', ctx.Xi0)
fprintf('Xi1                        = %.9g N/m\n', ctx.Xi1)
fprintf('Xi2                        = %.9g N/m\n', ctx.Xi2)
fprintf('Xi3                        = %.9g\n', ctx.Xi3)
fprintf('p1                         = [%.6f %.6f %.6f] m\n', pred0.p1)
fprintf('wrap, t1                   = [%.6f %.6f %.6f] m\n', pred0.pWrap)
fprintf('pEnd, t1                   = [%.6f %.6f %.6f] m\n', pred0.pEnd)
fprintf('wrap y-z line fraction     = %.6f\n', ...
    pred0.routeInfo.wrapYZFraction)
[~, idxFullExtension] = max(ctx.phiD);
fprintf('wrap active at +10 deg     = %d\n', ...
    pred0.routeInfo.active(idxFullExtension))
fprintf('wrap release found         = %d\n', pred0.routeInfo.releaseFound)
if pred0.routeInfo.releaseFound
    fprintf('first inactive wrap angle   = %+.6f deg\n', ...
        pred0.routeInfo.releaseAngleD)

    releasePosition = find( ...
        pred0.routeInfo.sweepOrder == pred0.routeInfo.releaseIndex, 1);
    idxLastActive = pred0.routeInfo.sweepOrder(releasePosition - 1);
    idxFirstInactive = pred0.routeInfo.releaseIndex;

    fprintf('last active angle           = %+.6f deg\n', ...
        ctx.phiD(idxLastActive))
    fprintf('  aDirect, aWrap            = %+.6f, %+.6f deg\n', ...
        pred0.routeInfo.aDirectD(idxLastActive), ...
        pred0.routeInfo.aWrapEndD(idxLastActive))
    fprintf('first inactive angle        = %+.6f deg\n', ...
        ctx.phiD(idxFirstInactive))
    fprintf('  aDirect, aWrap            = %+.6f, %+.6f deg\n', ...
        pred0.routeInfo.aDirectD(idxFirstInactive), ...
        pred0.routeInfo.aWrapEndD(idxFirstInactive))
else
    fprintf('first inactive wrap angle   = NONE IN MODELED RANGE\n')
end
fprintf('max Contraction/KMAX       = %.6f\n', ...
    max(pred0.relativeContraction))
fprintf('strain_f range             = %.6f to %.6f\n', ...
    min(pred0.strain_f), max(pred0.strain_f))
fprintf('strain_p range             = %.6f to %.6f\n', ...
    min(pred0.strain_p), max(pred0.strain_p))
fprintf('maximum x0 off-axis torque = %.6f N m\n', max(pred0.offAxisTorque))
fprintf('extension route length     = %.6f m\n', pred0.extensionDistance)
fprintf('restLmt                    = %.6f m\n', pred0.restLmt)
fprintf('cRestLength                = %+.9f m\n', pred0.cRestLength)
fprintf('=======================================\n')

%% Torque
figure('Name','Flexor x0 torque','Color','w')
plot(ctx.phiD, pred0.TorqueZ, 'LineWidth', 2)
hold on
plot(ctx.humanAngleD, -ctx.humanTorqueAbs, ':k', 'LineWidth', 3)
yline(0, ':', 'LineWidth', 2)
box off
grid off
xlabel('\theta_k, deg')
ylabel('Torque, N m')
legend('x0 BPA', 'Human target', 'Location', 'best', 'Box', 'off')

%% Strain: requested three definitions, with the plotted minimum at zero.
figure('Name','Flexor x0 strain','Color','w')
plot(ctx.phiD, pred0.strain_f, '-', 'LineWidth', 2, ...
    'DisplayName', 'strain_f, includes Xi3')
hold on
plot(ctx.phiD, pred0.strain_p, '--', 'LineWidth', 2, ...
    'DisplayName', 'strain_p, excludes Xi3')
plot(ctx.phiD, pred0.Contraction, '-.', 'LineWidth', 2, ...
    'DisplayName', 'Contraction')
yline(0, ':k', 'Minimum strain', 'LineWidth', 2, ...
    'HandleVisibility', 'off')
yline(ctx.KMAX, ':', 'KMAX', 'LineWidth', 2, ...
    'HandleVisibility', 'off')
box off
grid off
xlabel('\theta_k, deg')
ylabel('Strain')
legend('Location', 'best', 'Box', 'off')

%% Off-axis baseline
figure('Name','Flexor x0 off-axis torque','Color','w')
plot(ctx.phiD, pred0.offAxisTorque, 'LineWidth', 2)
box off
grid off
xlabel('\theta_k, deg')
ylabel('sqrt(T_x^2 + T_y^2), N m')
title('Pointwise x0 off-axis baseline')

%% Route states: extension, immediately before release, and full flexion.
[~, idxExtension] = max(ctx.phiD);
[~, idxFlexion] = min(ctx.phiD);
idxBefore = idxExtension;

if pred0.routeInfo.releaseFound
    releasePosition = find( ...
        pred0.routeInfo.sweepOrder == pred0.routeInfo.releaseIndex, 1);
    if releasePosition > 1
        idxBefore = pred0.routeInfo.sweepOrder(releasePosition - 1);
    end
end

poseIndex = unique([idxExtension, idxBefore, idxFlexion], 'stable');
figure('Name','Flexor route transition','Color','w')
tiledlayout(1, numel(poseIndex), 'TileSpacing','compact', 'Padding','compact')

for k = 1:numel(poseIndex)
    ii = poseIndex(k);
    P = pred0.Location(:,:,ii);
    for j = ctx.CrossPoint:size(P,1)
        P(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), P(j,:));
    end

    ax = nexttile;
    hold(ax,'on')
    plot(ax, P(:,1), P(:,2), 'o-', 'LineWidth', 2, ...
        'DisplayName', 'BPA route')
    scatter(ax, P(1,1), P(1,2), 50, 'filled', ...
        'DisplayName', 'p1')

    if pred0.routeInfo.active(ii)
        scatter(ax, P(2,1), P(2,2), 50, 'filled', ...
            'DisplayName', 't1 wrap')
    end

    scatter(ax, P(3,1), P(3,2), 50, 'filled', ...
        'DisplayName', 'pEnd')
    axis(ax,'equal')
    box(ax,'off')
    grid(ax,'off')
    xlabel(ax,'x, m')
    ylabel(ax,'y, m')

    if ii == idxBefore && pred0.routeInfo.releaseFound
        title(ax, sprintf('wrap releases next, \\theta_k=%+.2f deg', ...
            ctx.phiD(ii)), 'Interpreter','tex')
    else
        title(ax, sprintf('\\theta_k=%+.2f deg', ctx.phiD(ii)), ...
            'Interpreter','tex')
    end

    if k == numel(poseIndex)
        legend(ax,'Location','best','Box','off', ...
            'FontName','Arial','FontSize',8,'FontWeight','bold')
    end
end

%% Collision and short surrogate smoke test
obj = @(x) objective_KneeFlexor20mm(x, ctx);
objconstr = @(x) objconstrExclusion(x, obj, geo, ctx, idxP2);
nonlcon = @(x) nonlconExclusion(x, geo, ctx, idxP2);

[c0, ~, collision0] = nonlcon(ctx.x0);
fprintf('\nCollision constraints, c <= 0 is feasible:\n')
fprintf('pEnd exclusion             = %+.9f m\n', c0(1))
fprintf('wrapped route vs tibia     = %+.9f m\n', c0(2))
fprintf('wrapped route vs femur     = %+.9f m\n', c0(3))
fprintf('series length              = %+.9f m\n', c0(4))
fprintf('wrap active at collision   = %d\n', collision0.wrapActive)

optsTest = optimoptions('surrogateopt', ...
    'Display', 'iter', ...
    'UseParallel', false, ...
    'MaxFunctionEvaluations', 30);

[xTest, fTest] = surrogateopt(objconstr, ctx.lb, ctx.ub, optsTest); %#ok<NASGU,ASGLU>
