%% Run this before the full optimizer.

clear
clear functions
clc
rehash

% Save the caller's current root defaults, then clear any marker-edge
% overrides while this script runs. Restore the saved values at the end.
oldLineMarkerEdgeColor = get(groot, 'defaultLineMarkerEdgeColor');
oldScatterMarkerEdgeColor = get(groot, 'defaultScatterMarkerEdgeColor');
set(groot, 'defaultLineMarkerEdgeColor', 'remove', ...
    'defaultScatterMarkerEdgeColor', 'remove');

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
geo = ctx.geo;
idxP2 = 4:6;

predOriginal = predictOriginalKneeFlexor20mm(ctx);
pred0 = predictKneeFlexor20mm(ctx.x0, ctx);
J0 = objective_KneeFlexor20mm(ctx.x0, ctx);

if ~predOriginal.ok
    error('Original no-wrap prediction failed: %s', predOriginal.failReason)
end

if ~pred0.ok
    error('x0 prediction failed: %s', pred0.failReason)
end

fprintf('\n========== FLEXOR X3 SANITY ==========\n')
fprintf('J0                         = %.9g\n', J0)
fprintf('Xi0                        = %.9g m\n', ctx.Xi0)
fprintf('Xi1                        = %.9g N/m\n', ctx.Xi1)
fprintf('Xi2                        = %.9g N/m\n', ctx.Xi2)
fprintf('Xi3                        = %.9g\n', ctx.Xi3)
fprintf('BPA count                  = %d\n', ctx.BPAcount)
fprintf('BPA radius source          = %s\n', pred0.bpaRadiusMode)
fprintf('BPA radius range           = %.6f to %.6f m\n', ...
    min(pred0.bpaRadius), max(pred0.bpaRadius))
if pred0.bpaRadiusMode == "bpaR"
    fprintf('radius iterations           = %d\n', pred0.bpaRadiusIteration)
    fprintf('radius update converged     = %d\n', pred0.bpaRadiusConverged)
    fprintf('maximum radius change       = %.9g m\n', pred0.bpaRadiusChange)
end
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
fprintf('maximum |original Tx|      = %.6f N m\n', ...
    max(abs(predOriginal.TorqueX)))
fprintf('maximum |wrapped x0 Tx|    = %.6f N m\n', ...
    max(abs(pred0.TorqueX)))
fprintf('maximum |original Ty|      = %.6f N m\n', ...
    max(abs(predOriginal.TorqueY)))
fprintf('maximum |wrapped x0 Ty|    = %.6f N m\n', ...
    max(abs(pred0.TorqueY)))
fprintf('extension route length     = %.6f m\n', pred0.extensionDistance)
fprintf('restLmt                    = %.6f m\n', pred0.restLmt)
fprintf('cRestLength                = %+.9f m\n', pred0.cRestLength)
fprintf('=======================================\n')

plt = plotStyle();
c = plt.hexclr;

%% Torque
figure('Name','Flexor x0 torque','Color','w')
plot(ctx.phiD, predOriginal.TorqueZ, '--', 'LineWidth', plt.lineW)
hold on
plot(ctx.phiD, pred0.TorqueZ, 'LineWidth', plt.lineW)
plot(ctx.humanAngleD, -ctx.humanTorqueAbs, ':k', 'LineWidth', plt.lineW)
yline(0, ':', 'LineWidth', plt.lineW)
box off
grid off
xlabel('\theta_k, deg')
ylabel('Torque, N m')
lgd = legend('Original BPA', 'Wrapped x0 BPA', 'Human target', ...
    'Location', 'best', 'Box', 'off');
styleAxis(gca, plt)
styleLegend(lgd, plt)

%% Strain: requested three definitions, with the plotted minimum at zero.
figure('Name','Flexor x0 strain','Color','w')
plot(ctx.phiD, pred0.strain_f, '-', 'LineWidth', plt.lineW, ...
    'DisplayName', 'strain_f, includes Xi3')
hold on
plot(ctx.phiD, pred0.strain_p, '--', 'LineWidth', plt.lineW, ...
    'DisplayName', 'strain_p, excludes Xi3')
plot(ctx.phiD, pred0.Contraction, '-.', 'LineWidth', plt.lineW, ...
    'DisplayName', 'Contraction')
yline(0, ':k', 'Minimum strain', 'LineWidth', plt.lineW, ...
    'HandleVisibility', 'off')
yline(ctx.KMAX, ':', 'KMAX', 'LineWidth', plt.lineW, ...
    'HandleVisibility', 'off')
box off
grid off
xlabel('\theta_k, deg')
ylabel('Strain')
lgd = legend('Location', 'best', 'Box', 'off');
styleAxis(gca, plt)
styleLegend(lgd, plt)

%% X-axis torque baseline
figure('Name','Flexor x0 X-axis torque','Color','w')
plot(ctx.phiD, predOriginal.TorqueX, '--', 'LineWidth', plt.lineW)
hold on
plot(ctx.phiD, pred0.TorqueX, 'LineWidth', plt.lineW)
box off
grid off
xlabel('\theta_k, deg')
ylabel('T_x, N m')
title('X-axis Torque')
lgd = legend('Original BPA', 'Wrapped x0 BPA', ...
    'Location', 'best', 'Box', 'off');
styleAxis(gca, plt)
styleLegend(lgd, plt)

%% Y-axis torque baseline
figure('Name','Flexor x0 Y-axis torque','Color','w')
plot(ctx.phiD, predOriginal.TorqueY, '--', 'LineWidth', plt.lineW)
hold on
plot(ctx.phiD, pred0.TorqueY, 'LineWidth', plt.lineW)
box off
grid off
xlabel('\theta_k, deg')
ylabel('T_y, N m')
title('Y-axis Torque')
lgd = legend('Original BPA', 'Wrapped x0 BPA', ...
    'Location', 'best', 'Box', 'off');
styleAxis(gca, plt)
styleLegend(lgd, plt)

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
    hRoute = plot(ax, P(:,1), P(:,2), '-', ...
        'Color', c{5}, 'LineWidth', plt.lineW, ...
        'MarkerEdgeColor', 'none', 'DisplayName', 'BPA route');
    hP1 = scatter(ax, P(1,1), P(1,2), plt.scattersz, ...
        'MarkerFaceColor', c{2}, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'p1');

    % Keep the wrap legend handle even when no wrap is drawn in this pose.
    wrapXY = [NaN, NaN];
    if pred0.routeInfo.active(ii)
        wrapXY = P(2,1:2);
    end
    hWrap = scatter(ax, wrapXY(1), wrapXY(2), plt.scattersz, ...
        'MarkerFaceColor', c{1}, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'pWrap');
    hEnd = scatter(ax, P(3,1), P(3,2), plt.scattersz, ...
        'MarkerFaceColor', c{6}, 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'pEnd');
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

    styleAxis(ax, plt)
    if k == numel(poseIndex)
        lgd = legend(ax, [hRoute, hP1, hWrap, hEnd], ...
            {'BPA route', 'p1', 'pWrap', 'pEnd'}, ...
            'Location', 'best', 'AutoUpdate', 'off');
        styleLegend(lgd, plt)
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

% Restore the marker-edge defaults that were active before this script.
set(groot, 'defaultLineMarkerEdgeColor', oldLineMarkerEdgeColor, ...
    'defaultScatterMarkerEdgeColor', oldScatterMarkerEdgeColor);

function plt = plotStyle()

[plt.hexclr, plt.rgbclr] = loadColors();
plt.lineW = 2;
plt.scattersz = 60;
plt.markersz = 6;
plt.fontN = 'Arial';
plt.axFontsz = 12;
plt.rulerFontsz = 10;
plt.lgdFontsz = 8;
plt.tickL = [0.025, 0.05];

end


function [hexclr, rgbclr] = loadColors()

% Use the existing palette, without a separate hard-coded fallback.
Colors
hexclr = { ...
    c{1}; ... % gold
    c{2}; ... % orange
    c{3}; ... % light orange
    c{4}; ... % pink
    c{5}; ... % magenta
    c{6}; ... % purple (magenta 2 in Colors.m)
    c{7}};    % indigo
rgbclr = d;   % matching RGB rows, in the same color order

end


function styleAxis(ax, plt, xLim)

if nargin >= 3 && ~isempty(xLim) && all(isfinite(xLim))
    xlim(ax, xLim)
end

grid(ax, 'off')
box(ax, 'off')
set(ax, ...
    'FontName', plt.fontN, ...
    'FontSize', plt.rulerFontsz, ...
    'FontWeight', 'bold', ...
    'LineWidth', plt.lineW, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', plt.tickL)

for k = 1:numel(ax.XAxis)
    ax.XAxis(k).LineWidth = plt.lineW;
    ax.XAxis(k).FontSize = plt.rulerFontsz;
    set(ax.XAxis(k).Label, 'FontName', plt.fontN, ...
        'FontSize', plt.axFontsz, 'FontWeight', 'bold')
end

for k = 1:numel(ax.YAxis)
    ax.YAxis(k).LineWidth = plt.lineW;
    ax.YAxis(k).FontSize = plt.rulerFontsz;
    set(ax.YAxis(k).Label, 'FontName', plt.fontN, ...
        'FontSize', plt.axFontsz, 'FontWeight', 'bold')
end

set(ax.Title, 'FontName', plt.fontN, ...
    'FontSize', plt.axFontsz, 'FontWeight', 'bold')

end


function styleLegend(lgd, plt)

if isempty(lgd) || ~isgraphics(lgd)
    return
end

set(lgd, 'FontName', plt.fontN, 'FontSize', plt.lgdFontsz, ...
    'FontWeight', 'bold', 'Box', 'off')

end
