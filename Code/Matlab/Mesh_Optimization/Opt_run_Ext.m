clear
clc
close all
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
    'MaxFunctionEvaluations', 1200);

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

%% display results
fprintf('\n========== OPTIMIZED DESIGN VALUES ==========\n')

fprintf('\nObjective value:\n')
fprintf('fBest = %.6g\n', fBest)

fprintf('\nEndpoint coordinates, m:\n')
fprintf('%-15s | %9s | %9s | %9s\n', 'Point set', 'x, m', 'y, m', 'z, m')
fprintf('%s\n', '----------------+-----------+-----------+----------')
fprintf('%-15s | %9.3f | %9.3f | %9.3f\n', ...
    'p1 original', ctx.x0(1:3))
fprintf('%-15s | %9.3f | %9.3f | %9.3f\n', ...
    'p1 optimized', predBest.p1)
fprintf('%-15s | %9.3f | %9.3f | %9.3f\n', ...
    'pEnd original', ctx.x0(4:6))
fprintf('%-15s | %9.3f | %9.3f | %9.3f\n', ...
    'pEnd optimized', predBest.pEnd)
fprintf('%-15s | %+9.3f | %+9.3f | %+9.3f\n', ...
    'p1 change', predBest.p1 - ctx.x0(1:3))
fprintf('%-15s | %+9.3f | %+9.3f | %+9.3f\n', ...
    'pEnd change', predBest.pEnd - ctx.x0(4:6))

fprintf('\nBPA/tendon lengths:\n')
fprintf('%-15s | %10s\n', 'Quantity', 'Value, m')
fprintf('%s\n', '----------------+------------')
fprintf('%-15s | %10.3f\n', 'rest', predBest.rest)
fprintf('%-15s | %10.3f\n', 'tendon', predBest.tendon)

fprintf('\nKMAX/kmax:\n')
fprintf('%-15s | %10s\n', 'Quantity', 'Value')
fprintf('%s\n', '----------------+------------')
fprintf('%-15s | %10.3f\n', 'KMAX', predBest.KMAX)
fprintf('%-15s | %10.3f m\n', 'kmax', predBest.kmax)

fprintf('\nConstraint check:\n')
fprintf('%-15s | %10s\n', 'Quantity', 'Value, m')
fprintf('%s\n', '----------------+------------')
fprintf('%-15s | %10.3f\n', 'maxPathLength0', predBest.maxPathLength0)
fprintf('%-15s | %10.3f\n', 'restLmt', predBest.restLmt)
fprintf('%-15s | %10.3f\n', 'cRestLength', predBest.cRestLength)

fprintf('=============================================\n')

restcheck = ...
    max(predBest.bpa.MuscleLength) ...
    - 2*predBest.bpa.FittingLength ...
    - predBest.bpa.TendonL ...
    - predBest.bpa.Xi0;
distanceChange = max(predBest.bpa.MuscleLength) - min(predBest.bpa.MuscleLength);
distancePossible = predBest.bpa.RestingL - predBest.bpa.Kmax;

fprintf('\nBPA travel check:\n')
fprintf('%-17s | %10s\n', 'Quantity', 'Value, m')
fprintf('%s\n', '------------------+------------')
fprintf('%-17s | %10.3f\n', 'restcheck', restcheck)
fprintf('%-17s | %10.3f\n', 'distanceChange', distanceChange)
fprintf('%-17s | %10.3f\n', 'distancePossible', distancePossible)

%% ================================================================
%  OPTIMIZED-DESIGN PLOTS
%  ================================================================

plt = plotStyle20mm();
xLim = [min(ctx.phiD), max(ctx.phiD)];

%% Interpolate selected human target onto solver knee-angle grid
humanBest = interp1( ...
    ctx.humanAngleD, ...
    ctx.humanTorque, ...
    ctx.phiD, ...
    'pchip', ...
    'extrap');

%% 1. Optimized torque against selected human target
figure('Name','Optimized extensor torque','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.colors)
plot(ax, ctx.phiD, predBest.TorqueZ, ...
    'Color', plotColor20mm(plt,1), 'LineWidth', plt.lineWidth)
plot(ax, ctx.humanAngleD, ctx.humanTorque, ...
    'Color', plotColor20mm(plt,7), 'LineWidth', plt.lineWidth)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Torque, N m')
lgd = legend(ax, '20 mm BPA prediction', char(ctx.targetName), ...
    'Location', 'best');
title(ax, 'Optimized extensor torque')
styleAxis20mm(ax, plt, xLim)
styleLegend20mm(lgd, plt)

%% 2. Torque margin relative to human target
% Positive = BPA exceeds target
% Negative = BPA falls below target
torqueMargin = (predBest.TorqueZ - humanBest)/humanBest;

figure('Name','Extensor torque margin','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.colors)
plot(ax, ctx.phiD, torqueMargin, ...
    'Color', plotColor20mm(plt,2), 'LineWidth', plt.lineWidth)
yline(ax, 0, '--', 'Color', plotColor20mm(plt,7), ...
    'LineWidth', plt.lineWidth)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Fractional torque margin')
title(ax, 'Optimized torque margin')
styleAxis20mm(ax, plt, xLim)

%% 3. Path lengths
figure('Name','Optimized extensor path lengths','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.colors)
plot(ax, ctx.phiD, predBest.pathLength, ...
    'Color', plotColor20mm(plt,3), 'LineWidth', plt.lineWidth)
plot(ax, ctx.phiD, predBest.pathLength0, '--', ...
    'Color', plotColor20mm(plt,7), 'LineWidth', plt.lineWidth)
yline(ax, predBest.restLmt, ':', 'Color', plotColor20mm(plt,5), ...
    'LineWidth', plt.lineWidth)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Length, m')
lgd = legend(ax, 'Deformed path length', ...
       'Undeformed path length', ...
       'rest+tendon+2fitting+Xi0', ...
       'Location', 'best');
title(ax, 'Optimized extensor path length')
styleAxis20mm(ax, plt, xLim)
styleLegend20mm(lgd, plt)

%% 4. BPA strain
figure('Name','Optimized extensor strain','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.colors)
plot(ax, ctx.phiD, predBest.strain_f, ...
    'Color', plotColor20mm(plt,4), 'LineWidth', plt.lineWidth)
plot(ax, ctx.phiD, predBest.strain_p, '--', ...
    'Color', plotColor20mm(plt,7), 'LineWidth', plt.lineWidth)
yline(ax, ctx.KMAX, ':', 'Color', plotColor20mm(plt,5), ...
    'LineWidth', plt.lineWidth)
yline(ax, ctx.minStrain, ':', 'Color', plotColor20mm(plt,6), ...
    'LineWidth', plt.lineWidth)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Strain')
lgd = legend(ax, 'strain_f includes Xi3', ...
       'strain_p excludes Xi3', ...
       'KMAX', ...
       'minStrain', ...
       'Location', 'best');
title(ax, 'Optimized extensor strain')
styleAxis20mm(ax, plt, xLim)
styleLegend20mm(lgd, plt)

%% 5. BPA force and effective moment arm
Fmag = vecnorm(predBest.bpa.F_p, 2, 2);
rEff = predBest.TorqueZ ./ max(Fmag, eps);

figure('Name','Optimized force and moment arm','Color','w')
ax = axes;
colororder(ax, plt.colors)

yyaxis(ax, 'left')
plot(ax, ctx.phiD, Fmag, ...
    'Color', plotColor20mm(plt,7), 'LineWidth', plt.lineWidth)
ylabel(ax, 'Total BPA force, N')

yyaxis(ax, 'right')
plot(ax, ctx.phiD, 1000*rEff, ...
    'Color', plotColor20mm(plt,2), 'LineWidth', plt.lineWidth)
ylabel(ax, 'Effective Z moment arm, mm')

xlabel(ax, 'Knee angle, deg')
title(ax, 'Optimized BPA force and effective moment arm')
styleAxis20mm(ax, plt, xLim)
ax.YAxis(1).Color = plotColor20mm(plt,7);
ax.YAxis(2).Color = plotColor20mm(plt,2);

%% 6. Xi3 bending / routing length loss
figure('Name','Optimized Xi3 length loss','Color','w')
ax = axes;
colororder(ax, plt.colors)
plot(ax, ctx.phiD, 1000*predBest.delta_L, ...
    'Color', plotColor20mm(plt,5), 'LineWidth', plt.lineWidth)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Xi3 length loss, mm')
title(ax, 'Optimized Xi3 length loss')
styleAxis20mm(ax, plt, xLim)

%% 7. Active routing contacts
routeInfo = predBest.routeInfo;

figure('Name','Optimized routing contact activation','Color','w')
ax = axes;
imagesc(ax, ctx.phiD, 1:ctx.routeRows, double(routeInfo.active))
axis(ax, 'xy')
colormap(ax, [1 1 1; plotColor20mm(plt,7)])
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Route point')
yticks(ax, 1:ctx.routeRows)
yticklabels(ax, compose('p%d',1:ctx.routeRows))
title(ax, 'Active routing contacts')
cb = colorbar(ax);
cb.FontName = plt.fontName;
cb.FontSize = plt.legendFontSize;
cb.FontWeight = 'bold';
styleAxis20mm(ax, plt, xLim)

%% 8. Route smoothness / discontinuity diagnostic

phiMid = ctx.phiD(2:end-1);

d2Path = diff(predBest.pathLength0, 2);
d2Loss = diff(predBest.delta_L, 2);
d2Torque = diff(predBest.TorqueZ, 2);

figure('Name','Optimized route smoothness diagnostic','Color','w')
ax = axes;
colororder(ax, plt.colors)

yyaxis(ax, 'left')
hold(ax, 'on')
plot(ax, phiMid, 1000*d2Path, ...
    'Color', plotColor20mm(plt,3), 'LineWidth', plt.lineWidth)
plot(ax, phiMid, 1000*d2Loss, '--', ...
    'Color', plotColor20mm(plt,5), 'LineWidth', plt.lineWidth)
ylabel(ax, 'Second difference in length, mm')

yyaxis(ax, 'right')
plot(ax, phiMid, d2Torque, ...
    'Color', plotColor20mm(plt,7), 'LineWidth', plt.lineWidth)
ylabel(ax, 'Second difference in torque, N m')

xlabel(ax, 'Knee angle, deg')
lgd = legend(ax, '\Delta^2 undeformed path', ...
       '\Delta^2 Xi3 loss', ...
       '\Delta^2 torque', ...
       'Location', 'best');
title(ax, 'Routing / torque smoothness diagnostic')
styleAxis20mm(ax, plt, xLim)
styleLegend20mm(lgd, plt)
ax.YAxis(1).Color = plotColor20mm(plt,3);
ax.YAxis(2).Color = plotColor20mm(plt,7);

%% 9. Wrap angles and effective bend radii
if isfield(routeInfo, 'wrap')

    figure('Name','Optimized wrap angles','Color','w')
    ax = axes;
    hold(ax, 'on')
    colororder(ax, plt.colors)

    hWrap = gobjects(numel(routeInfo.wrap.labels),1);

    for cWrap = 1:numel(routeInfo.wrap.labels)
        hWrap(cWrap) = plot(ax, ...
            ctx.phiD, ...
            routeInfo.wrap.angleDeg(cWrap,:), ...
            'Color', plotColor20mm(plt,cWrap), ...
            'LineWidth', plt.lineWidth);
    end

    addEliminationLines20mm(routeInfo, ctx.phiD, plt)
    xlabel(ax, 'Knee angle, deg')
    ylabel(ax, 'Wrap angle, deg')
    lgd = legend(ax, hWrap, routeInfo.wrap.labels, 'Location', 'best');
    title(ax, 'Optimized wrap angles by contact')
    styleAxis20mm(ax, plt, xLim)
    styleLegend20mm(lgd, plt)

    if isfield(routeInfo.wrap, 'radius')

        figure('Name','Optimized bend radii','Color','w')
        ax = axes;
        hold(ax, 'on')
        colororder(ax, plt.colors)

        hRadius = gobjects(numel(routeInfo.wrap.labels),1);

        for cWrap = 1:numel(routeInfo.wrap.labels)
            hRadius(cWrap) = plot(ax, ...
                ctx.phiD, ...
                1000*routeInfo.wrap.radius(cWrap,:), ...
                'Color', plotColor20mm(plt,cWrap), ...
                'LineWidth', plt.lineWidth);
        end

        addEliminationLines20mm(routeInfo, ctx.phiD, plt)
        xlabel(ax, 'Knee angle, deg')
        ylabel(ax, 'Effective bend radius, mm')
        lgd = legend(ax, hRadius, routeInfo.wrap.labels, ...
            'Location', 'best');
        title(ax, 'Optimized effective bend radii')
        styleAxis20mm(ax, plt, xLim)
        styleLegend20mm(lgd, plt)
    end
end

%% Plot full optimized geometry and p1:p9 route

routeInfo = predBest.routeInfo;

plotAnglesRouteD = fliplr([9 7 -6 -44 -58 -91 -120]);

thPlot = linspace(0,2*pi,200).';

figure( ...
    'Name','Optimized 9-point extensor route geometry', ...
    'Color','w')

tGeo = tiledlayout( ...
    3, 3, ...
    'TileSpacing','compact', ...
    'Padding','compact');
hGeoLegend = gobjects(9,1);

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

    ax = nexttile(tGeo,qPlot);
    hold(ax, 'on')
    colororder(ax, plt.colors)
    hGeo = gobjects(9,1);

    %% Femur cylinder clearance
    C = ctx.geo.femurCylCenter;
    R = ctx.geo.femurCylClearRadius;

    hGeo(1) = plot(ax, ...
        C(1)+R*cos(thPlot), ...
        C(2)+R*sin(thPlot), ...
        '-', ...
        'Color', plotColor20mm(plt,1), ...
        'LineWidth', plt.lineWidth);

    %% Femur straight-wall clearance
    hGeo(2) = plot(ax, ...
        [ctx.geo.femurLineX ctx.geo.femurLineX], ...
        ctx.geo.femurLineY, ...
        '-', ...
        'Color', plotColor20mm(plt,2), ...
        'LineWidth', plt.lineWidth);

    %% True normal-offset condyle clearance
    Q = ctx.geo.femurOffsetBoundary;

    hGeo(3) = plot(ax, ...
        [Q(:,1);Q(1,1)], ...
        [Q(:,2);Q(1,2)], ...
        '-', ...
        'Color', plotColor20mm(plt,3), ...
        'LineWidth', plt.lineWidth);
    if isfield(ctx.geo, 'femurCondyleClipY')
        plot(ax, [ctx.geo.femurCondyleClipX ctx.geo.femurCondyleClipX], ...
            ctx.geo.femurCondyleClipY, '-', ...
            'Color', hGeo(3).Color, ...
            'LineWidth', plt.lineWidth, ...
            'HandleVisibility', 'off')
    end

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

    hGeo(4) = plot(ax, ...
        Lf(:,1), ...
        Lf(:,2), ...
        '-', ...
        'Color', plotColor20mm(plt,4), ...
        'LineWidth', plt.lineWidth);

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

    hGeo(5) = plot(ax, ...
        Uf(:,1), ...
        Uf(:,2), ...
        '-', ...
        'Color', plotColor20mm(plt,5), ...
        'LineWidth', plt.lineWidth);

    %% Local p2/p8 bend radius lines
    tibiaLowerCenter = [ctx.geo.tibiaLowerCenter, 0];
    tibiaLowerCenter = RowVecTrans( ...
        ctx.T_ICR_t1(:,:,ii), ...
        tibiaLowerCenter);
    tibiaLowerCenter = RowVecTrans( ...
        ctx.T_Pam(:,:,ii), ...
        tibiaLowerCenter);

    radiusLineX = [ ...
        ctx.geo.femurCylCenter(1), P(2,1), ...
        NaN, ...
        tibiaLowerCenter(1), P(8,1)];

    radiusLineY = [ ...
        ctx.geo.femurCylCenter(2), P(2,2), ...
        NaN, ...
        tibiaLowerCenter(2), P(8,2)];

    hGeo(6) = plot(ax, ...
        radiusLineX, ...
        radiusLineY, ...
        '--', ...
        'Color', plotColor20mm(plt,6), ...
        'LineWidth', plt.lineWidth);

    %% Optimized route
    hGeo(7) = plot(ax, ...
        P(:,1), ...
        P(:,2), ...
        'o-', ...
        'Color', plotColor20mm(plt,7), ...
        'LineWidth', plt.lineWidth, ...
        'MarkerSize', plt.markerSize, ...
        'MarkerFaceColor', plotColor20mm(plt,7), ...
        'MarkerEdgeColor', 'none');

    %% Label active route points
    for j = 1:9

        if routeInfo.active(j,ii)

            text(ax, ...
                P(j,1), ...
                P(j,2), ...
                sprintf(' p%d',j), ...
                'FontName', plt.fontName, ...
                'FontSize', plt.legendFontSize, ...
                'FontWeight', 'bold', ...
                'Interpreter','none')
        end
    end


%% 10. Highlight optimized design endpoints
    hGeo(8) = scatter(ax, ...
        P(1,1), ...
        P(1,2), ...
        plt.scatterSize, ...
        plotColor20mm(plt,1), ...
        's', ...
        'filled', ...
        'MarkerEdgeColor', 'none');

    hGeo(9) = scatter(ax, ...
        P(9,1), ...
        P(9,2), ...
        plt.scatterSize, ...
        plotColor20mm(plt,2), ...
        'd', ...
        'filled', ...
        'MarkerEdgeColor', 'none');

    if qPlot == 1
        hGeoLegend = hGeo;
    end

    axis(ax, 'equal')

    xlabel(ax, 'Femur-frame x, m')
    ylabel(ax, 'Femur-frame y, m')

    title(ax, sprintf( ...
        '%s, phi = %.2f deg', ...
        routeStateTitle20mm(routeInfo, ii), ctx.phiD(ii)), ...
        'Interpreter', 'none')
    styleAxis20mm(ax, plt)

end

axLeg = nexttile(tGeo,numel(plotAnglesRouteD)+1);
makeRouteLegend20mm(axLeg, hGeoLegend, { ...
    'Femur cylinder clr', ...
    'Femur line clr', ...
    'Corrected condyle clr', ...
    'Tibia lower clr', ...
    'Tibia upper clr', ...
    'Local bend radii', ...
    'p1:p9 optimized route', ...
    'Optimized p1', ...
    'Optimized pEnd'}, plt);


function plt = plotStyle20mm()

[hexColors, rgbColors] = loadColors20mm();

plt.hexColors = hexColors;
plt.colors = rgbColors;
plt.lineWidth = 2;
plt.scatterSize = 60;
plt.markerSize = 6;
plt.fontName = availableSansFont20mm();
plt.axisFontSize = 12;
plt.rulerFontSize = 10;
plt.legendFontSize = 8;
plt.tickLength = [0.025, 0.05];

end


function [hexColors, rgbColors] = loadColors20mm()

% Colors.m is the shared accessible palette.  It defines hex colors in c
% and the matching RGB triplets in d; different MATLAB graphics calls accept
% different forms, so keep both available.
colorFile = which('Colors.m');

if isempty(colorFile)
    scriptDir = fileparts(mfilename('fullpath'));
    colorFile = fullfile(fileparts(scriptDir), 'Colors.m');
end

if isfile(colorFile)
    try
        c = cell(0,1);
        d = zeros(0,3);
        run(colorFile)

        if iscell(c) && isnumeric(d) && size(d,2) == 3 && ~isempty(c)
            hexColors = c(:);
            rgbColors = d;
            return
        end
    catch
    end
end

hexColors = { ...
    '#FFD700'; ...
    '#FFB14E'; ...
    '#FA8775'; ...
    '#EA5F94'; ...
    '#CD34B5'; ...
    '#9D02D7'; ...
    '#0000FF'};
rgbColors = [ ...
    1.0000, 0.8431, 0.0000; ...
    1.0000, 0.6941, 0.3059; ...
    0.9804, 0.5294, 0.4588; ...
    0.9176, 0.3725, 0.5804; ...
    0.8039, 0.2039, 0.7098; ...
    0.6157, 0.0078, 0.8431; ...
    0.0000, 0.0000, 1.0000];

end


function fontName = availableSansFont20mm()

preferredFonts = {'Liberation Sans', 'DejaVu Sans', 'Arial', 'Helvetica'};
fontName = preferredFonts{end};

try
    fonts = listfonts;

    for k = 1:numel(preferredFonts)
        if any(strcmpi(fonts, preferredFonts{k}))
            fontName = preferredFonts{k};
            return
        end
    end
catch
end

end


function c = plotColor20mm(plt, k)

idx = mod(k-1, size(plt.colors,1)) + 1;
c = plt.colors(idx,:);

end


function styleAxis20mm(ax, plt, xLim)

if nargin >= 3 && ~isempty(xLim) && all(isfinite(xLim))
    xlim(ax, xLim)
end

grid(ax, 'off')
set(ax, ...
    'FontName', plt.fontName, ...
    'FontSize', plt.axisFontSize, ...
    'FontWeight', 'bold', ...
    'LineWidth', plt.lineWidth, ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'TickLength', plt.tickLength)

for k = 1:numel(ax.XAxis)
    ax.XAxis(k).LineWidth = plt.lineWidth;
    ax.XAxis(k).FontSize = plt.rulerFontSize;
end

for k = 1:numel(ax.YAxis)
    ax.YAxis(k).LineWidth = plt.lineWidth;
    ax.YAxis(k).FontSize = plt.rulerFontSize;
end

ax.XLabel.FontName = plt.fontName;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontName = plt.fontName;
ax.YLabel.FontWeight = 'bold';
ax.Title.FontName = plt.fontName;
ax.Title.FontWeight = 'bold';

end


function styleLegend20mm(lgd, plt)

if isempty(lgd) || ~isgraphics(lgd)
    return
end

lgd.FontName = plt.fontName;
lgd.FontSize = plt.legendFontSize;
lgd.FontWeight = 'bold';

end


function label = routeStateTitle20mm(routeInfo, ii)

optionalRows = 2:8;
inactiveRows = optionalRows(~routeInfo.active(optionalRows,ii));

if isempty(inactiveRows)
    label = 'all';
    return
end

if isfield(routeInfo, 'eliminationOrder')
    orderedRows = routeInfo.eliminationOrder(:).';
    inactiveRows = [ ...
        orderedRows(ismember(orderedRows, inactiveRows)), ...
        setdiff(inactiveRows, orderedRows, 'stable')];
end

label = ['minus ', strjoin(compose('p%d', inactiveRows), ' ')];

end


function makeRouteLegend20mm(axLeg, hTemplate, labels, plt)

axis(axLeg, 'off')
hold(axLeg, 'on')

hLegend = gobjects(numel(labels),1);

for k = 1:numel(labels)
    hLegend(k) = plot(axLeg, NaN, NaN);

    if k <= numel(hTemplate) && isgraphics(hTemplate(k))
        if isprop(hTemplate(k), 'Color')
            hLegend(k).Color = hTemplate(k).Color;
        elseif isprop(hTemplate(k), 'CData')
            hLegend(k).Color = hTemplate(k).CData(1,:);
        end

        if isprop(hTemplate(k), 'LineStyle')
            hLegend(k).LineStyle = hTemplate(k).LineStyle;
        end

        if isprop(hTemplate(k), 'LineWidth')
            hLegend(k).LineWidth = hTemplate(k).LineWidth;
        end

        if isprop(hTemplate(k), 'Marker')
            hLegend(k).Marker = hTemplate(k).Marker;
        end

        if isprop(hTemplate(k), 'MarkerSize')
            hLegend(k).MarkerSize = hTemplate(k).MarkerSize;
        elseif isprop(hTemplate(k), 'SizeData')
            hLegend(k).MarkerSize = sqrt(hTemplate(k).SizeData);
        end
    end
end

lgd = legend(axLeg, hLegend, labels, ...
    'Location', 'northwest', ...
    'NumColumns', 1);
styleLegend20mm(lgd, plt)

end


function addEliminationLines20mm(routeInfo, phiD, plt)

if ~isfield(routeInfo, 'eliminatedAngleD')
    return
end

for j = 1:numel(routeInfo.eliminatedAngleD)
    phi = routeInfo.eliminatedAngleD(j);

    if isfield(routeInfo, 'eliminatedSweepIndex')
        ii = routeInfo.eliminatedSweepIndex(j);

        if isfinite(ii) && ii > 1 && ii <= numel(phiD)
            % Mark the last displayed sample before the repeated-row state.
            phi = phiD(ii-1);
        end
    end

    if isfinite(phi)
        xline(phi, ':', sprintf('p%d',j), ...
            'Color', plotColor20mm(plt,6), ...
            'LineWidth', plt.lineWidth, ...
            'HandleVisibility', 'off');
    end
end

end
