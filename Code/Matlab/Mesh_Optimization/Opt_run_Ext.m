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
[cBest, ~] = nonlconExt20mm(xBest, ctx);
relativeContractionBest = predBest.bpa.Contraction(:)/predBest.KMAX;

%% Run with adjusted seed
% Section commented out if unused.
% Earlier design, reconstructed from rounded printed values.
% This is a trial seed, not an exact recovery of that solution.
xSeed = [0.032, 0.148, 0.001, ...
         0.033, -0.217, 0.001, ...
         0.828, 0.074];

Jseed = objective_KneeExt20mm(xSeed, ctx);
cSeed = nonlconExt20mm(xSeed, ctx);

boundViolation = max([ ...
    ctx.lb(:) - xSeed(:); ...
    xSeed(:) - ctx.ub(:)]);

fprintf('Earlier-design objective = %.6g\n', Jseed);
fprintf('Maximum nonlinear constraint = %.6g\n', max(cSeed));
fprintf('Maximum bound violation = %.6g\n', boundViolation);

% Adjust only coordinates outside the existing bounds.
xStart = min(max(xSeed(:), ctx.lb(:)), ctx.ub(:)).';

% Use the current context and objective.
objRefine = @(x) objective_KneeExt20mm(x, ctx);
conRefine = @(x) nonlconExt20mm(x, ctx);

% Keep the new result separate from your current xBest.
[xRefined, fRefined, exitRefined] = patternsearch( ...
    objRefine, xStart, ...
    [], [], [], [], ctx.lb, ctx.ub, conRefine, optsP);

cRefined = conRefine(xRefined);

fprintf('Refined objective = %.6g\n', fRefined);
fprintf('Maximum nonlinear constraint = %.6g\n', max(cRefined));
fprintf('Maximum bound violation = %.6g\n', max([ ...
    ctx.lb(:) - xRefined(:); ...
    xRefined(:) - ctx.ub(:)]));
fprintf('Exit flag = %d\n', exitRefined);

xBest = xRefined;
fBest = fRefined;

predBest = predictKneeExt20mm(xBest, ctx);
[cBest, ~] = nonlconExt20mm(xBest, ctx);
relativeContractionBest = predBest.bpa.Contraction(:)/predBest.KMAX;

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
fprintf('%-15s | %10.6f\n', 'max constraint', max(cBest))
fprintf('%-15s | %10.6f\n', 'max Contr/KMAX', max(relativeContractionBest))

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

plt = plotStyle();
c = plt.hexclr;
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
colororder(ax, plt.rgbclr)
plot(ax, ctx.phiD, predBest.TorqueZ, ...
    'Color', c{1}, 'LineWidth', plt.lineW)
plot(ax, ctx.humanAngleD, ctx.humanTorque, ...
    'Color', c{7}, 'LineWidth', plt.lineW)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Torque, N m')
lgd = legend(ax, '20 mm BPA prediction', char(ctx.targetName), ...
    'Location', 'best');
title(ax, 'Optimized extensor torque')
styleAxis(ax, plt, xLim)
styleLegend(lgd, plt)

%% 2. Torque margin relative to human target
% Positive = BPA exceeds target
% Negative = BPA falls below target
torqueMargin = (predBest.TorqueZ(:) - humanBest(:))./humanBest;

figure('Name','Extensor torque margin','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.rgbclr)
plot(ax, ctx.phiD, torqueMargin, ...
    'Color', c{2}, 'LineWidth', plt.lineW)
yline(ax, 0, '--', 'Color', c{7}, ...
    'LineWidth', plt.lineW)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Fractional torque margin')
title(ax, 'Optimized torque margin')
styleAxis(ax, plt, xLim)

%% 3. Path lengths
figure('Name','Optimized extensor path lengths','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.rgbclr)
plot(ax, ctx.phiD, predBest.pathLength, ...
    'Color', c{3}, 'LineWidth', plt.lineW)
plot(ax, ctx.phiD, predBest.pathLength0, '--', ...
    'Color', c{7}, 'LineWidth', plt.lineW)
yline(ax, predBest.restLmt, ':', 'Color', c{5}, ...
    'LineWidth', plt.lineW)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Length, m')
lgd = legend(ax, 'Deformed path length', ...
       'Undeformed path length', ...
       'rest+tendon+2fitting+Xi0', ...
       'Location', 'best');
title(ax, 'Optimized extensor path length')
styleAxis(ax, plt, xLim)
styleLegend(lgd, plt)

torqueMargin = (predBest.TorqueZ(:) - humanBest(:)) ./ humanBest(:);
[minMargin, iWorst] = min(torqueMargin);

shortfallFraction = max(0, ...
    abs(humanBest(:)) - abs(predBest.TorqueZ(:))) ./ ctx.torqueScale;

fprintf('Minimum torque margin = %+.3f%% at %.3f deg\n', ...
    100*minMargin, ctx.phiD(iWorst));
fprintf('Weighted worst-shortfall penalty = %.6g\n', ...
    1e5*max(shortfallFraction));
fprintf('Re-evaluated objective = %.6g\n', ...
    objective_KneeExt20mm(xBest, ctx));

%% 4. BPA strain
contractionBest = predBest.bpa.Contraction(:);

figure('Name','Optimized extensor strain','Color','w')
ax = axes;
hold(ax, 'on')
colororder(ax, plt.rgbclr)
plot(ax, ctx.phiD, predBest.strain_f, ...
    'Color', c{4}, 'LineWidth', plt.lineW)
plot(ax, ctx.phiD, predBest.strain_p, '--', ...
    'Color', c{7}, 'LineWidth', plt.lineW)
plot(ax, ctx.phiD, contractionBest, '-.', ...
    'Color', c{2}, 'LineWidth', plt.lineW)
yline(ax, ctx.KMAX, ':', 'Color', c{5}, ...
    'LineWidth', plt.lineW)
yline(ax, 0, ':', 'Color', c{6}, ...
    'LineWidth', plt.lineW)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Strain')
lgd = legend(ax, 'strain_f includes Xi3', ...
       'strain_p excludes Xi3', ...
       'Contraction', ...
       'KMAX', ...
       'Minimum strain = 0', ...
       'Location', 'best');
title(ax, 'Optimized extensor strain')
styleAxis(ax, plt, xLim)
styleLegend(lgd, plt)

%% 5. BPA force and effective moment arm
Fmag = vecnorm(predBest.bpa.F_p, 2, 2);
rEff = predBest.TorqueZ ./ max(Fmag, eps);

figure('Name','Optimized force and moment arm','Color','w')
ax = axes;
colororder(ax, plt.rgbclr)

yyaxis(ax, 'left')
plot(ax, ctx.phiD, Fmag, ...
    'Color', c{7}, 'LineWidth', plt.lineW)
ylabel(ax, 'Total BPA force, N')

yyaxis(ax, 'right')
plot(ax, ctx.phiD, 1000*rEff, ...
    'Color', c{2}, 'LineWidth', plt.lineW)
ylabel(ax, 'Effective Z moment arm, mm')

xlabel(ax, 'Knee angle, deg')
title(ax, 'Optimized BPA force and effective moment arm')
styleAxis(ax, plt, xLim)
ax.YAxis(1).Color = c{7};
ax.YAxis(2).Color = c{2};

%% 6. Xi3 bending / routing length loss
figure('Name','Optimized Xi3 length loss','Color','w')
ax = axes;
colororder(ax, plt.rgbclr)
plot(ax, ctx.phiD, 1000*predBest.delta_L, ...
    'Color', c{5}, 'LineWidth', plt.lineW)
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Xi3 length loss, mm')
title(ax, 'Optimized Xi3 length loss')
styleAxis(ax, plt, xLim)

%% 7. Active routing contacts
routeInfo = predBest.routeInfo;

figure('Name','Optimized routing contact activation','Color','w')
ax = axes;
imagesc(ax, ctx.phiD, 1:ctx.routeRows, double(routeInfo.active))
axis(ax, 'xy')
colormap(ax, [1 1 1; plt.rgbclr(7,:)])
xlabel(ax, 'Knee angle, deg')
ylabel(ax, 'Route point')
yticks(ax, 1:ctx.routeRows)
yticklabels(ax, compose('p%d',1:ctx.routeRows))
title(ax, 'Active routing contacts')
cb = colorbar(ax);
cb.FontName = plt.fontN;
cb.FontSize = plt.lgdFontsz;
cb.FontWeight = 'bold';
styleAxis(ax, plt, xLim)

%% 8. Route smoothness / discontinuity diagnostic

phiMid = ctx.phiD(2:end-1);

d2Path = diff(predBest.pathLength0, 2);
d2Loss = diff(predBest.delta_L, 2);
d2Torque = diff(predBest.TorqueZ, 2);

figure('Name','Optimized route smoothness diagnostic','Color','w')
ax = axes;
colororder(ax, plt.rgbclr)

yyaxis(ax, 'left')
hold(ax, 'on')
plot(ax, phiMid, 1000*d2Path, ...
    'Color', c{3}, 'LineWidth', plt.lineW)
plot(ax, phiMid, 1000*d2Loss, '--', ...
    'Color', c{5}, 'LineWidth', plt.lineW)
ylabel(ax, 'Second difference in length, mm')

yyaxis(ax, 'right')
plot(ax, phiMid, d2Torque, ...
    'Color', c{7}, 'LineWidth', plt.lineW)
ylabel(ax, 'Second difference in torque, N m')

xlabel(ax, 'Knee angle, deg')
lgd = legend(ax, '\Delta^2 undeformed path', ...
       '\Delta^2 Xi3 loss', ...
       '\Delta^2 torque', ...
       'Location', 'best');
title(ax, 'Routing / torque smoothness diagnostic')
styleAxis(ax, plt, xLim)
styleLegend(lgd, plt)
ax.YAxis(1).Color = c{3};
ax.YAxis(2).Color = c{7};

%% 9. Wrap angles and effective bend radii
if isfield(routeInfo, 'wrap')

    figure('Name','Optimized wrap angles','Color','w')
    ax = axes;
    hold(ax, 'on')
    colororder(ax, plt.rgbclr)

    hWrap = gobjects(numel(routeInfo.wrap.labels),1);

    for cWrap = 1:numel(routeInfo.wrap.labels)
        hWrap(cWrap) = plot(ax, ...
            ctx.phiD, ...
            routeInfo.wrap.angleDeg(cWrap,:), ...
            'Color', c{cWrap}, ...
            'LineWidth', plt.lineW);
    end

    addEliminationLines(routeInfo, ctx.phiD, plt)
    xlabel(ax, 'Knee angle, deg')
    ylabel(ax, 'Wrap angle, deg')
    lgd = legend(ax, hWrap, routeInfo.wrap.labels, 'Location', 'best');
    title(ax, 'Optimized wrap angles by contact')
    styleAxis(ax, plt, xLim)
    styleLegend(lgd, plt)

    if isfield(routeInfo.wrap, 'radius')

        figure('Name','Optimized bend radii','Color','w')
        ax = axes;
        hold(ax, 'on')
        colororder(ax, plt.rgbclr)

        hRadius = gobjects(numel(routeInfo.wrap.labels),1);

        for cWrap = 1:numel(routeInfo.wrap.labels)
            hRadius(cWrap) = plot(ax, ...
                ctx.phiD, ...
                1000*routeInfo.wrap.radius(cWrap,:), ...
                'Color', c{cWrap}, ...
                'LineWidth', plt.lineW);
        end

        addEliminationLines(routeInfo, ctx.phiD, plt)
        xlabel(ax, 'Knee angle, deg')
        ylabel(ax, 'Effective bend radius, mm')
        lgd = legend(ax, hRadius, routeInfo.wrap.labels, ...
            'Location', 'best');
        title(ax, 'Optimized effective bend radii')
        styleAxis(ax, plt, xLim)
        styleLegend(lgd, plt)
    end
end

%% Plot full optimized geometry and p1:p9 route

routeInfo = predBest.routeInfo;

% Plot full flexion, the solver frame immediately before each unique
% elimination event, and full extension.
transitionIdx = find(any( ...
    routeInfo.active(:,1:end-1) & ~routeInfo.active(:,2:end), 1));
plotIdx = unique([1, transitionIdx, numel(ctx.phiD)], 'stable');
nPoseTiles = numel(plotIdx);

if nPoseTiles == 9
    nTileRows = 3;
    nTileCols = 4;
else
    nTileRows = 3;
    nTileCols = 3;
end

if nPoseTiles > nTileRows*nTileCols - 1
    error('Too many route poses for the requested tiled layout.')
end

thPlot = linspace(0,2*pi,200).';

figure( ...
    'Name','Optimized 9-point extensor route geometry', ...
    'Color','w')

tGeo = tiledlayout( ...
    nTileRows, nTileCols, ...
    'TileSpacing','compact', ...
    'Padding','compact');
hGeoLegend = gobjects(9,1);

for qPlot = 1:nPoseTiles

    ii = plotIdx(qPlot);

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
    colororder(ax, plt.rgbclr)
    hGeo = gobjects(9,1);

    %% Femur cylinder clearance
    C = ctx.geo.femurCylCenter;
    R = ctx.geo.femurCylClearRadius;

    hGeo(1) = plot(ax, ...
        C(1)+R*cos(thPlot), ...
        C(2)+R*sin(thPlot), ...
        '-', ...
        'Color', c{1}, ...
        'LineWidth', plt.lineW);

    %% Femur straight-wall clearance
    hGeo(2) = plot(ax, ...
        [ctx.geo.femurLineX ctx.geo.femurLineX], ...
        ctx.geo.femurLineY, ...
        '-', ...
        'Color', c{2}, ...
        'LineWidth', plt.lineW);

    %% True normal-offset condyle clearance
    Q = ctx.geo.femurOffsetBoundary;

    hGeo(3) = plot(ax, ...
        [Q(:,1);Q(1,1)], ...
        [Q(:,2);Q(1,2)], ...
        '-', ...
        'Color', c{3}, ...
        'LineWidth', plt.lineW);
    if isfield(ctx.geo, 'femurCondyleClipY')
        plot(ax, [ctx.geo.femurCondyleClipX ctx.geo.femurCondyleClipX], ...
            ctx.geo.femurCondyleClipY, '-', ...
            'Color', hGeo(3).Color, ...
            'LineWidth', plt.lineW, ...
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
        'Color', c{4}, ...
        'LineWidth', plt.lineW);

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
        'Color', c{5}, ...
        'LineWidth', plt.lineW);

    %% Local p2/p8 bend radius lines
    tibiaLowerCenter = [ctx.geo.tibiaLowerCenter, 0];
    tibiaLowerCenter = RowVecTrans( ...
        ctx.T_ICR_t1(:,:,ii), ...
        tibiaLowerCenter);
    tibiaLowerCenter = RowVecTrans( ...
        ctx.T_Pam(:,:,ii), ...
        tibiaLowerCenter);

    radiusLineX = NaN;
    radiusLineY = NaN;

    if routeInfo.active(2,ii)
        radiusLineX = [radiusLineX, ...
            ctx.geo.femurCylCenter(1), P(2,1), NaN];
        radiusLineY = [radiusLineY, ...
            ctx.geo.femurCylCenter(2), P(2,2), NaN];
    end

    if routeInfo.active(8,ii)
        radiusLineX = [radiusLineX, ...
            tibiaLowerCenter(1), P(8,1)];
        radiusLineY = [radiusLineY, ...
            tibiaLowerCenter(2), P(8,2)];
    end

    hGeo(6) = plot(ax, ...
        radiusLineX, ...
        radiusLineY, ...
        '--', ...
        'Color', c{6}, ...
        'LineWidth', plt.lineW);

    %% Optimized route
    hGeo(7) = plot(ax, ...
        P(:,1), ...
        P(:,2), ...
        'o-', ...
        'Color', c{5}, ...
        'LineWidth', plt.lineW, ...
        'MarkerSize', plt.markersz, ...
        'MarkerFaceColor', c{5}, ...
        'MarkerEdgeColor', 'none');

    %% Label active route points
    for j = 1:9

        if routeInfo.active(j,ii)

            text(ax, ...
                P(j,1), ...
                P(j,2), ...
                sprintf(' p%d',j), ...
                'FontName', plt.fontN, ...
                'FontSize', plt.lgdFontsz, ...
                'FontWeight', 'bold', ...
                'Interpreter','none')
        end
    end


%% 10. Highlight optimized design endpoints
            hGeo(8) = scatter(ax, P(1,1), P(1,2), plt.scattersz, ...
                'Marker', 's', ...
                'MarkerFaceColor', c{1}, ...
                'MarkerEdgeColor', 'none');
        
            hGeo(9) = scatter(ax, P(9,1), P(9,2), plt.scattersz, ...
                'Marker', 'd', ...
                'MarkerFaceColor', c{2}, ...
                'MarkerEdgeColor', 'none');

    if qPlot == 1
        hGeoLegend = hGeo;
    end

    axis(ax, 'equal')

    xlabel(ax, 'Femur-frame x, m')
    ylabel(ax, 'Femur-frame y, m')

    if qPlot == 1 || qPlot == nPoseTiles
        tileTitle = sprintf( ...
            '\\theta_k = %.1f^\\circ', ctx.phiD(ii));
    else
        removedNext = find( ...
            routeInfo.active(:,ii) & ~routeInfo.active(:,ii+1));
        removedText = strjoin( ...
            cellstr(compose('-p%d', removedNext)), ', ');
        tileTitle = sprintf( ...
            '%s, \\theta_k = %.1f^\\circ', ...
            removedText, ctx.phiD(ii));
    end

    title(ax, tileTitle, 'Interpreter', 'tex')
    styleAxis(ax, plt)

end

axLeg = nexttile(tGeo,nTileRows*nTileCols);
makeRouteLegend(axLeg, hGeoLegend, { ...
    'Femur cylinder clr', ...
    'Femur line clr', ...
    'Corrected condyle clr', ...
    'Tibia lower clr', ...
    'Tibia upper clr', ...
    'Local bend radii', ...
    'p1:p9 optimized route', ...
    'Optimized p1', ...
    'Optimized pEnd'}, plt);



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


function label = routeStateTitle(routeInfo, ii)

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



function makeRouteLegend(axLeg, hTemplate, labels, plt)

axis(axLeg, 'off')
hold(axLeg, 'on')
hLegend = gobjects(numel(labels),1);

for k = 1:numel(labels)
    hLegend(k) = plot(axLeg, NaN, NaN, 'MarkerEdgeColor', 'none');

    if k <= numel(hTemplate) && isgraphics(hTemplate(k))
        if isprop(hTemplate(k), 'Color')
            hLegend(k).Color = hTemplate(k).Color;
        elseif isprop(hTemplate(k), 'CData')
            hLegend(k).Color = hTemplate(k).CData(1,:);
            hLegend(k).LineStyle = 'none'; % scatter symbols have no line
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
            hLegend(k).MarkerSize = sqrt(hTemplate(k).SizeData(1));
        end

        if isprop(hTemplate(k), 'MarkerFaceColor')
            faceColor = hTemplate(k).MarkerFaceColor;
            if isequal(faceColor, 'flat') || isequal(faceColor, 'auto')
                faceColor = hLegend(k).Color;
            end
            hLegend(k).MarkerFaceColor = faceColor;
        end
        hLegend(k).MarkerEdgeColor = 'none';
    end
end

lgd = legend(axLeg, hLegend, labels, ...
    'Location', 'northwest', 'NumColumns', 1, 'AutoUpdate', 'off');
styleLegend(lgd, plt)

end


function addEliminationLines(routeInfo, phiD, plt)

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
            'Color', plt.hexclr{6}, ...
            'LineWidth', plt.lineW, ...
            'HandleVisibility', 'off');
    end
end

end
