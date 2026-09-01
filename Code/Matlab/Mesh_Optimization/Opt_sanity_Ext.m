%% Run this first instead of the full optimizer as a sanity check.

clearvars
clc
close all
clear buildDistalRingLocation20mm buildKneeExtContext20mm
rehash

% Clear old global overrides; markers below specify their edges explicitly.
set(groot, 'defaultLineMarkerEdgeColor', 'remove', ...
    'defaultScatterMarkerEdgeColor', 'remove');

fprintf('Using route builder: %s\n', which('buildDistalRingLocation20mm'))
fprintf('Using context file:  %s\n', which('buildKneeExtContext20mm'))

ctx = buildKneeExtContext20mm();

fprintf('\n========== EXTENSOR CONTEXT SANITY ==========\n')
fprintf('Target: %s\n', ctx.targetName)
fprintf('BPA diameter: %d mm\n', ctx.Dia)
fprintf('BPA count: %d\n', ctx.BPAcount)
fprintf('CrossPoint: %d\n', ctx.CrossPoint)
fprintf('Xi0 = %.3f m\n', ctx.Xi0)
fprintf('Xi1 = %.6g N/m\n', ctx.Xi1)
fprintf('Xi2 = %.6g N/m\n', ctx.Xi2)
fprintf('Xi3 = %.3f\n', ctx.Xi3)
fprintf('Initial tendon = %.3f m\n', ctx.initialTendon0)
fprintf('Initial geometry-dependent tendon maximum = %.3f m\n', ...
    ctx.initialTendonMax)
fprintf('Inflated BPA clearance diameter = %.3f mm\n', 2*1000*ctx.geo.bpaRadius)
fprintf('Inflated BPA clearance radius = %.3f mm\n', 1000*ctx.geo.bpaRadius)
fprintf('Route rows: %d total = 5 femur + 4 t1\n', ctx.routeRows)
fprintf('Initial max undeformed Lmt = %.3f m at %.1f deg\n', ...
    ctx.initialMaxLmt0, ctx.initialMaxLmtAngleD)
fprintf('Initial rest = max(Lmt0) - 2*fitting - tendon - Xi0 = %.3f m\n', ...
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

[c0, ~] = nonlconExt20mm(ctx.x0, ctx);
relativeContraction0 = pred0.bpa.Contraction(:)/pred0.KMAX;

fprintf('\np1   = [% .3f % .3f % .3f] m\n', pred0.p1)
fprintf('pEnd = [% .3f % .3f % .3f] m  [route row p9]\n', pred0.pEnd)
fprintf('rest   = %.3f m\n', pred0.rest)
fprintf('tendon = %.3f m\n', pred0.tendon)
fprintf('KMAX   = %.3f\n', pred0.KMAX)
fprintf('kmax   = %.3f m\n', pred0.kmax)

fprintf('\nPath-length check:\n')
fprintf('maxPathLength0 = %.3f m\n', pred0.maxPathLength0)
fprintf('maxPathAngleD0 = %.1f deg\n', pred0.maxPathAngleD0)
fprintf('restLmt        = %.3f m\n', pred0.restLmt)
fprintf('cRestLength    = %.3f m\n', pred0.cRestLength)
fprintf('max constraint = %+.6g\n', max(c0))

if isfield(pred0, 'maxPathLengthDeformed')
    fprintf('maxPathLength deformed diagnostic = %.3f m\n', ...
        pred0.maxPathLengthDeformed)
end

if isfield(pred0, 'maxPathAngleDDeformed')
    fprintf('maxPathAngleD deformed diagnostic = %.1f deg\n', ...
        pred0.maxPathAngleDDeformed)
end

fprintf('\nStrain check:\n')
fprintf('max(strain_f) = %.3f  [includes Xi3, force strain]\n', max(pred0.strain_f))
fprintf('min(strain_f) = %.3f\n', min(pred0.strain_f))
fprintf('max(strain_p) = %.3f  [excludes Xi3, measured-comparison strain]\n', max(pred0.strain_p))
fprintf('min(strain_p) = %.3f\n', min(pred0.strain_p))
fprintf('max(Contraction/KMAX) = %.6f\n', max(relativeContraction0))

[~, idxExtension] = max(ctx.phiD);
reserveAtExtension = pred0.KMAX - pred0.strain_f(idxExtension);

fprintf('Reserve contraction at +10 deg = %.3f strain fraction\n', reserveAtExtension)
fprintf('(No equality constraint forces strain_f(+10 deg) = KMAX.)\n')
fprintf('====================================\n')


%% Routing elimination / discontinuity diagnostics
routeInfo = pred0.routeInfo;

fprintf('\n========== ROUTING CONTACT ELIMINATION ==========\n')
fprintf('Sweep increment = %.1f deg\n', routeInfo.sweepStepD)
fprintf('Sweep start     = %.1f deg\n', routeInfo.sweepD(1))
fprintf('Solver max      = %.1f deg\n', max(ctx.phiD))
fprintf('Solver min      = %.1f deg\n\n', min(ctx.phiD))
fprintf('%-5s | %-5s | %11s | %14s | %7s | %7s | %7s\n', ...
    'Point', 'Frame', 'Initial deg', 'Eliminated deg', 'x, m', 'y, m', 'z, m')
fprintf('%s\n', ...
    '------+-------+-------------+----------------+---------+---------+---------')

for j = 1:9
    q = routeInfo.addedPoint(j,:);

    if isnan(routeInfo.eliminatedAngleD(j))
        elimText = 'active';
    else
        elimText = sprintf('% .1f', routeInfo.eliminatedAngleD(j));
    end

    fprintf('p%-4d | %-5s | %11.1f | %14s | %7.3f | %7.3f | %7.3f\n', ...
        j, char(routeInfo.pointFrame(j)), routeInfo.addedAngleD(j), ...
        elimText, q(1), q(2), q(3))
end

fprintf('\nElimination order: ')
if isfield(routeInfo, 'eliminationOrder')
    for k = 1:numel(routeInfo.eliminationOrder)
        j = routeInfo.eliminationOrder(k);
        fprintf('p%d ', j)
    end
end
fprintf('\n')

fprintf('\n========== ACTIVE POINTS OVER SOLVER RANGE ==========\n')
fprintf('%8s | %s\n', 'Phi deg', 'Active route rows')
fprintf('%s\n', '---------+------------------')
prevActive = false(9,1);
for ii = 1:ctx.N
    a = routeInfo.active(:,ii);
    if ii == 1 || any(a ~= prevActive)
        fprintf('%8.1f |', ctx.phiD(ii))
        fprintf(' p%d', find(a))
        fprintf('\n')
    end
    prevActive = a;
end

if isfield(routeInfo, 'angleCull')

    sampleAnglesD = [9 7 -6 -44 -58 -91 -120];
    sampleIdx = zeros(size(sampleAnglesD));

    for kk = 1:numel(sampleAnglesD)
        [~, sampleIdx(kk)] = min(abs(ctx.phiD - sampleAnglesD(kk)));
    end

    fprintf('\n========== ABSOLUTE-ANGLE ELIMINATION MARGINS ==========\n')
    fprintf('Femur margins use y-flipped clockwise atan2d; t1 margins use native counterclockwise atan2d.\n')
    fprintf('No principal-angle wrapping is applied to the comparison.\n')
    fprintf('Positive margin means the angle gate allows removal. Collision gates are disabled for this route-state test.\n')
    fprintf('Columns after Point are nearest solver knee-angle samples, deg.\n')
    fprintf('%-8s |', 'Point')

    for kk = 1:numel(sampleIdx)
        fprintf(' %8.1f |', ctx.phiD(sampleIdx(kk)))
    end

    fprintf('\n')

    cullRows = [5 4 3 2 6 7 8];

    for rr = 1:numel(cullRows)
        j = cullRows(rr);
        fprintf('p%-7d |', j)
        fprintf(' %8.3f |', routeInfo.angleCull.marginD(j,sampleIdx))
        fprintf('\n')
    end

    fprintf('====================================================\n')
end

raw = routeInfo.raw;
loc = pred0.Location;

fprintf('\n========== MAX POINT-TO-POINT JUMPS ==========\n')
fprintf('%-5s | %12s | %12s | %10s | %15s | %12s | %10s\n', ...
    'Point', 'Raw jump m', 'Raw from deg', 'Raw to deg', ...
    'Location jump m', 'Loc from deg', 'Loc to deg')
fprintf('%s\n', ...
    '------+--------------+--------------+------------+-----------------+--------------+------------')
for j = 1:9
    Rj = squeeze(raw(j,:,:)).';
    Lj = squeeze(loc(j,:,:)).';

    dRaw = vecnorm(diff(Rj,1,1),2,2);
    dLoc = vecnorm(diff(Lj,1,1),2,2);

    [maxRaw, iRaw] = max(dRaw);
    [maxLoc, iLoc] = max(dLoc);

    fprintf('p%-4d | %12.3f | %12.1f | %10.1f | %15.3f | %12.1f | %10.1f\n', ...
        j, ...
        maxRaw, ctx.phiD(iRaw), ctx.phiD(iRaw+1), ...
        maxLoc, ctx.phiD(iLoc), ctx.phiD(iLoc+1))
end

[dPathJump, iPathJump] = max(abs(diff(pred0.pathLength0)));
[dBendJump, iBendJump] = max(abs(diff(pred0.bendMeasure)));
[dLossJump, iLossJump] = max(abs(diff(pred0.delta_L)));

fprintf('\n========== LARGEST CURVE JUMPS ==========\n')
fprintf('%-18s | %12s | %8s | %8s\n', ...
    'Signal', 'Max jump m', 'From deg', 'To deg')
fprintf('%s\n', ...
    '-------------------+--------------+----------+---------')
fprintf('%-18s | %12.3f | %8.1f | %8.1f\n', 'undeformed path', ...
    dPathJump, ctx.phiD(iPathJump), ctx.phiD(iPathJump+1))
fprintf('%-18s | %12.3f | %8.1f | %8.1f\n', 'bendMeasure', ...
    dBendJump, ctx.phiD(iBendJump), ctx.phiD(iBendJump+1))
fprintf('%-18s | %12.3f | %8.1f | %8.1f\n', 'X3 delta_L', ...
    dLossJump, ctx.phiD(iLossJump), ctx.phiD(iLossJump+1))
fprintf('Interference gates are disabled; route elimination is angle-only in this run.\n')
fprintf('===============================================\n')

if isfield(routeInfo, 'wrap')

    wrap = routeInfo.wrap;

    fprintf('\n========== WRAP ANGLES BY CONTACT ==========\n')
    fprintf('Summary table: max single-step wrap-angle change by contact.\n')
    fprintf('%-24s | %12s | %8s | %8s | %13s\n', ...
        'Contact', 'Max jump deg', 'From deg', 'To deg', 'Max angle deg')
    fprintf('%s\n', ...
        '-------------------------+--------------+----------+----------+---------------')

    for cWrap = 1:numel(wrap.labels)

        aD = wrap.angleDeg(cWrap,:);
        [maxJump, iJump] = max(abs(diff(aD)));

        fprintf('%-24s | %12.1f | %8.1f | %8.1f | %13.1f\n', ...
            wrap.labels{cWrap}, maxJump, ...
            ctx.phiD(iJump), ctx.phiD(iJump+1), max(aD))
    end

    sampleAnglesD = [9 7 -6 -44 -58 -91 -120];
    sampleIdx = zeros(size(sampleAnglesD));

    for kk = 1:numel(sampleAnglesD)
        [~, sampleIdx(kk)] = min(abs(ctx.phiD - sampleAnglesD(kk)));
    end

    fprintf('\nWrap angle samples, deg:\n')
    fprintf('Columns after Contact are nearest solver knee-angle samples, deg.\n')
    fprintf('%-24s |', 'Contact')
    for kk = 1:numel(sampleIdx)
        fprintf(' %8.1f |', ctx.phiD(sampleIdx(kk)))
    end
    fprintf('\n')

    for cWrap = 1:numel(wrap.labels)
        fprintf('%-24s |', wrap.labels{cWrap})
        fprintf(' %8.1f |', wrap.angleDeg(cWrap,sampleIdx))
        fprintf('\n')
    end

    if isfield(wrap, 'radius')
        fprintf('\nEffective bend radius samples, mm:\n')
        fprintf('Columns after Contact are nearest solver knee-angle samples, deg.\n')
        fprintf('%-24s |', 'Contact')
        for kk = 1:numel(sampleIdx)
            fprintf(' %8.1f |', ctx.phiD(sampleIdx(kk)))
        end
        fprintf('\n')

        for cWrap = 1:numel(wrap.labels)
            fprintf('%-24s |', wrap.labels{cWrap})
            fprintf(' %8.3f |', 1000*wrap.radius(cWrap,sampleIdx))
            fprintf('\n')
        end

    end

    fprintf('============================================\n')
end

plt = plotStyle();
c = plt.hexclr;

%% Plot full geometry and p1:p9 route
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

figure('Name','9-point extensor route geometry','Color','w')
tGeo = tiledlayout(nTileRows,nTileCols, ...
    'TileSpacing','compact','Padding','compact');
hGeoLegend = gobjects(7,1);

for qPlot = 1:nPoseTiles

    ii = plotIdx(qPlot);

    Praw = routeInfo.raw(:,:,ii);
    P = Praw;

    % routeInfo.raw keeps p6:p9 in native t1. Convert them t1 -> ICR ->
    % femur for a common-frame geometry plot.
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), P(j,:));
        P(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
    end

    ax = nexttile(tGeo,qPlot);
    hold on
    hGeo = gobjects(7,1);

    % --- femur cylinder clearance
    C = ctx.geo.femurCylCenter;
    R = ctx.geo.femurCylClearRadius;
    hGeo(1) = plot(C(1)+R*cos(thPlot), C(2)+R*sin(thPlot), '-', 'LineWidth', plt.lineW);

    % --- femur wall / line clearance
    hGeo(2) = plot([ctx.geo.femurLineX ctx.geo.femurLineX], ctx.geo.femurLineY, '-', ...
        'LineWidth', plt.lineW);

    % --- corrected condyle clearance
    Q = ctx.geo.femurOffsetBoundary;
    hGeo(3) = plot([Q(:,1);Q(1,1)], [Q(:,2);Q(1,2)], '-', 'LineWidth', plt.lineW);
    if isfield(ctx.geo, 'femurCondyleClipY')
        plot([ctx.geo.femurCondyleClipX ctx.geo.femurCondyleClipX], ...
            ctx.geo.femurCondyleClipY, '-', ...
            'Color', hGeo(3).Color, 'LineWidth', plt.lineW, ...
            'HandleVisibility', 'off')
    end

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
    hGeo(4) = plot(Lf(:,1), Lf(:,2), '-', 'LineWidth', plt.lineW);

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
    hGeo(5) = plot(Uf(:,1), Uf(:,2), '-', 'LineWidth', plt.lineW);

    % --- local p2/p8 bend radius lines
    tibiaLowerCenter = [ctx.geo.tibiaLowerCenter, 0];
    tibiaLowerCenter = RowVecTrans(ctx.T_ICR_t1(:,:,ii), tibiaLowerCenter);
    tibiaLowerCenter = RowVecTrans(ctx.T_Pam(:,:,ii), tibiaLowerCenter);
    radiusLineX = [];
    radiusLineY = [];

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
    if isempty(radiusLineX)
        radiusLineX = NaN;
        radiusLineY = NaN;
    end
    hGeo(6) = plot(radiusLineX, radiusLineY, '--', 'LineWidth', plt.lineW);

    % --- route
    hGeo(7) = plot(P(:,1), P(:,2), 'o-', 'Color', c{5}, ...
        'LineWidth', plt.lineW, 'MarkerSize', plt.markersz, ...
        'MarkerFaceColor', c{5}, 'MarkerEdgeColor', 'none');

    if qPlot == 1
        hGeoLegend = hGeo;
    end

    for j = 1:9
        if routeInfo.active(j,ii)
            text(P(j,1), P(j,2), sprintf(' p%d',j), ...
                'FontName', plt.fontN, 'FontSize', plt.lgdFontsz, ...
                'FontWeight', 'bold', 'Interpreter', 'none')
        end
    end

    axis equal
    grid off
    xlabel('Femur-frame x, m')
    ylabel('Femur-frame y, m')
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

    title(tileTitle, 'Interpreter', 'tex')
    styleAxis(ax, plt)

end

axLeg = nexttile(tGeo,nTileRows*nTileCols);
makeRouteLegend(axLeg, hGeoLegend, { ...
    'Femur cylinder clr', 'Femur line clr', ...
    'Corrected condyle clr', ...
    'Tibia lower clr', 'Tibia upper clr', 'Local bend radii', ...
    'p1:p9 route'}, plt);

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

figure('Name','Extensor torque sanity','Color','w')
hold on
plot(ctx.phiD, pred0.TorqueZ, 'LineWidth', plt.lineW)
plot(Hang, T1, '--', 'LineWidth', plt.lineW)
plot(Hang, T2, '-.', 'LineWidth', plt.lineW)
plot(Hang, T3, ':', 'LineWidth', plt.lineW)
grid off
xlabel('Knee angle, deg')
ylabel('Torque, N m')
lgd = legend('20 mm BPA prediction', Name1, Name2, Name3, ...
    'Location', 'best')
title('Initial extensor torque sanity check')
styleAxis(gca, plt)
styleLegend(lgd, plt)

%% Plot path lengths
figure('Name','Extensor path length sanity','Color','w')
hold on
plot(ctx.phiD, pred0.pathLength, 'LineWidth', plt.lineW)
plot(ctx.phiD, pred0.pathLength0, '--', 'LineWidth', plt.lineW)
yline(pred0.restLmt, ':', 'LineWidth', plt.lineW)
grid off
xlabel('Knee angle, deg')
ylabel('Length, m')
lgd = legend('Deformed path length', 'Undeformed path length', 'rest+tendon+2fitting+Xi0', ...
    'Location', 'best')
title('Extensor path-length sanity check')
styleAxis(gca, plt)
styleLegend(lgd, plt)

if isfield(routeInfo, 'wrap')

    figure('Name','Extensor wrap-angle sanity','Color','w')
    hold on

    hWrap = gobjects(numel(routeInfo.wrap.labels),1);

    for cWrap = 1:numel(routeInfo.wrap.labels)
        hWrap(cWrap) = plot(ctx.phiD, routeInfo.wrap.angleDeg(cWrap,:), ...
            'LineWidth', plt.lineW);
    end

    addEliminationLines(routeInfo, ctx.phiD, plt)
    grid off
    xlabel('Knee angle, deg')
    ylabel('Wrap angle, deg')
    lgd = legend(hWrap, routeInfo.wrap.labels, 'Location', 'best')
    title('Extensor wrap angles by contact')
    styleAxis(gca, plt)
    styleLegend(lgd, plt)

    if isfield(routeInfo.wrap, 'radius')
        figure('Name','Extensor wrap-radius sanity','Color','w')
        hold on

        hRadius = gobjects(numel(routeInfo.wrap.labels),1);

        for cWrap = 1:numel(routeInfo.wrap.labels)
            hRadius(cWrap) = plot(ctx.phiD, ...
                1000*routeInfo.wrap.radius(cWrap,:), 'LineWidth', plt.lineW);
        end

        addEliminationLines(routeInfo, ctx.phiD, plt)
        grid off
        xlabel('Knee angle, deg')
        ylabel('Effective bend radius, mm')
        lgd = legend(hRadius, routeInfo.wrap.labels, 'Location', 'best')
        title('Extensor bend radii by contact')
        styleAxis(gca, plt)
        styleLegend(lgd, plt)
    end
end

%% Plot strain
figure('Name','Extensor strain sanity','Color','w')
hold on
plot(ctx.phiD, pred0.strain_f, 'LineWidth', plt.lineW)
plot(ctx.phiD, pred0.strain_p, '--', 'LineWidth', plt.lineW)
plot(ctx.phiD, pred0.bpa.Contraction(:), '-.', 'LineWidth', plt.lineW)
yline(ctx.KMAX, ':', 'LineWidth', plt.lineW)
yline(0, ':', 'LineWidth', plt.lineW)
grid off
xlabel('Knee angle, deg')
ylabel('Strain')
lgd = legend('strain_f includes Xi3', 'strain_p excludes Xi3', 'Contraction', ...
    'KMAX', 'Minimum strain = 0', ...
    'Location', 'best')
title('Extensor strain sanity check')
styleAxis(gca, plt)
styleLegend(lgd, plt)

Fmag = vecnorm(pred0.bpa.F_p, 2, 2);

rEff = pred0.TorqueZ ./ max(Fmag, eps);

figure('Name','Extensor force and moment arm sanity','Color','w')
yyaxis left
plot(ctx.phiD, Fmag, 'LineWidth', plt.lineW)
ylabel('Total BPA force, N')

yyaxis right
plot(ctx.phiD, 1000*rEff, 'LineWidth', plt.lineW)
ylabel('Effective Z moment arm, mm')

xlabel('Knee angle, deg')
grid off
styleAxis(gca, plt)
figure('Name','Extensor X3 length loss sanity','Color','w')
plot(ctx.phiD, pred0.delta_L*1000, 'LineWidth', plt.lineW)
xlabel('Knee angle, deg')
ylabel('X3 length loss, mm')
grid off
styleAxis(gca, plt)

%% Optional tiny optimizer smoke test
% Disabled until the route-state model is stable enough for optimizer calls
% to diagnose optimizer behavior instead of routing-transition behavior.
runSmallOptimizer = true;

if runSmallOptimizer
    obj = @(x) objective_KneeExt20mm(x, ctx);
    nonlcon = @(x) nonlconExt20mm(x, ctx);
    objconstr = @(x) objconstrExt20mm(x, obj, nonlcon);

    optsTest = optimoptions('surrogateopt', ...
        'Display', 'iter', ...
        'UseParallel', true, ...
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


% Leave subsequent scripts with MATLAB's default marker-edge settings.
set(groot, 'defaultLineMarkerEdgeColor', 'remove', ...
    'defaultScatterMarkerEdgeColor', 'remove');

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
            'LineWidth', plt.lineW, ...
            'HandleVisibility', 'off')
    end
end

end
