function Export_AppendixRouteGeo()
% Export_AppendixRouteGeo  Dissertation Appendix C route-geometry figure.
%
% Reproduces the tiledlayout "Optimized 9-point extensor route geometry"
% block from Opt_run_Ext.m (lines 493-756) together with its local helpers
% plotStyle/loadColors/styleAxis/styleLegend/makeRouteLegend
% (Opt_run_Ext.m lines 760-914), loading the dated full-workspace result
% mat (the 2026-09-25 design of record) instead of re-running the
% optimizer. Opt_run_Ext.m itself cannot run headless top-to-bottom: the
% optimizer sections execute first and the file carries no export call.
%
% Run headless:
%   matlab -batch "cd('D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization'); Export_AppendixRouteGeo"
%
% Outputs (Figures/94-AppendixC):
%   route_geometry_ext20mm_20260925.pdf  (vector)
%   route_geometry_ext20mm_20260925.png  (200 dpi raster twin)
%   route_geometry_ext20mm_20260925.alt.txt (alt text sidecar)
%
% Figure standard compliance (Ben): Colors.m accessible palette via
% loadColors(); 10 pt text floor enforced by overriding plt.lgdFontsz
% (Opt_run_Ext.m uses 8 pt for this figure's legend and point labels).

%% Paths
thisDir = fileparts(mfilename('fullpath'));   % ...\Code\Matlab\Mesh_Optimization
codeMatlab = fileparts(thisDir);              % ...\Code\Matlab

% Colors.m lives at the Code\Matlab root; RowVecTrans lives in
% Functions\ModernRobotics (the legacy Previous Optimization Code copy is
% deliberately NOT added, so it cannot shadow).
addpath(codeMatlab);
addpath(fullfile(codeMatlab, 'Functions'));
addpath(fullfile(codeMatlab, 'Functions', 'ModernRobotics'));

root = codeMatlab;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', thisDir)
    end
    root = parent;
end

figDir = fullfile(root, 'Documentation', 'Reports and Papers', ...
    'Dissertation', 'Figures', '94-AppendixC');
if ~exist(figDir, 'dir')
    mkdir(figDir)
end
pdfTarget = fullfile(figDir, 'route_geometry_ext20mm_20260925.pdf');
pngTarget = fullfile(figDir, 'route_geometry_ext20mm_20260925.png');
altTarget = fullfile(figDir, 'route_geometry_ext20mm_20260925.alt.txt');

%% Load the dated design-of-record result (full-workspace save)
matFile = fullfile(thisDir, 'Results', 'Vas_Pam_20mm_Result_20260925.mat');
S = load(matFile, 'ctx', 'predBest');
ctx = S.ctx;
predBest = S.predBest;
fprintf('Loaded %s (xBest-driven predBest + ctx, routeInfo %d bytes class %s)\n', ...
    matFile, numel(fieldnames(S.predBest.routeInfo)), class(S.predBest.routeInfo));

%% Style (helpers ported verbatim from Opt_run_Ext.m:760-914)
plt = plotStyle();
% 10 pt floor (Ben's figure standards): Opt_run_Ext.m's 8 pt legend and
% route-point labels are raised to 10 pt for the dissertation export.
plt.lgdFontsz = 10;
c = plt.hexclr;

%% Plot full optimized geometry and p1:p9 route
% Block below is Opt_run_Ext.m lines 493-756 verbatim (figure handle
% captured as figGeo for the export calls; nothing else changed).

routeInfo = predBest.routeInfo;

% Plot full flexion, the solver frame immediately before each unique
% elimination event, and full extension.
transitionIdx = find(any( ...
    routeInfo.active(:,1:end-1) & ~routeInfo.active(:,2:end), 1));
plotIdx = unique([1, transitionIdx, numel(ctx.phiD)], 'stable');
nPoseTiles = numel(plotIdx);

% Pose tiles flow four per row; the legend gets its own fifth column
% spanning all tile rows instead of consuming a pose tile.
nPosesPerRow = 4;
nTileRows = ceil(nPoseTiles/nPosesPerRow);
nTileCols = nPosesPerRow + 1;

if nPoseTiles > nTileRows*nPosesPerRow
    error('Too many route poses for the requested tiled layout.')
end

thPlot = linspace(0,2*pi,200).';

figGeo = figure( ...
    'Name','Optimized 9-point extensor route geometry', ...
    'Color','w', ...
    'Position', [40, 40, 1900, 250+560*nTileRows]);

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

    % Pose tiles flow four per row; the legend reserves the fifth column
    % (spanning all tile rows). Opt_run_Ext.m's nexttile(tGeo,qPlot) feeds
    % row-major linear indices, so with more than four poses the fifth
    % pose lands in the legend column and the legend's spanning nexttile
    % call then DELETES that pose's axes (verified by tmp_probe_nexttile:
    % tile 5 axes is destroyed and replaced). Map each pose to its
    % row/column tile explicitly so all nPoseTiles poses render and the
    % legend call below stays exactly as written in Opt_run_Ext.m.
    rowPlot = floor((qPlot-1)/nPosesPerRow) + 1;
    colPlot = mod(qPlot-1, nPosesPerRow) + 1;
    ax = nexttile(tGeo, (rowPlot-1)*nTileCols + colPlot);
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

% Legend occupies the fifth column, spanning every tile row, so it never
% pushes the poses into an extra row.
axLeg = nexttile(tGeo, nPosesPerRow+1, [nTileRows, 1]);
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

%% Exports
exportgraphics(figGeo, pdfTarget, 'ContentType', 'vector');
exportgraphics(figGeo, pngTarget, 'Resolution', 200);

altText = [ ...
    'Multi-panel route-geometry figure for the optimized 9-point ' newline ...
    'extensor BPA routing of the 20 mm vasti design of record ' newline ...
    '(Vas_Pam_20mm_Result_20260925.mat). Each panel draws one sagittal' newline ...
    'cross-section in the femur frame at a key knee angle: full flexion,' newline ...
    'the solver frame immediately before each contact-elimination event,' newline ...
    'and full extension. In every panel the magenta circled polyline is' newline ...
    'the optimized p1-to-p9 route; clearance envelopes are drawn for the' newline ...
    'femur cylinder (gold), femur straight wall (orange), corrected' newline ...
    'condyle normal-offset boundary (light orange), and the two tibia' newline ...
    'clearance circles (pink and magenta), with purple dashed local bend' newline ...
    'radius lines at p2 and p8. A gold square marks the optimized origin' newline ...
    'p1 and an orange diamond the optimized distal endpoint pEnd. Panel' newline ...
    'titles name the knee angle and which contact points (p3 to p8) are' newline ...
    'removed in the following frame; a shared legend column explains all' newline ...
    'nine entries.'];
fid = fopen(altTarget, 'w');
fprintf(fid, '%s', altText);
fclose(fid);

fprintf('Exported:\n  %s\n  %s\n  %s\n', pdfTarget, pngTarget, altTarget);
assert(isfile(pdfTarget) && dir(pdfTarget).bytes > 0, 'PDF export failed')
assert(isfile(pngTarget) && dir(pngTarget).bytes > 0, 'PNG export failed')
fprintf('Export_AppendixRouteGeo: done (%d pose tiles, %d tile rows).\n', ...
    nPoseTiles, nTileRows)

end

%% Local helpers: verbatim port of Opt_run_Ext.m lines 760-914
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
