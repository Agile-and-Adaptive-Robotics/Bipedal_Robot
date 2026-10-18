function regen_figs()
% regen_figs  Regenerate the five dissertation figures that FAILED the
% 2026-10-15 print-legibility QC (figure_qc.md), as NEW *_r2.png files in
% Figures/30-torque. Originals are untouched; no original .m or .mat is
% modified.
%
% Ancestors (plotting logic copied from, NOT run):
%   - Vas_Pam_20mm_cam_winner_1..3: Opt_run_ExtCam.m lines 306-352
%     ("Winner diagnostic figures"), data from the dated full-workspace
%     mat minted by save_cam_result_20261014.m
%     (Vas_Pam_20mm_Result_cam_20261014_0416.mat).
%   - route_check_run8 / run8b: the tiledlayout route block of
%     Knee_Extensor_20mm.m lines 770-975 (11-entry legend version) and
%     the headless template Export_AppendixRouteGeo.m. The run-8 design
%     is recomputed from the winner values in tmp/optcam_run8_20261013.log
%     lines 2255-2258; run-8b is the final design of record (predBest in
%     the dated mat).
%
% QC fixes applied (figure_qc.md FAIL verdicts):
%   cam 1: orange-vs-yellow near-overlap -> indigo/magenta/gold with
%          dash AND marker differences; gold given black marker edges.
%   cam 2: gray dashed vs gray dotted -> black solid vs gold dashed,
%          direct right-edge labels, ylim headroom so +5% leaves the axis.
%   cam 3: identical gray dotted 0/KMAX -> black solid vs gold dash-dot,
%          direct labels, ylim headroom.
%   Fonts: cam figures built at 7.5 in native width (0.83x at the 6.2 in
%          print assumption) with 13-14 pt text -> >= 10.8 pt effective.
%   route: 8 panels at 0.38x print scale -> 4 informative panels (full
%          flexion overview + three shared-window distal zooms) on a
%          6.5x9.2 in canvas at ~0.95x print scale, 11-12.5 pt text;
%          route families separated by hue + dash + marker shape.

root = 'D:/Github/Bipedal_Robot';
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

moDir  = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization');
figDir = fullfile(root, 'Documentation', 'Reports and Papers', ...
    'Dissertation', 'Figures', '30-torque');

matFile = fullfile(moDir, 'Results', 'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S = load(matFile, 'ctx', 'predBest', 'marginBest', 'requiredBest', ...
    'routeInfoBase', 'xFixed');
ctx  = S.ctx;
phiD = ctx.phiD(:);
xLim = [min(phiD), max(phiD)];

%% Palette: Tol accessible set from Code/Matlab/Colors.m (c hex / d rgb)
Colors;                                  % defines c, d in this workspace
COL.cam      = d(7, :);   % indigo       -- cam design family
COL.black    = [0 0 0];   % reference/context
COL.gold     = d(1, :);   % gold         -- thresholds / KMAX / p1
COL.orange   = d(2, :);   % orange       -- femur line / pEnd
COL.lorange  = d(3, :);   % light orange -- condyle boundary
COL.pink     = d(4, :);   % pink         -- tibia lower circle
COL.mag      = d(5, :);   % magenta      -- human target / tibia upper
COL.purple   = d(6, :);   % purple       -- bend radii

%% ===================== Cam winner figures 1-3 =========================
% 7.5 x 5.0 in native; at the QC's 6.2 in print width scale = 0.827, so
% 13 pt ticks/legend -> 10.8 pt effective, 14 pt labels -> 11.6 pt.
FS.tick = 13; FS.lab = 14; FS.leg = 13; FS.tit = 14; FS.txt = 13;

mStride = 8;   % marker stride on the dense phiD curves (100 poses)
xR = xLim(2) - 0.025*diff(xLim);   % right-edge anchor for direct labels

% ---- Figure 1: torque vs human target + required (+5%) ----------------
f1 = figure('Color','w','Units','inches','Position',[0 0 7.5 5.0]);
ax1 = axes(f1); hold(ax1, 'on')
plot(ax1, phiD, S.predBest.TorqueZ, '-', 'Color', COL.cam, 'LineWidth', 3);
plot(ax1, ctx.humanAngleD, ctx.humanTorque, '-', 'Color', COL.mag, ...
    'LineWidth', 2.2, 'Marker', 'o', 'MarkerSize', 5.5, ...
    'MarkerFaceColor', COL.mag, 'MarkerEdgeColor', 'none');
plot(ax1, phiD, S.requiredBest, '--', 'Color', COL.gold, 'LineWidth', 3, ...
    'Marker', 'v', 'MarkerSize', 6, 'MarkerIndices', 1:mStride:numel(phiD), ...
    'MarkerFaceColor', COL.gold, 'MarkerEdgeColor', COL.black);
xlim(ax1, xLim); grid(ax1, 'off'); box(ax1, 'off')
xlabel(ax1, 'Knee angle, deg'); ylabel(ax1, 'Torque, N m')
lg1 = legend(ax1, {'Cam BPA prediction', char(ctx.targetName), ...
    'Required (+5%)'}, 'Location', 'best');
title(ax1, 'Cam-optimized extensor torque')
styleAx(ax1, FS); styleLeg(lg1, FS)

% ---- Figure 2: torque margin with separated reference lines -----------
f2 = figure('Color','w','Units','inches','Position',[0 0 7.5 5.0]);
ax2 = axes(f2); hold(ax2, 'on')
hMargin = plot(ax2, phiD, 100*S.marginBest, '-', 'Color', COL.cam, ...
    'LineWidth', 3);
hZero = yline(ax2, 0, '-', 'Color', COL.black, 'LineWidth', 1.8);
hReq  = yline(ax2, 100*ctx.requiredTorqueMargin, '--', 'Color', COL.gold, ...
    'LineWidth', 3);
xlim(ax2, xLim); grid(ax2, 'off'); box(ax2, 'off')
% Margin vs the plain target blows up to ~265% in deep flexion where the
% human target is small; a full-range axis buries the 0 and +5%
% reference lines together at the bottom (the exact QC defect). Cap the
% axis at 120% so the informative band -- min margin ~+38% vs the +5%
% requirement -- fills the panel; the deep-flexion blowup runs off the
% top edge (data unchanged, presentation-only clip).
ylim(ax2, [-10, 120])
xlabel(ax2, 'Knee angle, deg'); ylabel(ax2, 'Torque margin, %')
lg2 = legend(ax2, [hMargin, hZero, hReq], ...
    {'vs plain target', 'zero', 'required +5%'}, 'Location', 'best');
title(ax2, 'Cam torque margin relative to human')
% Direct gold label on the +5% line so it stays identified even in
% grayscale at small size (the black 0 line is labeled by the legend).
text(ax2, xR, 100*ctx.requiredTorqueMargin + 2.5, '+5% required', ...
    'Color', COL.gold, 'FontWeight', 'bold', 'FontSize', FS.txt, ...
    'HorizontalAlignment', 'right', 'FontName', 'Arial');
styleAx(ax2, FS); styleLeg(lg2, FS)

% ---- Figure 3: strain over the RoM ------------------------------------
f3 = figure('Color','w','Units','inches','Position',[0 0 7.5 5.0]);
ax3 = axes(f3); hold(ax3, 'on')
plot(ax3, phiD, S.predBest.strain_f, '-', 'Color', COL.cam, ...
    'LineWidth', 3, 'Marker', 'o', 'MarkerSize', 5, ...
    'MarkerIndices', 1:mStride:numel(phiD), 'MarkerFaceColor', COL.cam, ...
    'MarkerEdgeColor', 'none');
plot(ax3, phiD, S.predBest.strain_p, '--', 'Color', COL.mag, ...
    'LineWidth', 3, 'Marker', '^', 'MarkerSize', 5.5, ...
    'MarkerIndices', 1:mStride:numel(phiD), 'MarkerFaceColor', COL.mag, ...
    'MarkerEdgeColor', 'none');
yline(ax3, 0, '-', 'Color', COL.black, 'LineWidth', 2);
yline(ax3, ctx.KMAX, '-.', 'Color', COL.gold, 'LineWidth', 3);
xlim(ax3, xLim); grid(ax3, 'off'); box(ax3, 'off')
ylim(ax3, [-0.02, 1.12*ctx.KMAX])
xlabel(ax3, 'Knee angle, deg'); ylabel(ax3, 'Strain')
lg3 = legend(ax3, {'strain_f (incl Xi3)', 'strain_p', '0 (no stretch)', ...
    'KMAX'}, 'Location', 'best');
title(ax3, 'Cam BPA strain over the RoM')
text(ax3, xR, 0.004, '0 (no stretch)', 'Color', COL.black, 'FontWeight', ...
    'bold', 'FontSize', FS.txt, 'HorizontalAlignment', 'right', ...
    'FontName', 'Arial');
text(ax3, xR, ctx.KMAX - 0.009, 'KMAX', 'Color', COL.gold, 'FontWeight', ...
    'bold', 'FontSize', FS.txt, 'HorizontalAlignment', 'right', ...
    'FontName', 'Arial');
styleAx(ax3, FS); styleLeg(lg3, FS)

figs = {f1, f2, f3};
for k = 1:3
    outK = fullfile(figDir, sprintf('Vas_Pam_20mm_cam_winner_%d_r2.png', k));
    exportgraphics(figs{k}, outK, 'Resolution', 200);
    fprintf('wrote %s\n', outK);
end

%% ===================== Route-check figures ============================
% 6.5 x 9.2 in native; at the QC's 6.2 in print width scale = 0.954, so
% 11 pt ticks -> 10.5 pt effective, 11.5 pt labels -> 11 pt.
RF.tick = 11; RF.lab = 11.5; RF.tit = 12.5; RF.leg = 11; RF.plab = 11;
RF.sgt = 14;

% run-8 winner geometry verbatim from tmp/optcam_run8_20261013.log
% lines 2255-2258 (the paste block Opt_run_ExtCam printed for run 8).
xCamRun8 = [0.116323, 0.015040, ...    % p6Adj
            0.131323, -0.018974, ...   % p7Adj
            0.066359, -0.018981, ...   % p8Adj
            0.090440];                 % tendonAdj
[xFull8, ctxCam8] = camToExtDesign20mm(xCamRun8, ctx, S.xFixed);
pred8 = predictKneeExt20mm(xFull8, ctxCam8);

% Shared distal-zoom window over poses {1, p3-release, end} for BOTH
% designs, sized to the ACTIVE route points (rows 3:9) plus the cam CAD
% candidates. Eliminated rows repeat neighbor positions far from the
% contact cluster (e.g. cam p3:p5 collapse to the +0.14 m region at
% extension), which would re-inflate the window to the full leg span --
% the exact blob defect the QC failed. Active-only keeps the window on
% the p3..p9 contact zone, shared by both figures for comparability.
iOver = 1;
iEnd  = numel(phiD);
iRel3 = transPoseFor(S.predBest.routeInfo, 3);   % p3 release frame
winPoses = [iOver, iRel3, iEnd];
assert(iRel3 == transPoseFor(pred8.routeInfo, 3), ...
    'p3 release frame differs between run8 and run8b')
pts = [ ...
    activeFramePts(S.routeInfoBase, 3:9, winPoses, ctx); ...
    activeFramePts(S.predBest.routeInfo, 3:9, winPoses, ctx); ...
    activeFramePts(pred8.routeInfo, 3:9, winPoses, ctx); ...
    candFramePts(S.predBest.routeInfo, 6:8, winPoses, ctx); ...
    candFramePts(pred8.routeInfo, 6:8, winPoses, ctx)];
zoomWin = squareWindow(pts, 0.18);
fprintf('shared distal zoom window: x [%.3f %.3f] y [%.3f %.3f] m\n', ...
    zoomWin(1), zoomWin(2), zoomWin(3), zoomWin(4));

make_route_check_figure(S.predBest.routeInfo, ...
    'Final cam design of record (run 8 + p8 clip edge, +2.05 mm tendon)', ...
    fullfile(figDir, 'route_check_run8b_r2.png'), ...
    ctx, phiD, S.routeInfoBase, zoomWin, COL, RF, iRel3);

make_route_check_figure(pred8.routeInfo, ...
    'Run-8 cam winner (pre-adjustment p8, tendon 0.0904 m)', ...
    fullfile(figDir, 'route_check_run8_r2.png'), ...
    ctx, phiD, S.routeInfoBase, zoomWin, COL, RF, iRel3);

fprintf('regen_figs: done\n');
end

%% ======================== local functions =============================

function styleAx(ax, FS)
% Arial everywhere, no italics, QC font floor (see header).
set(ax, 'FontName', 'Arial', 'FontSize', FS.tick, 'FontAngle', 'normal', ...
    'FontWeight', 'normal', 'LineWidth', 1.2)
set(ax.XAxis, 'FontName', 'Arial', 'FontSize', FS.tick)
set(ax.YAxis, 'FontName', 'Arial', 'FontSize', FS.tick)
xlabel(ax, ax.XLabel.String, 'FontName', 'Arial', 'FontSize', FS.lab, ...
    'FontAngle', 'normal', 'FontWeight', 'normal')
ylabel(ax, ax.YLabel.String, 'FontName', 'Arial', 'FontSize', FS.lab, ...
    'FontAngle', 'normal', 'FontWeight', 'normal')
title(ax, ax.Title.String, 'FontName', 'Arial', 'FontSize', FS.tit, ...
    'FontAngle', 'normal', 'FontWeight', 'normal', ...
    'Interpreter', ax.Title.Interpreter)
end

function styleLeg(lgd, FS)
if isempty(lgd) || ~isgraphics(lgd)
    return
end
set(lgd, 'FontName', 'Arial', 'FontSize', FS.leg, 'FontAngle', 'normal', ...
    'FontWeight', 'normal', 'Box', 'off')
end

function ipose = transPoseFor(ri, pointRow)
% Pose index of the frame immediately before `pointRow` is eliminated
% (NaN if never eliminated).
tr = find(any(ri.active(:, 1:end-1) & ~ri.active(:, 2:end), 1));
ipose = NaN;
for ii = tr(:).'
    if any(find(ri.active(:, ii) & ~ri.active(:, ii+1)) == pointRow)
        ipose = ii;
    end
end
end

function P = activeFramePts(ri, rows, poses, ctx)
% Like femurFramePts but keeps only rows ACTIVE at each pose, so
% eliminated repeats (which park far from the contact cluster) cannot
% inflate the zoom window.
P = zeros(0, 2);
for ii = poses(:).'
    rowsUse = rows(ri.active(rows, ii) == 1);
    Q = femurFramePts(ri, 1:9, ii, ctx);
    P = [P; Q(rowsUse, :)]; %#ok<AGROW>
end
end

function P = femurFramePts(ri, rows, poses, ctx)
% Route points (t1-frame rows 6:9 converted t1 -> ICR -> femur) for the
% given poses, concatenated as [x y] rows.
P = zeros(0, 2);
for ii = poses(:).'
    Q = ri.raw(:, :, ii);
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), Q(j, :));
        Q(j, :) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end
    P = [P; Q(rows, 1:2)]; %#ok<AGROW>
end
end

function P = candFramePts(ri, rows, poses, ctx)
P = zeros(0, 2);
for ii = poses(:).'
    Q = ri.candidateRaw(:, :, ii);
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), Q(j, :));
        Q(j, :) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end
    P = [P; Q(rows, 1:2)]; %#ok<AGROW>
end
end

function win = squareWindow(P, padFrac)
x0 = min(P(:,1)); x1 = max(P(:,1));
y0 = min(P(:,2)); y1 = max(P(:,2));
wx = (1 + padFrac) * max(x1 - x0, 1e-3);
wy = (1 + padFrac) * max(y1 - y0, 1e-3);
s  = max(wx, wy);
cx = (x0 + x1) / 2; cy = (y0 + y1) / 2;
win = [cx - s/2, cx + s/2, cy - s/2, cy + s/2];
end

function make_route_check_figure(riCam, designTitle, outPng, ctx, phiD, ...
    riBase, zoomWin, COL, RF, iRel3)
% One route-check figure: (1,1) full-flexion overview; (1,2)/(2,1)/(2,2)
% shared-window distal zooms at full flexion, the p3 release frame, and
% full extension. Stock route = black context; cam route = indigo
% squares dashed; cam CAD points = indigo open diamonds (hue + dash +
% marker shape all differ, CVD-safe per Ben's 2026-09-21 standard).

iOver = 1;
iEnd  = numel(phiD);
poses = struct('ii', {iOver, iOver, iRel3, iEnd}, ...
    'zoom', {false, true, true, true}, ...
    'note', {'full flexion, overview', 'full flexion, distal zoom', ...
             'p3 releases next, distal zoom', 'full extension, distal zoom'});

fig = figure('Color', 'w', 'Units', 'inches', 'Position', [0 0 6.5 9.2]);
tGeo = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
thPlot = linspace(0, 2*pi, 200).';

for q = 1:4
    ii = poses(q).ii;
    ax = nexttile(tGeo, q);
    hold(ax, 'on')

    P = riBase.raw(:, :, ii);              % stock route -> femur frame
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), P(j, :));
        P(j, :) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end
    Pc = riCam.raw(:, :, ii);              % cam route -> femur frame
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), Pc(j, :));
        Pc(j, :) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end
    Pk = riCam.candidateRaw(:, :, ii);     % cam CAD points -> femur frame
    for j = 6:9
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), Pk(j, :));
        Pk(j, :) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end

    % --- clearance envelopes (Tol palette, ancestor styling) ----------
    C = ctx.geo.femurCylCenter; R = ctx.geo.femurCylClearRadius;
    plot(ax, C(1)+R*cos(thPlot), C(2)+R*sin(thPlot), '-', ...
        'Color', COL.gold, 'LineWidth', 1.4);
    plot(ax, [ctx.geo.femurLineX ctx.geo.femurLineX], ...
        ctx.geo.femurLineY, '-', 'Color', COL.orange, 'LineWidth', 1.4);
    Qb = ctx.geo.femurOffsetBoundary;
    plot(ax, [Qb(:,1); Qb(1,1)], [Qb(:,2); Qb(1,2)], '-', ...
        'Color', COL.lorange, 'LineWidth', 1.4);
    if isfield(ctx.geo, 'femurCondyleClipY')
        plot(ax, [ctx.geo.femurCondyleClipX ctx.geo.femurCondyleClipX], ...
            ctx.geo.femurCondyleClipY, '-', 'Color', COL.lorange, ...
            'LineWidth', 1.4);
    end
    L = [ctx.geo.tibiaLowerCenter(1) + ...
        ctx.geo.tibiaLowerClearRadius*cos(thPlot), ...
        ctx.geo.tibiaLowerCenter(2) + ...
        ctx.geo.tibiaLowerClearRadius*sin(thPlot), zeros(numel(thPlot),1)];
    Lf = zeros(size(L));
    for kk = 1:size(L,1)
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), L(kk,:));
        Lf(kk,:) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end
    plot(ax, Lf(:,1), Lf(:,2), '-', 'Color', COL.pink, 'LineWidth', 1.4);
    U = [ctx.geo.tibiaUpperCenter(1) + ...
        ctx.geo.tibiaUpperClearRadius*cos(thPlot), ...
        ctx.geo.tibiaUpperCenter(2) + ...
        ctx.geo.tibiaUpperClearRadius*sin(thPlot), zeros(numel(thPlot),1)];
    Uf = zeros(size(U));
    for kk = 1:size(U,1)
        qICR = RowVecTrans(ctx.T_ICR_t1(:, :, ii), U(kk,:));
        Uf(kk,:) = RowVecTrans(ctx.T_Pam(:, :, ii), qICR);
    end
    plot(ax, Uf(:,1), Uf(:,2), '-', 'Color', COL.mag, 'LineWidth', 1.4);

    % --- bend-radius lines (p2 about femur cylinder, p8 about tibia) ---
    tibC = RowVecTrans(ctx.T_ICR_t1(:, :, ii), [ctx.geo.tibiaLowerCenter, 0]);
    tibC = RowVecTrans(ctx.T_Pam(:, :, ii), tibC);
    rX = NaN; rY = NaN;
    if riCam.active(2, ii)
        rX = [rX, C(1), Pc(2,1), NaN]; rY = [rY, C(2), Pc(2,2), NaN];
    end
    if riCam.active(8, ii)
        rX = [rX, tibC(1), Pc(8,1)]; rY = [rY, tibC(2), Pc(8,2)];
    end
    plot(ax, rX, rY, '--', 'Color', COL.purple, 'LineWidth', 1.4);

    % --- stock route (black context) -----------------------------------
    plot(ax, P(:,1), P(:,2), 'o-', 'Color', COL.black, 'LineWidth', 1.8, ...
        'MarkerSize', 5.5, 'MarkerFaceColor', COL.black, ...
        'MarkerEdgeColor', 'none');

    % --- cam design route (indigo squares, dashed) ----------------------
    plot(ax, Pc(:,1), Pc(:,2), 's--', 'Color', COL.cam, 'LineWidth', 2.2, ...
        'MarkerSize', 6, 'MarkerFaceColor', COL.cam, ...
        'MarkerEdgeColor', 'none');

    % --- cam CAD candidate points p6:p8 (open diamonds) -----------------
    plot(ax, Pk(6:8,1), Pk(6:8,2), 'd', 'Color', COL.cam, ...
        'LineWidth', 1.6, 'MarkerSize', 7.5, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', COL.cam, 'LineStyle', 'none');

    % --- endpoints ------------------------------------------------------
    scatter(ax, Pc(1,1), Pc(1,2), 60, 'Marker', 's', ...
        'MarkerFaceColor', COL.gold, 'MarkerEdgeColor', 'none');
    scatter(ax, Pc(9,1), Pc(9,2), 60, 'Marker', 'd', ...
        'MarkerFaceColor', COL.orange, 'MarkerEdgeColor', 'none');

    % --- point labels (parity-staggered so tight clusters do not
    % collide: odd points label above-right, even below-right of the
    % marker; ~6 pt vertical split between neighbors) ---------------------
    spanWin = zoomWin(2) - zoomWin(1);
    if poses(q).zoom
        for j = 1:9
            if riCam.active(j, ii)
                yPar = 1; if mod(j, 2) == 0, yPar = -1; end
                text(ax, Pc(j,1) + 0.024*spanWin, ...
                    Pc(j,2) + yPar*0.036*spanWin, sprintf('p%d', j), ...
                    'Color', COL.cam, 'FontName', 'Arial', ...
                    'FontSize', RF.plab, 'FontWeight', 'bold', ...
                    'Interpreter', 'none', 'Clipping', 'on', ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle')
            end
        end
    else
        text(ax, Pc(1,1) + 0.024*spanWin, Pc(1,2), 'p1', 'Color', COL.cam, ...
            'FontName', 'Arial', 'FontSize', RF.plab, 'FontWeight', ...
            'bold', 'Interpreter', 'none', 'Clipping', 'on', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
        text(ax, Pc(9,1) + 0.024*spanWin, Pc(9,2), 'p9 (pEnd)', ...
            'Color', COL.cam, 'FontName', 'Arial', 'FontSize', RF.plab, ...
            'FontWeight', 'bold', 'Interpreter', 'none', ...
            'Clipping', 'on', 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle')
    end

    axis(ax, 'equal')
    if poses(q).zoom
        xlim(ax, zoomWin(1:2)); ylim(ax, zoomWin(3:4));
    end
    xlabel(ax, 'Femur-frame x, m'); ylabel(ax, 'Femur-frame y, m')
    % Two-line compact titles: a single wide line collides with the
    % neighboring panel's title at 3 in column width.
    title(ax, {sprintf('\\theta_k = %+.1f^{\\circ}', phiD(ii)), ...
        poses(q).note}, 'Interpreter', 'tex')
    set(ax, 'FontName', 'Arial', 'FontSize', RF.tick, 'FontAngle', ...
        'normal', 'FontWeight', 'normal', 'LineWidth', 1.2)
    set(ax.XAxis, 'FontName', 'Arial', 'FontSize', RF.tick)
    set(ax.YAxis, 'FontName', 'Arial', 'FontSize', RF.tick)
    xlabel(ax, ax.XLabel.String, 'FontName', 'Arial', ...
        'FontSize', RF.lab, 'FontAngle', 'normal')
    ylabel(ax, ax.YLabel.String, 'FontName', 'Arial', ...
        'FontSize', RF.lab, 'FontAngle', 'normal')
    title(ax, ax.Title.String, 'FontName', 'Arial', 'FontSize', RF.tit, ...
        'FontAngle', 'normal', 'FontWeight', 'normal', 'Interpreter', 'tex')
end

% --- shared legend, full-width third row -------------------------------
axLeg = nexttile(tGeo, 5, [1, 2]);
axis(axLeg, 'off'); hold(axLeg, 'on')
h = gobjects(11, 1);
h(1)  = plot(axLeg, NaN, NaN, '-',  'Color', COL.gold);
h(2)  = plot(axLeg, NaN, NaN, '-',  'Color', COL.orange);
h(3)  = plot(axLeg, NaN, NaN, '-',  'Color', COL.lorange);
h(4)  = plot(axLeg, NaN, NaN, '-',  'Color', COL.pink);
h(5)  = plot(axLeg, NaN, NaN, '-',  'Color', COL.mag);
h(6)  = plot(axLeg, NaN, NaN, '--', 'Color', COL.purple);
h(7)  = plot(axLeg, NaN, NaN, 'o-', 'Color', COL.black, 'LineWidth', 1.8, ...
    'MarkerSize', 6, 'MarkerFaceColor', COL.black, 'MarkerEdgeColor', 'none');
h(8)  = plot(axLeg, NaN, NaN, 's--', 'Color', COL.cam, 'LineWidth', 2.2, ...
    'MarkerSize', 6, 'MarkerFaceColor', COL.cam, 'MarkerEdgeColor', 'none');
h(9)  = plot(axLeg, NaN, NaN, 'd', 'Color', COL.cam, 'LineWidth', 1.6, ...
    'MarkerSize', 7.5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', ...
    COL.cam, 'LineStyle', 'none');
h(10) = scatter(axLeg, NaN, NaN, 60, 'Marker', 's', ...
    'MarkerFaceColor', COL.gold, 'MarkerEdgeColor', 'none');
h(11) = scatter(axLeg, NaN, NaN, 60, 'Marker', 'd', ...
    'MarkerFaceColor', COL.orange, 'MarkerEdgeColor', 'none');
lgd = legend(axLeg, h, { ...
    'Femur cylinder clearance', ...
    'Femur straight-wall clearance', ...
    'Condyle normal-offset boundary', ...
    'Tibia lower clearance circle', ...
    'Tibia upper clearance circle', ...
    'Local bend radii (p2, p8)', ...
    'Stock 9-point route (context)', ...
    'Cam design route p1:p9', ...
    'Cam CAD points p6:p8 (open)', ...
    'Cam origin p1', ...
    'Cam distal endpoint p9 (pEnd)'}, ...
    'Location', 'north', 'NumColumns', 2, 'AutoUpdate', 'off');
styleLeg(lgd, RF)

title(tGeo, designTitle, 'FontName', 'Arial', 'FontSize', RF.sgt, ...
    'FontAngle', 'normal', 'FontWeight', 'normal')

exportgraphics(fig, outPng, 'Resolution', 200);
fprintf('wrote %s\n', outPng);
end
