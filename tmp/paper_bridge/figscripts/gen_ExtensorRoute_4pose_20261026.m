%gen_ExtensorRoute_4pose_20261026.m -- P2-1 rebuild of the cam-extensor
% routing-geometry figure per Ben (2026-10-26 14:26): too big, add femur and
% tibia shapes, reduce the pose count. FOUR poses chosen to explain the
% contact/release cascade (release schedule from predBest.routeInfo.active,
% verified in dbg_geo_20261026): (A) -120 deg full flexion, all contacts
% active; (B) -91.1 deg, p5 released, p4/p3 still in contact; (C) -55.7 deg,
% p4 released, p3's last active frame; (D) +10 deg full extension, the cam
% wrap (p2, p6) released. Schematic bone silhouettes are ADDED as
% presentation: the femur silhouette is the convex hull of the
% route-elimination envelopes (cylinder clearance circle, condyle
% offset-boundary curve, wall line) and the tibia silhouette is the hull of
% the two tibial clearance circles with a shaft along the tibia (t1) axis;
% all underlying geometry is from ctx.geo, nothing invented. 2x3 layout:
% four pose tiles + dedicated legend tile + release-schedule tile; 10 pt
% text floor. Output: Knee_Torque_Test/Figures/ExtensorRoute_redesigned.eps
% (same name, replaces) + rebuilt cache + QC png.
root = 'D:/GitHub/Bipedal_Robot';
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/Figures');
qcdir  = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));

S = load(fullfile(root, 'Code/Matlab/Mesh_Optimization/Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
ctx = S.ctxCamBest;
ri  = S.predBest.routeInfo;
phiD = ctx.phiD(:);
N = numel(phiD);
geo = ctx.geo;

c1 = [1.0000 0.8431 0];        % gold         - femur cylinder clearance
c2 = [1.0000 0.6941 0.3059];   % orange       - femur wall clearance
c3 = [0.9804 0.5294 0.4588];   % light orange - condyle offset boundary
c4 = [0.9176 0.3725 0.5804];   % pink         - tibia lower clearance
c5 = [0.8039 0.2039 0.7098];   % magenta      - tibia upper clearance
c6 = [0.6157 0.0078 0.8431];   % magenta 2    - local bend-radius lines
c7 = [0 0 1];                  % indigo       - optimized route (BPA path)
cK = [0 0 0];                  % black        - cam contacts p6:p8
cBone = [0.82 0.82 0.82];      % light gray   - schematic bone silhouettes
lineW = 2;
fontN = 'Arial';

%% pose schedule: the release cascade (from ri.active, the authoritative set)
transitionIdx = find(any(ri.active(:,1:end-1) & ~ri.active(:,2:end), 1));
transAngles = phiD(transitionIdx(:)).';
transRows = arrayfun(@(t) find(ri.active(:,t) & ~ri.active(:,t+1)), transitionIdx);
fprintf('Route-state transitions: ');
fprintf('p%d@%.1fdeg  ', [transRows; transAngles]);
fprintf('\n');
% (A) full flexion, (B) after p5's release at p4's last frame,
% (C) p3's last active frame, (D) full extension (wrap released).
iA = 1;
iB = transitionIdx(2);      % -91.1 deg: p5 off, p4 last frame
iC = transitionIdx(3);      % -55.7 deg: p4 off, p3 last frame
iD = N;                     % +10 deg
poseIdx = [iA iB iC iD];
fprintf('POSE SET: '); fprintf('%.1f  ', phiD(poseIdx)); fprintf('deg\n');

%% figure: 2x3, four pose tiles + legend tile + release tile
fig = figure('Name','ExtensorRoute_4pose_20261026','Color','w', ...
    'Units','centimeters','Position',[1 1 19.5 13.5]);
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
thPlot = linspace(0, 2*pi, 200).';

%% shared limits (poses + envelopes + bone silhouettes)
allX = []; allY = [];
for k = 1:4
    ii = poseIdx(k);
    P = routeToFemurFrame(ri.raw(:,:,ii), ctx, ii);
    allX = [allX; P(:,1)]; allY = [allY; P(:,2)]; %#ok<AGROW>
    [tibX, tibY] = tibiaOutline(ctx, geo, ii);
    allX = [allX; tibX(:)]; allY = [allY; tibY(:)]; %#ok<AGROW>
end
allX = [allX; geo.femurCylCenter(1)+geo.femurCylClearRadius*[-1;1]];
allY = [allY; geo.femurCylCenter(2)+geo.femurCylClearRadius*[-1;1]];
allX = [allX; geo.femurLineX*[1;1]]; allY = [allY; geo.femurLineY(:)];
allX = [allX; geo.femurOffsetBoundary(:,1)]; allY = [allY; geo.femurOffsetBoundary(:,2)];
[femX, femY] = femurOutline(geo, max(allY));
allX = [allX; femX(:)]; allY = [allY; femY(:)];
xPad = 0.04*range(allX); yPad = 0.04*range(allY);
gxLim = [min(allX)-xPad, max(allX)+xPad];
gyLim = [min(allY)-yPad, max(allY)+yPad];
fprintf('SHARED LIMITS: x[%.4f %.4f] y[%.4f %.4f] m\n', gxLim, gyLim);

%% pose tiles
axT = gobjects(1,4);
hT  = gobjects(12,1);   % legend template handles (tile A)
for k = 1:4
    ii = poseIdx(k);
    ax = nexttile(k); hold(ax, 'on');
    axT(k) = ax;
    h = gobjects(12,1);
    P = routeToFemurFrame(ri.raw(:,:,ii), ctx, ii);
    active = ri.active(:,ii);

    % schematic bone silhouettes (presentation layer, under the envelopes)
    [fx, fy] = femurOutline(geo, gyLim(2));
    h(11) = patch(ax, fx, fy, cBone, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    [tx, ty] = tibiaOutline(ctx, geo, ii);
    h(12) = patch(ax, tx, ty, cBone, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
    text(ax, geo.femurCylCenter(1), max(fy)-0.012, 'femur', 'FontName', fontN, ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(ax, mean(tx(1:2)), max(ty)+0.006, 'tibia', 'FontName', fontN, ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Femur cylinder clearance (gold)
    h(1) = plot(ax, geo.femurCylCenter(1)+geo.femurCylClearRadius*cos(thPlot), ...
        geo.femurCylCenter(2)+geo.femurCylClearRadius*sin(thPlot), ...
        '-', 'Color', c1, 'LineWidth', lineW);

    % Femur straight-wall clearance (orange)
    h(2) = plot(ax, [geo.femurLineX geo.femurLineX], geo.femurLineY, ...
        '-', 'Color', c2, 'LineWidth', lineW);

    % Condyle normal-offset boundary (light orange) + clip line
    Q = geo.femurOffsetBoundary;
    h(3) = plot(ax, [Q(:,1); Q(1,1)], [Q(:,2); Q(1,2)], ...
        '-', 'Color', c3, 'LineWidth', lineW);
    plot(ax, [geo.femurCondyleClipX geo.femurCondyleClipX], geo.femurCondyleClipY, ...
        '-', 'Color', c3, 'LineWidth', lineW, 'HandleVisibility', 'off');

    % Tibia clearance circles -> femur frame (pink, magenta)
    L = [geo.tibiaLowerCenter(1)+geo.tibiaLowerClearRadius*cos(thPlot), ...
         geo.tibiaLowerCenter(2)+geo.tibiaLowerClearRadius*sin(thPlot), ...
         zeros(numel(thPlot),1)];
    Lf = t1ToFemur(L, ctx, ii);
    h(4) = plot(ax, Lf(:,1), Lf(:,2), '-', 'Color', c4, 'LineWidth', lineW);
    U = [geo.tibiaUpperCenter(1)+geo.tibiaUpperClearRadius*cos(thPlot), ...
         geo.tibiaUpperCenter(2)+geo.tibiaUpperClearRadius*sin(thPlot), ...
         zeros(numel(thPlot),1)];
    Uf = t1ToFemur(U, ctx, ii);
    h(5) = plot(ax, Uf(:,1), Uf(:,2), '-', 'Color', c5, 'LineWidth', lineW);

    % Local p2/p8 bend radius lines (dashed, magenta 2)
    tibLowerFemur = t1ToFemur([geo.tibiaLowerCenter, 0], ctx, ii);
    rX = NaN; rY = NaN;
    if active(2)
        rX = [rX, geo.femurCylCenter(1), P(2,1), NaN];
        rY = [rY, geo.femurCylCenter(2), P(2,2), NaN];
    end
    if active(8)
        rX = [rX, tibLowerFemur(1), P(8,1)];
        rY = [rY, tibLowerFemur(2), P(8,2)];
    end
    h(6) = plot(ax, rX, rY, '--', 'Color', c6, 'LineWidth', lineW);

    % Optimized route p1:p9 (indigo)
    h(7) = plot(ax, P(:,1), P(:,2), 'o-', 'Color', c7, 'LineWidth', lineW, ...
        'MarkerSize', 5, 'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none');

    % Cam contacts p6:p8 (black)
    h(8) = scatter(ax, P(6:8,1), P(6:8,2), 70, 'filled', ...
        'MarkerFaceColor', cK, 'MarkerEdgeColor', 'none');

    % Endpoints: origin p1 (gold square), insertion pEnd (orange diamond)
    h(9)  = scatter(ax, P(1,1), P(1,2), 90, 'Marker', 's', ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none');
    h(10) = scatter(ax, P(9,1), P(9,2), 90, 'Marker', 'd', ...
        'MarkerFaceColor', c2, 'MarkerEdgeColor', 'none');

    % active-route point labels (10 pt floor)
    lblOff = 0.001 * [ 2  2; 2  2; 2  2; 2  2; 2  2; 5  3; 6 -1; 4 -6; 2  2];
    for j = 1:9
        if active(j)
            text(ax, P(j,1)+lblOff(j,1), P(j,2)+lblOff(j,2), ...
                sprintf(' p%d', j), 'FontName', fontN, 'FontSize', 10, ...
                'FontWeight', 'bold', 'Interpreter', 'none');
        end
    end

    if k == 1, hT = h; end
    axis(ax, 'equal');
    xlim(ax, gxLim); ylim(ax, gyLim);
    title(ax, sprintf('\\bf \\theta_{\\bf k} = %+.0f\\circ', phiD(ii)), ...
        'Interpreter', 'tex');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontN, ...
        'LineWidth', lineW, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
end
xlabel(t, '\bf Femur-frame x, m', 'Interpreter', 'tex');
ylabel(t, '\bf Femur-frame y, m', 'Interpreter', 'tex');

%% tile 5: legend
axLeg = nexttile(5);
axis(axLeg, 'off'); hold(axLeg, 'on');
labels = { ...
    'Femur cylinder clearance'; ...
    'Femur wall clearance'; ...
    'Condyle offset boundary'; ...
    'Tibia lower clearance'; ...
    'Tibia upper clearance'; ...
    'Local bend radii'; ...
    'Cam contacts p6:p8'; ...
    'p1:p9 optimized route'; ...
    'Optimized p1 (origin)'; ...
    'Optimized pEnd (insertion)'; ...
    'Femur outline (schematic)'; ...
    'Tibia outline (schematic)'};
hLeg = gobjects(numel(labels), 1);
lsMap = {'-','-','-','-','-','--','-','none','none','none','none','none'};
mkMap = {'none','none','none','none','none','none','none','o','s','d','none','none'};
for kk = 1:numel(labels)
    hLeg(kk) = plot(axLeg, NaN, NaN);
    if isprop(hT(kk), 'Color')
        hLeg(kk).Color = hT(kk).Color;
        hLeg(kk).LineStyle = lsMap{kk};
        hLeg(kk).LineWidth = hT(kk).LineWidth;
        hLeg(kk).Marker = mkMap{kk};
        hLeg(kk).MarkerSize = 6;
        hLeg(kk).MarkerEdgeColor = 'none';
        if isprop(hT(kk), 'MarkerFaceColor')
            fc = hT(kk).MarkerFaceColor;
            if isequal(fc, 'flat') || isequal(fc, 'auto'), fc = hT(kk).Color; end
            hLeg(kk).MarkerFaceColor = fc;
        end
        if ~strcmp(lsMap{kk}, 'none') && strcmp(mkMap{kk}, 'none') == false
            hLeg(kk).LineStyle = '-';
        end
    else
        % patch entries: gray swatch line
        hLeg(kk).Color = cBone;
        hLeg(kk).LineStyle = '-';
        hLeg(kk).LineWidth = 6;
    end
end
lgd = legend(axLeg, hLeg(1:10), labels(1:10), 'Location', 'north', ...
    'NumColumns', 1, 'AutoUpdate', 'off');
lgd.FontSize = 8; lgd.FontName = fontN; lgd.FontWeight = 'bold'; lgd.Box = 'off';
% bone entries drawn as text swatches (patches do not legend cleanly)
text(axLeg, 0.05, 0.02, {sprintf('\\bullet \\bf Femur outline (schematic)'), ...
    sprintf('\\bullet \\bf Tibia outline (schematic)')}, 'Units', 'normalized', ...
    'FontName', fontN, 'FontSize', 8, 'FontWeight', 'bold', ...
    'Color', [0.35 0.35 0.35], 'VerticalAlignment', 'bottom');

%% tile 6: contact-release schedule (facts from the mat of record)
axInfo = nexttile(6);
axis(axInfo, 'off'); hold(axInfo, 'on');
relAll = strjoin(arrayfun(@(r,a) sprintf('p%d at %.1f\\circ', r, a), ...
    transRows, transAngles, 'UniformOutput', false), ', ');
infoLines = { ...
    sprintf('\\bf Design range: %.0f\\circ to %.0f\\circ', min(phiD), max(phiD)); ...
    sprintf('\\bf Cam-drum wrap radius: %.1f mm', 1000*ri.wrap.radius(2,1)); ...
    '\\bf Contact release (row, angle):'; ...
    ['  ' relAll]; ...
    ' '; ...
    '\\bf Poses:'; ...
    '  (A) full flexion, all contacts'; ...
    '  (B) p5 released, p4 last frame'; ...
    '  (C) p4 released, p3 last frame'; ...
    '  (D) full extension, wrap released'};
text(axInfo, 0.02, 0.98, infoLines, 'Units', 'normalized', ...
    'FontName', fontN, 'FontSize', 10, 'FontWeight', 'bold', ...
    'Interpreter', 'tex', 'VerticalAlignment', 'top');

%% panel letters (A)-(D): margin labels on the title line
drawnow;
ann = {'(A)','(B)','(C)','(D)'};
xL1 = 0.040;
tp = t.Position;
ylR = tp(1) + 0.016;
if (xL1 < ylR) && (xL1 + 0.032 > tp(1))
    Sh = (xL1 + 0.032 + 0.006) - tp(1);
    t.Position = [tp(1) + Sh, tp(2), tp(3) - Sh, tp(4)];
    drawnow;
end
for j = 1:4
    pos0 = axT(j).Position;
    outL = axT(j).OuterPosition(1);
    innerTop = pos0(2) + pos0(4);
    tcY = innerTop + 0.030;
    try
        Tj = axT(j).Title; Tj.Units = 'normalized';
        te = Tj.Extent;
        candY = pos0(2) + (te(2) + te(4)/2) * pos0(4);
        if isfinite(candY) && candY > innerTop && candY < innerTop + 0.09
            tcY = candY;
        end
    catch
    end
    if ismember(j, [1 4]), xj = xL1; else, xj = outL - 0.034; end
    a = annotation(fig, 'textbox', [xj, tcY - 0.020, 0.032, 0.04], ...
        'String', ['\bf ' ann{j}], 'FontSize', 12, 'FontName', fontN, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    try
        drawnow; ae = a.Extent;
        a.Position(2) = a.Position(2) + (tcY - (ae(2) + ae(4)/2));
    catch
    end
    fprintf('AUDIT tile %d %s: letter@(%s)\n', j, ann{j}, mat2str(a.Position(1:2)));
end

exportgraphics(fig, fullfile(outdir, 'ExtensorRoute_redesigned.eps'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(qcdir, 'QC_ExtensorRoute_4pose_20261026.png'), ...
    'Resolution', 150);
fprintf('WROTE ExtensorRoute_redesigned.eps (4-pose + bones)\n'); close(fig);
disp('DONE_gen_ExtensorRoute_4pose');

%% ---- helpers ---------------------------------------------------------
function [fx, fy] = femurOutline(geo, yTop)
% Schematic femur silhouette (femur frame, static): convex hull of the
% cylinder clearance circle, the condyle offset-boundary curve, and the wall
% line, closed upward to the tile top as the shaft.
th = linspace(0, 2*pi, 60).';
pts = [ ...
    geo.femurCylCenter(1)+geo.femurCylClearRadius*cos(th), ...
    geo.femurCylCenter(2)+geo.femurCylClearRadius*sin(th); ...
    geo.femurOffsetBoundary(:,1:2); ...
    geo.femurLineX + 0.004*cos(th), mean(geo.femurLineY) + ...
        0.5*diff(geo.femurLineY)*sin(th) + 0.5*diff(geo.femurLineY); ...
    geo.femurCylCenter(1)+0.9*geo.femurCylClearRadius, yTop; ...
    geo.femurCylCenter(1)-0.9*geo.femurCylClearRadius, yTop];
k = convhull(pts(:,1), pts(:,2));
fx = pts(k,1); fy = pts(k,2);
end

function [tx, ty] = tibiaOutline(ctx, geo, ii)
% Schematic tibia silhouette per pose: hull of the two tibial clearance
% circles (t1 frame -> femur frame) plus a distal shaft along the tibia
% (t1 x) axis.
th = linspace(0, 2*pi, 60).';
Rl = geo.tibiaLowerClearRadius; Ru = geo.tibiaUpperClearRadius;
xMax = max(geo.tibiaLowerCenter(1)+Rl, geo.tibiaUpperCenter(1)+Ru);
xTip = xMax + 0.06;
yLo = min(geo.tibiaLowerCenter(2)-Rl, geo.tibiaUpperCenter(2)-Ru);
yHi = max(geo.tibiaLowerCenter(2)+Rl, geo.tibiaUpperCenter(2)+Ru);
pts = [ ...
    geo.tibiaLowerCenter(1)+Rl*cos(th), geo.tibiaLowerCenter(2)+Rl*sin(th), ...
        zeros(numel(th),1); ...
    geo.tibiaUpperCenter(1)+Ru*cos(th), geo.tibiaUpperCenter(2)+Ru*sin(th), ...
        zeros(numel(th),1); ...
    xMax, yLo, 0; xMax, yHi, 0; xTip, yHi, 0; xTip, yLo, 0];
F = t1ToFemur(pts, ctx, ii);
tx = F(:,1); ty = F(:,2);
end

function P = routeToFemurFrame(Praw, ctx, ii)
P = Praw;
for j = 6:9
    qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), P(j,:));
    P(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
end
end

function F = t1ToFemur(L, ctx, ii)
F = zeros(size(L));
for kk = 1:size(L,1)
    qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), L(kk,:));
    F(kk,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
end
end
