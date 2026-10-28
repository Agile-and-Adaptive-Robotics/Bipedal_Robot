%gen_ExtensorRoute_A4.m — F2 of EDIT_SPEC_20261020: replots the redesigned
%(cam-drum) extensor route as a 3x3 A4-portrait figure from the cam result mat
%of record. 8 knee-pose tiles + 1 legend/info tile. Campaign style (spec §0):
%Ariel[sic]-family Arial 12 bold axes/titles, LineWidth 2, legend 8 pt, panel
%letters (A)-(H) as margin labels left of each tile centered on the title
%line. Distinct hues per element class: route = indigo, femur envelopes =
%gold/orange/light-orange, tibia envelopes = pink/magenta, cam contacts =
%black markers with magenta-2 dashed bend-radius lines.
%
% Route geometry source (no invented geometry): drawing logic ported from
% Code/Matlab/Mesh_Optimization/Export_AppendixRouteGeo.m (the established
% route-geometry figure, itself a port of Opt_run_Ext.m lines 493-756):
% envelopes from ctx.geo, route rows from routeInfo.raw with p6:p9 converted
% t1 -> ICR -> femur via RowVecTrans(T_ICR_t1, then T_Pam).
%
% Exports: ProofFinal/Figures/extensor_cam_route_bones.eps (vector) +
% tmp/paper_bridge/figscripts/QC_ExtensorRoute_A4.png (150 dpi).

root = 'D:/GitHub/Bipedal_Robot';
outdir = 'D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Archive/appendix_c_corrections_20261028/figures30';
qcdir = 'D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Archive/appendix_c_corrections_20261028';
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));

matFile = fullfile(root, 'Code/Matlab/Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S = load(matFile);
ctx = S.ctxCamBest;
ri  = S.predBest.routeInfo;
phiD = ctx.phiD(:);
N = numel(phiD);
geo = ctx.geo;
fprintf('Loaded %s\n', matFile);
fprintf('phiD: %d poses, [%.1f .. %.1f] deg; route raw %s, active %s\n', ...
    N, min(phiD), max(phiD), mat2str(size(ri.raw)), mat2str(size(ri.active)));

%% Campaign palette (Code/Matlab/Colors.m, Tol accessible set)
c1 = [1.0000 0.8431 0];        % gold         - femur cylinder clearance
c2 = [1.0000 0.6941 0.3059];   % orange       - femur wall clearance
c3 = [0.9804 0.5294 0.4588];   % light orange - condyle offset boundary
c4 = [0.9176 0.3725 0.5804];   % pink         - tibia lower clearance
c5 = [0.8039 0.2039 0.7098];   % magenta      - tibia upper clearance
c6 = [0.6157 0.0078 0.8431];   % magenta 2    - local bend-radius lines
c7 = [0 0 1];                  % indigo       - optimized route (BPA path)
cK = [0 0 0];                  % black        - cam contacts p6:p8
lineW = 2;
fontN = 'Arial';

%% Pose schedule: 8 poses.
% Preference 1 (prior-figure convention, Export_AppendixRouteGeo): full
% flexion, the frame immediately before each unique elimination transition,
% full extension. Preference 2 (spec floor): 8 evenly spaced poses when the
% transition count does not give exactly 8 tiles.
transitionIdx = find(any(ri.active(:,1:end-1) & ~ri.active(:,2:end), 1));
cand = unique([1, transitionIdx, N], 'stable');
transAngles = phiD(transitionIdx(:)).';
transRows = arrayfun(@(t) find(ri.active(:,t) & ~ri.active(:,t+1)), transitionIdx);
fprintf('Route-state transitions (authoritative, from ri.active): ');
fprintf('p%d@%.1fdeg  ', [transRows; transAngles]);
fprintf('\n');
fprintf('eliminatedAngleD field of record: ');
rowsE = find(isfinite(ri.eliminatedAngleD));
for kk = 1:numel(rowsE)
    fprintf('p%d@%.1f  ', rowsE(kk), ri.eliminatedAngleD(rowsE(kk)));
end
fprintf('\n');

if numel(cand) == 8
    poseIdx = cand;
    schedNote = 'poses: flexion + pre-elimination frames + extension';
else
    targets = linspace(min(phiD), max(phiD), 8);
    poseIdx = zeros(1,8);
    for k = 1:8
        [~, poseIdx(k)] = min(abs(phiD - targets(k)));
    end
    poseIdx = unique(poseIdx, 'stable');
    schedNote = sprintf('poses: 8 evenly spaced over design range (%d transition frames found)', numel(cand)-2);
end
nPoses = numel(poseIdx);
fprintf('POSE SCHEDULE (%s): ', schedNote);
fprintf('%.1f  ', phiD(poseIdx));
fprintf('deg\n');
if nPoses ~= 8
    error('gen_ExtensorRoute_A4: pose schedule gave %d poses, need 8.', nPoses);
end

%% Figure: A4 portrait print size, 3x3 tiles
fig = figure('Name','ExtensorRoute_redesigned_A4','Color','w', ...
    'Units','centimeters','Position',[1 1 17 24]);
t = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
thPlot = linspace(0, 2*pi, 200).';

%% Pre-pass: shared axis limits over all pose geometry + static envelopes
allX = []; allY = [];
for k = 1:nPoses
    ii = poseIdx(k);
    P = routeToFemurFrame(ri.raw(:,:,ii), ctx, ii);
    allX = [allX; P(:,1)]; allY = [allY; P(:,2)]; %#ok<AGROW>
    for ccIdx = 1:2
        if ccIdx == 1
            Cc = geo.tibiaLowerCenter; Rc = geo.tibiaLowerClearRadius;
        else
            Cc = geo.tibiaUpperCenter; Rc = geo.tibiaUpperClearRadius;
        end
        L = [Cc(1)+Rc*cos(thPlot), Cc(2)+Rc*sin(thPlot), zeros(numel(thPlot),1)];
        Lf = t1ToFemur(L, ctx, ii);
        allX = [allX; Lf(:,1)]; allY = [allY; Lf(:,2)]; %#ok<AGROW>
    end
end
allX = [allX; geo.femurCylCenter(1)+geo.femurCylClearRadius*[-1;1]];
allY = [allY; geo.femurCylCenter(2)+geo.femurCylClearRadius*[-1;1]];
allX = [allX; geo.femurLineX*[1;1]]; allY = [allY; geo.femurLineY(:)];
allX = [allX; geo.femurOffsetBoundary(:,1)]; allY = [allY; geo.femurOffsetBoundary(:,2)];
xPad = 0.04*range(allX); yPad = 0.04*range(allY);
gxLim = [min(allX)-xPad, max(allX)+xPad];
gyLim = [min(allY)-yPad, max(allY)+yPad];
fprintf('SHARED LIMITS: x[%.4f %.4f] y[%.4f %.4f] m\n', gxLim, gyLim);

%% Pose tiles
axT = gobjects(1,nPoses);
hT  = gobjects(1,nPoses);   % legend template handles (tile 1)
for k = 1:nPoses
    ii = poseIdx(k);
    ax = nexttile(k); hold(ax, 'on');
    axT(k) = ax;

    P = routeToFemurFrame(ri.raw(:,:,ii), ctx, ii);
    active = ri.active(:,ii);

    h = gobjects(10,1);

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

    % Tibia lower clearance circle -> femur frame (pink)
    L = [geo.tibiaLowerCenter(1)+geo.tibiaLowerClearRadius*cos(thPlot), ...
         geo.tibiaLowerCenter(2)+geo.tibiaLowerClearRadius*sin(thPlot), ...
         zeros(numel(thPlot),1)];
    Lf = t1ToFemur(L, ctx, ii);
    h(4) = plot(ax, Lf(:,1), Lf(:,2), '-', 'Color', c4, 'LineWidth', lineW);

    % Tibia upper clearance circle -> femur frame (magenta)
    U = [geo.tibiaUpperCenter(1)+geo.tibiaUpperClearRadius*cos(thPlot), ...
         geo.tibiaUpperCenter(2)+geo.tibiaUpperClearRadius*sin(thPlot), ...
         zeros(numel(thPlot),1)];
    Uf = t1ToFemur(U, ctx, ii);
    h(5) = plot(ax, Uf(:,1), Uf(:,2), '-', 'Color', c5, 'LineWidth', lineW);

    % Local p2/p8 bend radius lines (dashed, magenta 2) — prior-figure rule
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

    % Optimized route p1:p9 (indigo, filled circle markers, no edge)
    h(7) = plot(ax, P(:,1), P(:,2), 'o-', 'Color', c7, 'LineWidth', lineW, ...
        'MarkerSize', 5, 'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none');

    % Cam contacts p6:p8 (black) — the clip-on-cam tibia part surface points
    h(8) = scatter(ax, P(6:8,1), P(6:8,2), 70, 'filled', ...
        'MarkerFaceColor', cK, 'MarkerEdgeColor', 'none');

    % Endpoints: femur-side origin p1 (gold square), insertion pEnd (orange diamond)
    h(9)  = scatter(ax, P(1,1), P(1,2), 90, 'Marker', 's', ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none');
    h(10) = scatter(ax, P(9,1), P(9,2), 90, 'Marker', 'd', ...
        'MarkerFaceColor', c2, 'MarkerEdgeColor', 'none');

    % Label active route points (10 pt floor, Mesh_Optimization journal
    % standard). Per-point label offsets declutter the tight p6:p8 cam-contact
    % cluster (QC_ExtensorRoute_A4 round 1: labels overlapped there).
    lblOff = 0.001 * [ 2  2; 2  2; 2  2; 2  2; 2  2; ...
                       5  3; 6 -1; 4 -6; 2  2];
    for j = 1:9
        if active(j)
            text(ax, P(j,1)+lblOff(j,1), P(j,2)+lblOff(j,2), ...
                sprintf(' p%d', j), ...
                'FontName', fontN, 'FontSize', 10, 'FontWeight', 'bold', ...
                'Interpreter', 'none');
        end
    end

    if k == 1, hT = h; end

    axis(ax, 'equal');
    xlim(ax, gxLim); ylim(ax, gyLim);
    title(ax, sprintf('\\bf \\theta_{k} = %+.0f\\circ', phiD(ii)), ...
        'Interpreter', 'tex');
    set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', fontN, ...
        'LineWidth', lineW, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
end



axLeg=nexttile(9);axis(axLeg,'off');
%% Panel letters (A)-(H): margin labels LEFT of each pose tile, vertically
% centered on the tile TITLE line (campaign rule, spec §0 / gen_ExtPin6Panel).
drawnow;
ann = arrayfun(@(k) sprintf('(%c)', 'A' + k - 1), 1:nPoses, 'UniformOutput', false);
xL1 = 0.040;                                   % column-1 target left edge
tp = t.Position;
ylR = tp(1) + 0.016;                           % rotated y-label band right edge
if (xL1 < ylR) && (xL1 + 0.032 > tp(1))        % y-label overlaps: open the margin
    Sh = (xL1 + 0.032 + 0.006) - tp(1);
    t.Position = [tp(1) + Sh, tp(2), tp(3) - Sh, tp(4)];
    drawnow;
end
fprintf('AUDIT margin: layout left=%.4f  (A)/(D)/(G) left=%.4f  gap=%.4f\n', ...
    t.Position(1), xL1, xL1 - (t.Position(1) + 0.016));
for j = 1:nPoses
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
    if j == 1 || j == 4 || j == 7
        xj = xL1;
    else
        xj = outL - 0.034;
    end
    a = annotation(fig, 'textbox', [xj, tcY - 0.020, 0.032, 0.04], ...
        'String', ['\bf ' ann{j}], 'FontSize', 12, 'FontName', fontN, ...
        'FontWeight', 'bold', 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    try
        drawnow;
        ae = a.Extent;
        a.Position(2) = a.Position(2) + (tcY - (ae(2) + ae(4)/2));
    catch
    end
    fprintf('AUDIT tile %d %s: letter@(%.3f,%.3f) outerL=%.3f titleC=%.4f\n', ...
        j, ann{j}, a.Position(1), a.Position(2), outL, tcY);
end

for kb=1:nPoses, appendix_c_bone_background(axT(kb),ctx,poseIdx(kb)); end
%% Exports
exportgraphics(fig, fullfile(outdir, 'extensor_cam_route_bones.eps'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(qcdir, 'QC_ExtensorRoute_A4.png'), ...
    'Resolution', 150);
fprintf('WROTE extensor_cam_route_bones.eps + QC_ExtensorRoute_A4.png\n');
exportgraphics(fig,fullfile(outdir,'extensor_cam_route_bones.pdf'),'ContentType','vector');
close(fig);
disp('DONE_gen_ExtensorRoute_A4');

%% ---- helpers ---------------------------------------------------------
function P = routeToFemurFrame(Praw, ctx, ii)
% Native-frame route rows -> one common femur frame: p6:p9 are native t1
% coordinates; convert t1 -> ICR -> femur (Export_AppendixRouteGeo rule).
P = Praw;
for j = 6:9
    qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), P(j,:));
    P(j,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
end
end

function F = t1ToFemur(L, ctx, ii)
% Row-wise t1 -> ICR -> femur for an N x 3 array.
F = zeros(size(L));
for kk = 1:size(L,1)
    qICR = RowVecTrans(ctx.T_ICR_t1(:,:,ii), L(kk,:));
    F(kk,:) = RowVecTrans(ctx.T_Pam(:,:,ii), qICR);
end
end
