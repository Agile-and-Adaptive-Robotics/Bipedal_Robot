%gen_Ext10mm_52cm_2panel_20261023.m -- Figure 8 becomes a two-panel MATLAB
% figure (Ben's display-item restructure, 2026-10-23): the supplement's
% pass-through photo moves in as panel (A), the existing single-panel torque
% plot (gen_Ext10mm_52cm.m machinery, initial extensor terms) regenerates as
% panel (B). Letters (A)/(B) per the margin-left convention (constants from
% gen_ExtPin6Panel.m; letters straddle the box top, tiles carry no titles).
% Output: Knee_Torque_Test/Figures/Ext10mm_52cm_20261023.eps + QC png.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('EXTX3_TRANS', '');

c1 = [1 0.843 0];        % #FFD700 gold    - Hybrid
c5 = [0.804 0.204 0.710];% #CD34B5 magenta - Optimized prediction
c7 = [0 0 1];            % #0000FF indigo  - Measured
cgray = [0.4 0.4 0.4];   % Original dashed
sz = 60;

[fB, bpaB] = minimizeExt(-0.0109, 1.02e5, 1.06e4, 0.123, 1); %#ok<ASGLU>
fprintf('bio ext GoF (initial coeffs): %.4f/%.4f/%.4f\n', fB(1,1), fB(1,2), fB(1,3));
b = bpaB(1);

photo = imread(fullfile(outdir, 'KneeExt_10mm_passthru.png'));

figT = figure('Name','Ext10mm_52cm_2panel','Color','w');
figT.Position = [70 70 980 460];
tl = tiledlayout(figT, 1, 2, 'TileSpacing','loose','Padding','loose');
axT = gobjects(1,2);

% (A) photo
axT(1) = nexttile(1);
imshow(photo, 'Parent', axT(1));
axis(axT(1), 'off');

% (B) torque plot (identical series to the previous single panel)
axT(2) = nexttile(2); hold(axT(2), 'on');
scatter(axT(2), b.A_h, b.M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName', 'Hybrid');
scatter(axT(2), b.Aexp, b.Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName', 'Measured');
plot(axT(2), b.Ak, b.M, '--', 'Color', cgray, 'LineWidth', 2, 'DisplayName', 'Original');
plot(axT(2), b.Ak, b.M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
xlabel(axT(2), '\bf \theta_{k} , \circ', 'Interpreter', 'tex');
ylabel(axT(2), '\bf Torque, N\cdotm', 'Interpreter', 'tex');
set(axT(2), 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
    'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XLim', [-125 15], 'XTick', [-120 -90 -60 -30 0], 'YLim', [0 20], ...
    'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
lg = legend(axT(2)); lg.Location = 'northeast'; lg.FontSize = 8;

% letters (margin-left convention; photo panel letter straddles its box top)
ann = {'(A)','(B)'};
drawnow;
ylL = 0.008; ylR = 0.024;
xL1 = 0.040;
if (xL1 < ylR) && (xL1 + 0.032 > ylL)
    Sx = (xL1 + 0.032 + 0.008) - ylL;
    tp = tl.Position; tl.Position = [tp(1) + Sx, tp(2), tp(3) - Sx, tp(4)];
    drawnow;
    ylL = ylL + Sx; ylR = ylR + Sx;
end
fprintf('AUDIT margin: yLabel band x[%.4f %.4f]  (A) left=%.4f  gap=%.4f\n', ...
    ylL, ylR, xL1, xL1 - ylR);
apos = zeros(2,2); outL = zeros(2,1);
for j = 1:2
    pos0 = axT(j).Position;
    outL(j) = axT(j).OuterPosition(1);
    tcY = pos0(2) + pos0(4);
    if j == 1, xj = xL1; else, xj = outL(j) - 0.037; end
    a = annotation(figT,'textbox',[xj, tcY - 0.020, 0.032, 0.04], ...
        'String',['\bf ' ann{j}], 'FontSize',12,'FontName','Arial', ...
        'FontWeight','bold','EdgeColor','none','HorizontalAlignment','left', ...
        'VerticalAlignment','top');
    try
        drawnow;
        ae = a.Extent;
        a.Position(2) = a.Position(2) + (tcY - (ae(2) + ae(4)/2));
    catch
    end
    apos(j,:) = a.Position(1:2);
    fprintf('AUDIT tile %d %s: letter@(%.3f,%.3f) outerL=%.3f\n', ...
        j, ann{j}, apos(j,1), apos(j,2), outL(j));
end

exportgraphics(figT, fullfile(outdir, 'Ext10mm_52cm_20261023.eps'), 'ContentType','vector');
exportgraphics(figT, fullfile(qcdir, 'QC_Ext10mm_52cm_2panel_20261023.png'), 'Resolution',150);
fprintf('WROTE Ext10mm_52cm_20261023.eps\n'); close(figT);
disp('DONE_Ext10mm_2panel');
