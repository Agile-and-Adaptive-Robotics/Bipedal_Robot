%gen_ExtPin5Panel_20261026.m -- record-derived extensor pinned-knee figure.
% Ben's membership correction (2026-10-26) + mat audit: the extensor
% identification of record (minimizeExt10mmX3_results.mat) ran on a
% FIVE-test pool, allBPA = [1 2 5 6 8], numBPA = 5, scores_cv 5x3: training
% {43.6, 43.2-tendon}, validation {40.0, 40.0-tendon, 46.5} (row-32
% holdCols values 1,2,8). The 45.7 cm test was NOT part of the
% identification and is dropped from the paper figure per the no-unused-
% tests ruling. Five data panels + the legend in the sixth tile; series and
% styling identical to gen_ExtPin6Panel.m at the adopted row-32 terms.
% Torque label M with knee-axis subscript. Output: Figures/
% ExtPin_group_20261026.eps + QC png.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');

c1 = [1 0.843 0];         % gold        - Hybrid
c3 = [0.980 0.529 0.459]; % light orange - Unoptimized
c5 = [0.804 0.204 0.710]; % magenta     - Optimized
c7 = [0 0 1];             % indigo      - Measured
sz = 60;

%% adopted extensor terms of record (row 32, ledger-verified)
S = load('minimizeExt10mmX3_results.mat');
g = S.filtered_results(32, S.xCols);
ledger = [-6.378022e-03, 3.997873e+04, 1.473375e+04, 1.582887e-01];
if any(abs(g(1:4) - ledger) > 1e-6 * abs(ledger))
    error('Row 32 mismatch vs ledger: got [%s]', sprintf(' %.6e', g));
end
fprintf('ADOPTED row32: X0=%.6e m  X1=%.6e  X2=%.6e N/m  X3=%.6f\n', g(1), g(2), g(3), g(4));

%% evaluator calls: call 1:9 and subset by index (tests arg unreliable)
tests = 1:9;
[f0, b0] = minimizeExtX3(0, Inf, Inf, 0, tests);          % zero-correction
[f1, b1] = minimizeExtX3(g(1), g(2), g(3), g(4), tests);  % adopted terms
pool = [1 2 5 6 8];   % record pool: 40.0, 40.0t, 43.6, 43.2t, 46.5
fprintf('\n=== Extensor pool GoF unoptimized [RMSE FVU Max] ===\n');
for i = pool
    fprintf('  rest=%.4f m: %8.4f %8.4f %8.4f\n', b0(i).rest, f0(i,:));
end
fprintf('=== Extensor pool GoF optimized [RMSE FVU Max] ===\n');
for i = pool
    fprintf('  rest=%.4f m: %8.4f %8.4f %8.4f\n', b1(i).rest, f1(i,:));
end
trG = [5 6]; vaG = [1 2 8];
fprintf('MEANS training {43.6,43.2t}: unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(trG,:),1), mean(f1(trG,:),1));
fprintf('MEANS validation {40,40t,46.5}: unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(vaG,:),1), mean(f1(vaG,:),1));

%% figure: 2x3, five data panels, legend in the sixth tile
pord = [5 6 1 2 8];
pttl = {'44 cm', '43 cm Tendon', '40 cm', '40 cm Tendon', '47 cm'};
psub = {'Training','Training','Validation','Validation','Validation'};
sty = {'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 30],'XTick',[-120 -90 -60 -30 0 30],'GridLineStyle','none'};
xmax = 30;
for i = pool
    xmax = max(xmax, max(b1(i).Aexp(b1(i).Aexp < 180)));
end
if xmax <= 30
    xup = 30;
elseif xmax <= 40
    xup = 40;
else
    xup = ceil(xmax / 5) * 5;
end
fprintf('X-LIMIT: max measured angle across pool = %.1f deg -> XLim upper %d\n', xmax, xup);
sty{find(strcmp(sty, 'XLim')) + 1} = [-120 xup];

fig = figure('Name','ExtPin_5panel_20261026','Color','w'); fig.Position = [60 60 1150 780];
t = tiledlayout(2,3,'TileSpacing','loose','Padding','loose');
allY = [];
axT = gobjects(1,5);
for k = 1:5
    i = pord(k);
    axT(k) = nexttile(k); hold on
    scatter(b1(i).A_h, b1(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName','Hybrid');
    scatter(b1(i).Aexp, b1(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
    plot(b1(i).Ak, b1(i).M, '--', 'Color', c3, 'LineWidth', 2, 'DisplayName','Unoptimized');
    plot(b1(i).Ak, b1(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Optimized');
    title(['\bf ' pttl{k}]); subtitle(['\bf ' psub{k}], 'FontName','Arial');
    set(gca, sty{:});
    allY = [allY; b1(i).M; b1(i).M_p(:,3); b1(i).Mexp; b1(i).M_h]; %#ok<AGROW>
end
% legend in the empty sixth tile
lgAx = nexttile(6);
lg = legend(lgAx, 'Location','northeast'); lg.FontSize = 10;
title(lgAx, '\bf Legend', 'FontName','Arial');
set(lgAx, 'XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[],'Box','off');
for k = 1:5
    axT(k).YLim = [0, ceil(max(allY)/5)*5];
end
ylabel(t, '\bf M_{\bf k}, N\cdotm'); xlabel(t, '\bf \theta_{\bf k} , \circ');
ann = {'(A)','(B)','(C)','(D)','(E)'};
drawnow;
ylL = 0.008; ylR = 0.024;
xL1 = 0.040;
if (xL1 < ylR) && (xL1 + 0.032 > ylL)
    Sx = (xL1 + 0.032 + 0.008) - ylL;
    tp = t.Position; t.Position = [tp(1) + Sx, tp(2), tp(3) - Sx, tp(4)];
    drawnow;
    ylL = ylL + Sx; ylR = ylR + Sx;
end
for j = 1:5
    pos0 = axT(j).Position;
    outL = axT(j).OuterPosition(1);
    innerTop = pos0(2) + pos0(4);
    tcY = innerTop + 0.030;
    try
        Tj = axT(j).Title; Tj.Units = 'normalized';
        te = Tj.Extent;
        cand = pos0(2) + (te(2) + te(4)/2) * pos0(4);
        if isfinite(cand) && cand > innerTop && cand < innerTop + 0.09
            tcY = cand;
        end
    catch
    end
    if ismember(j, [1 4]), xj = xL1; else, xj = outL - 0.037; end
    a = annotation(fig,'textbox',[xj, tcY - 0.020, 0.032, 0.04], ...
        'String',['\bf ' ann{j}], 'FontSize',12,'FontName','Arial', ...
        'FontWeight','bold','EdgeColor','none','HorizontalAlignment','left', ...
        'VerticalAlignment','top');
    try
        drawnow;
        ae = a.Extent;
        a.Position(2) = a.Position(2) + (tcY - (ae(2) + ae(4)/2));
    catch
    end
    fprintf('AUDIT tile %d %s: letter@(%s) titleC=%.4f\n', j, ann{j}, mat2str(a.Position(1:2)), tcY);
end
exportgraphics(fig, fullfile(outdir,'ExtPin_group_20261026.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(qcdir,'QC_ExtPin_5panel_20261026.png'), 'Resolution',150);
fprintf('WROTE ExtPin_group_20261026.eps (5-panel + legend tile)\n'); close(fig);
disp('DONE_ExtPin5Panel_20261026');
