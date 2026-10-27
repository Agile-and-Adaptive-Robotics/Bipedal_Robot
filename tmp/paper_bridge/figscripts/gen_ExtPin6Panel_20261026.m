%gen_ExtPin6Panel_20261026.m -- Ben's overruling (2026-10-26, evening):
% the SIX-test pool stands (training {43.6, 43.2-tendon, 45.7}, validation
% {40.0, 40.0-tendon, 46.5}); this restores the six-panel figure layout of
% gen_ExtPin6Panel.m (panels A-C training incl. 45.7, D-F validation,
% legend on panel C) at the adopted row-32 terms, writing to the dated
% asset name ExtPin_group_20261026.eps (same file Main.tex references).
% Torque label M with bold-k knee subscript per the notation ruling.
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

%% adopted extensor terms of record (canonical mat, row 32, ledger-verified)
S = load('minimizeExt10mmX3_results.mat');
g = [];
if isfield(S, 'x_actual')
    g = S.x_actual(32, 1:4);
elseif isfield(S, 'results_sort_actual')
    g = S.results_sort_actual(32, 4:7);
end
ledger = [-6.378022e-03, 3.997873e+04, 1.473375e+04, 1.582887e-01];
if isempty(g) || any(abs(g(1:4) - ledger) > 1e-6 * abs(ledger))
    error('Row 32 mismatch vs ledger: got [%s]', sprintf(' %.6e', g));
end
fprintf('ADOPTED row32: X0=%.6e m  X1=%.6e  X2=%.6e N/m  X3=%.6f\n', g(1), g(2), g(3), g(4));

%% evaluator calls (tests arg unreliable; call 1:9 and subset)
tests = 1:9;
[f0, b0] = minimizeExtX3(0, Inf, Inf, 0, tests);
[f1, b1] = minimizeExtX3(g(1), g(2), g(3), g(4), tests);
pool = [1 2 5 6 7 8];   % globals: 40.0, 40.0t, 43.6, 43.2t, 45.7, 46.5
trG = [5 6 7]; vaG = [1 2 8];
fprintf('MEANS training {43.6,43.2t,45.7}: unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(trG,:),1), mean(f1(trG,:),1));
fprintf('MEANS validation {40,40t,46.5}:    unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(vaG,:),1), mean(f1(vaG,:),1));

%% 3x2 figure: 43.6, 43.2t, 45.7 training; 40.0, 40.0t, 46.5 validation
pord = [5 6 7 1 2 8];
pttl = {'44 cm', '43 cm Tendon', '46 cm', '40 cm', '40 cm Tendon', '47 cm'};
psub = {'Training','Training','Training','Validation','Validation','Validation'};
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

fig = figure('Name','ExtPin_6panel_20261026','Color','w'); fig.Position = [60 60 1150 780];
t = tiledlayout(2,3,'TileSpacing','loose','Padding','loose');
allY = [];
axT = gobjects(1,6);
for k = 1:6
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
    if k == 3, lg = legend('Location','northeast'); lg.FontSize = 8; lgAx = gca; end
end
for k = 1:6
    axT(k).YLim = [0, ceil(max(allY)/5)*5];
end
ylabel(t, '\bf M_{\bf k}, N\cdotm'); xlabel(t, '\bf \theta_{\bf k} , \circ');
ann = {'(A)','(B)','(C)','(D)','(E)','(F)'};
drawnow;
ylL = 0.008; ylR = 0.024;
xL1 = 0.040;
if (xL1 < ylR) && (xL1 + 0.032 > ylL)
    Sx = (xL1 + 0.032 + 0.008) - ylL;
    tp = t.Position; t.Position = [tp(1) + Sx, tp(2), tp(3) - Sx, tp(4)];
    drawnow;
    ylL = ylL + Sx; ylR = ylR + Sx;
end
for j = 1:6
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
    fprintf('AUDIT tile %d %s: letter@(%s)\n', j, ann{j}, mat2str(a.Position(1:2)));
end
lgAll = findall(fig, 'Type', 'legend');
fprintf('AUDIT legend: count=%d  location=%s  on tile-3 axes=%d\n', ...
    numel(lgAll), lgAll(1).Location, lgAx == axT(3));
exportgraphics(fig, fullfile(outdir,'ExtPin_group_20261026.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(qcdir,'QC_ExtPin_6panel_20261026.png'), 'Resolution',150);
fprintf('WROTE ExtPin_group_20261026.eps (6-panel restoration)\n'); close(fig);
disp('DONE_ExtPin6Panel_20261026');
