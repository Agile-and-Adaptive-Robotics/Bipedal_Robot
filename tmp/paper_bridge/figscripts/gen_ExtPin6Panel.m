%gen_ExtPin6Panel.m — Ben's 2026-10-18 correction: the 3x2 six-panel figure is
% the PINNED EXTENSOR group (his online fig 5), showing the six pool tests of
% the nine (training {43.6, 43.2-tendon, 45.7}, validation {40.0, 40.0-tendon,
% 46.5}; 42, 42-tendon, 48 not used) at the adopted extensor terms (row 32).
% Series per panel: Hybrid gold, Measured blue, Unoptimized light-orange dashed,
% Optimized magenta. ONE legend (panel A). Replaces ExtPin_group.eps.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');

c1 = [1 0.843 0];         % gold        - Hybrid
c3 = [0.980 0.529 0.459]; % light orange - Unoptimized
c5 = [0.804 0.204 0.710]; % magenta     - Optimized
c7 = [0 0 1];             % blue        - Measured
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

%% evaluator calls — NOTE: the tests argument does NOT subset reliably (a
%% [1 2 5 6 7 8] call returned globals 1-6); call 1:9 and subset by index.
tests = 1:9;
[f0, b0] = minimizeExtX3(0, Inf, Inf, 0, tests);          % zero-correction
[f1, b1] = minimizeExtX3(g(1), g(2), g(3), g(4), tests);  % adopted terms
pool = [1 2 5 6 7 8];   % globals: 40.0, 40.0t, 43.6, 43.2t, 45.7, 46.5
lab = {'40.0 cm (val)', '40.0 cm tendon (val)', '43.6 cm (tr)', '43.2 cm tendon (tr)', '45.7 cm (tr)', '46.5 cm (val)'};
fprintf('\n=== Extensor pool GoF unoptimized [RMSE FVU Max] ===\n');
for m = 1:6, i = pool(m); fprintf('%-24s %8.4f %8.4f %8.4f  rest=%.4f m\n', lab{m}, f0(i,:), b0(i).rest); end
fprintf('=== Extensor pool GoF optimized [RMSE FVU Max] ===\n');
for m = 1:6, i = pool(m); fprintf('%-24s %8.4f %8.4f %8.4f\n', lab{m}, f1(i,:)); end
trG = [5 6 7]; vaG = [1 2 8];   % global test indices
fprintf('MEANS training {43.6,43.2t,45.7}: unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(trG,:),1), mean(f1(trG,:),1));
fprintf('MEANS validation {40,40t,46.5}:    unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(vaG,:),1), mean(f1(vaG,:),1));

%% 3x2 figure (panels: 43.6, 43.2t, 45.7 training; 40.0, 40.0t, 46.5 validation)
pord = [5 6 7 1 2 8];
pttl = {'44 cm', '43 cm Tendon', '46 cm', '40 cm', '40 cm Tendon', '47 cm'};
psub = {'Training','Training','Training','Validation','Validation','Validation'};
sty = {'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 30],'XTick',[-120 -90 -60 -30 0 30],'GridLineStyle','none'};

% x upper limit from the data: 30 unless a pool test's measured points run past it
xmax = 30;
for i = pool
    xmax = max(xmax, max(b1(i).Aexp(b1(i).Aexp < 180)));  % guard vs bad samples
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

fig = figure('Name','ExtPin_6panel','Color','w'); fig.Position = [60 60 1150 780];
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
    if k == 3, lg = legend('Location','northeast'); lg.FontSize = 8; end
end
for ax = reshape(findall(fig, 'Type', 'axes'), 1, [])
    ax.YLim = [0, ceil(max(allY)/5)*5];
end
ylabel(t, '\bf Torque, N\cdotm'); xlabel(t, '\bf \theta_{k} , \circ');
% panel letters anchored inside each tile's actual axes box
ann = {'(A)','(B)','(C)','(D)','(E)','(F)'};
drawnow;
for j = 1:6
    pos = axT(j).Position;
    annotation(fig,'textbox',[pos(1) + 0.004, pos(2) + pos(4) - 0.058, 0.035 0.05], ...
        'String',['\bf ' ann{j}], 'FontSize',12,'FontName','Arial', ...
        'EdgeColor','none','HorizontalAlignment','left');
end
exportgraphics(fig, fullfile(outdir,'ExtPin_group.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(qcdir,'QC_ExtPin_6panel.png'), 'Resolution',150);
fprintf('WROTE ExtPin_group.eps (6-panel replacement)\n'); close(fig);
disp('DONE_ExtPin6Panel');
