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

% curve-identity audit: the plotted "Unoptimized" series b1(i).M must equal
% the zero-correction prediction b0(i).M_p(:,3) (verified 0.0000 on the
% flexor family; stale .M is set at load time and is terms-independent)
dm = 0;
for i = pool
    dm = max(dm, max(abs(b1(i).M(:) - b0(i).M_p(:,3))));
end
fprintf('AUDIT series: max |b1.M - b0.M_p(:,3)| over pool = %.4f (expect ~0)\n', dm);

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
    if k == 3, lg = legend('Location','northeast'); lg.FontSize = 8; lgAx = gca; end
end
for ax = reshape(findall(fig, 'Type', 'axes'), 1, [])
    ax.YLim = [0, ceil(max(allY)/5)*5];
end
ylabel(t, '\bf Torque, N\cdotm'); xlabel(t, '\bf \theta_{k} , \circ');
% panel letters per Ben's FlxPin_group.fig convention: margin labels LEFT of
% each tile's outer box; VERTICALLY CENTERED ON THE PANEL TITLE (Ben
% 2026-10-19: letters at the same y level as the figure title). (A)/(D) sit a
% further ~14 thou -x with the left margin opened if the y-label is in the way.
ann = {'(A)','(B)','(C)','(D)','(E)','(F)'};
drawnow;
% rotated y-label band, measured off the render (layout Text exposes no
% Extent/Units): x ~ [0.008, 0.024] figure-normalized
ylL = 0.008; ylR = 0.024;
xL1 = 0.040;                                   % (A)/(D) target left edge
if (xL1 < ylR) && (xL1 + 0.032 > ylL)          % y-label overlaps: open the margin
    S = (xL1 + 0.032 + 0.008) - ylL;
    tp = t.Position; t.Position = [tp(1) + S, tp(2), tp(3) - S, tp(4)];
    drawnow;
    ylL = ylL + S; ylR = ylR + S;
end
fprintf('AUDIT margin: yLabel band x[%.4f %.4f]  (A)/(D) left=%.4f  gap=%.4f\n', ...
    ylL, ylR, xL1, xL1 - ylR);
apos = zeros(6,2); pos0 = zeros(6,4); outL = zeros(6,1); ttlC = zeros(6,1); ltrC = nan(6,1);
for j = 1:6
    pos0(j,:) = axT(j).Position;
    outL(j) = axT(j).OuterPosition(1);
    innerTop = pos0(j,2) + pos0(j,4);
    tcY = innerTop + 0.030;                    % fallback if title extent unreadable
    try
        Tj = axT(j).Title; Tj.Units = 'normalized';
        te = Tj.Extent;                        % axes-inner-box normalized (y > 1 above the box)
        cand = pos0(j,2) + (te(2) + te(4)/2) * pos0(j,4);
        if isfinite(cand) && cand > innerTop && cand < innerTop + 0.09
            tcY = cand;
        end
    catch
    end
    ttlC(j) = tcY;
    if ismember(j, [1 4]), xj = xL1; else, xj = outL(j) - 0.037; end
    a = annotation(fig,'textbox',[xj, tcY - 0.020, 0.032, 0.04], ...
        'String',['\bf ' ann{j}], 'FontSize',12,'FontName','Arial', ...
        'FontWeight','bold','EdgeColor','none','HorizontalAlignment','left', ...
        'VerticalAlignment','top');
    try
        drawnow;
        ae = a.Extent;                         % exact centering onto the title line
        a.Position(2) = a.Position(2) + (tcY - (ae(2) + ae(4)/2));
        ltrC(j) = ae(2) + ae(4)/2;
    catch
    end
    apos(j,:) = a.Position(1:2);
end

% audit: print the final drawn state (Ben's fixes, verified not assumed)
lgAll = findall(fig, 'Type', 'legend');
fprintf('AUDIT legend: count=%d  location=%s  on tile-3 axes=%d\n', ...
    numel(lgAll), lgAll(1).Location, lgAx == axT(3));
for j = 1:6
    innerTop = pos0(j,2) + pos0(j,4);
    if ismember(j, [1 4])
        anchored = abs(apos(j,1) - xL1) < 1e-9;          % (A)/(D) at the opened margin
    else
        anchored = abs(apos(j,1) + 0.032 - (outL(j) - 0.005)) < 5e-3;  % others end just before tile
    end
    titleAligned = abs(ltrC(j) - ttlC(j)) < 5e-3;        % letter centered on the title line
    fprintf('AUDIT tile %d %s: XLim=[%g %g] letter@(%.3f,%.3f) letterR=%.3f outerL=%.3f titleC=%.4f letterC=%.4f anchored=%d titleAligned=%d\n', ...
        j, ann{j}, axT(j).XLim(1), axT(j).XLim(2), apos(j,1), apos(j,2), ...
        apos(j,1) + 0.032, outL(j), ttlC(j), ltrC(j), anchored, titleAligned);
    % letter box in data coordinates of THIS tile: any plotted point under it?
    bx = axT(j).XLim(1) + ([apos(j,1), apos(j,1)+0.032] - pos0(j,1)) / pos0(j,3) * diff(axT(j).XLim);
    by = axT(j).YLim(1) + ([apos(j,2), apos(j,2)+0.04] - pos0(j,2)) / pos0(j,4) * diff(axT(j).YLim);
    i = pord(j);
    xs = [b1(i).A_h(:); b1(i).Aexp(:); b1(i).Ak(:); b1(i).Ak(:)];
    ys = [b1(i).M_h(:); b1(i).Mexp(:); b1(i).M(:); b1(i).M_p(:,3)];
    n_ov = sum(xs > bx(1) & xs < bx(2) & ys > by(1) & ys < by(2));
    fprintf('AUDIT tile %d letter data-box x[%.0f %.0f] y[%.1f %.1f] points-under-letter=%d\n', ...
        j, bx(1), bx(2), by(1), by(2), n_ov);
end
exportgraphics(fig, fullfile(outdir,'ExtPin_group.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(qcdir,'QC_ExtPin_6panel.png'), 'Resolution',150);
fprintf('WROTE ExtPin_group.eps (6-panel replacement)\n'); close(fig);
disp('DONE_ExtPin6Panel');
