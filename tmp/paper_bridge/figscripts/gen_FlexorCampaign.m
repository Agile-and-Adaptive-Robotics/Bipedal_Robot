%gen_FlexorCampaign.m — Ben's 2026-10-18 results-figure directive.
% Produces: (P1) FlxPin_preopt.eps  pre-optimized 2x2 (Hunt Eq, Unoptimized,
%               Measured, Hybrid) over the four cross-validation fold tests;
%           (P2) FlxPin_optimized.eps  optimized 2x2 (Unoptimized dashed,
%               Optimized, Measured filled) with Training/Validation subtitles;
%           (P3) FlxPin_group.eps (REPLACES online fig 5) 3x2: all five pinned
%               tests with Hybrid/Measured/Unoptimized/Optimized + legend cell;
%           (P4) FlxGrp.eps (REPLACES) biomimetic 1x2 at the adopted Xi terms.
% Fold of record (chapter tab:FlxPinGof caption): training {48.5, 41.0} cm,
% validation {45.7, 40.6-tendon} cm, 47.9 cm not used.
% Style: PlotPrompt palette + Plot_KneeFlx_all conventions (Arial 12 bold,
% LineWidth 2, minor ticks, legend 8 pt, XLim [-120 15]). All GoF printed.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');   % default off

% palette (PlotPrompt c{1..8})
c1 = [1 0.843 0];         % gold        - Hybrid
c3 = [0.980 0.529 0.459]; % light orange - Unoptimized (improved BPA char alone)
c4 = [0.918 0.373 0.580]; % pink        - reduced-pressure points
c5 = [0.804 0.204 0.710]; % magenta     - Optimized
c7 = [0 0 1];             % blue        - Measured (filled, no edge)
c8 = [0 0 0];             % black       - human model
cgray = [0.4 0.4 0.4];    % gray dashed - Hunt equation (older force surface)
sz = 60;

%% ---- 1. adopted terms of record (canonical mat, row 77, ledger-verified) ----
S = load('minimizeFlxPin10_results.mat');
if isfield(S, 'filtered_results') && isfield(S, 'xCols')
    g = S.filtered_results(77, S.xCols);
elseif isfield(S, 'results_sort_actual')
    g = S.results_sort_actual(77, 4:6);
else
    error('Cannot locate adopted flexor terms in minimizeFlxPin10_results.mat');
end
ledger = [3.937370e-03, 3.997873e+04, 1.473375e+04];
if any(abs(g(1:3) - ledger) > 1e-6 * ledger)
    error('Row 77 mismatch vs ledger: got [%.6e %.6e %.6e]', g(1), g(2), g(3));
end
fprintf('ADOPTED row77: X0=%.6e m  X1=%.6e N/m  X2=%.6e N/m\n', g(1), g(2), g(3));

%% ---- 2. evaluator calls over all five pinned tests ----
tests5 = [1 2 3 4 5];
labels5 = {'48 cm (48.5)', '46 cm (45.7)', '47 cm (47.9)', '40 cm Tendon (40.6)', '41 cm (41.0)'};
[f0, bpa0] = minimizeFlxPin(0, Inf, Inf, tests5);   % zero-correction = Unoptimized
[f1, bpa1] = minimizeFlxPin(g(1), g(2), g(3), tests5); % optimized terms
fprintf('\n=== GoF zero-correction (Unoptimized) [RMSE FVU Max] ===\n');
for i = 1:5, fprintf('%-22s %8.4f %8.4f %8.4f\n', labels5{i}, f0(i,:)); end
fprintf('=== GoF corrected (Optimized) [RMSE FVU Max] ===\n');
for i = 1:5, fprintf('%-22s %8.4f %8.4f %8.4f\n', labels5{i}, f1(i,:)); end
tr = [1 5]; va = [2 4];
fprintf('MEANS training {48.5,41.0}: unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(tr,:),1), mean(f1(tr,:),1));
fprintf('MEANS validation {45.7,40.6t}: unopt %.4f/%.4f/%.4f  opt %.4f/%.4f/%.4f\n', ...
    mean(f0(va,:),1), mean(f1(va,:),1));

%% ---- 3. Hunt-equation (older force surface) series + GoF ----
load('ForceStrainForFit.mat', 'z');
z(isnan(z)) = 0;
Xp = linspace(0, 620, 20); Xp = Xp(2:20);
Yp = linspace(0, 1, 30);
Mold = cell(5, 1); gofH = zeros(5, 3);
for i = 1:5
    Yq = bpa0(i).strain ./ ((bpa0(i).rest - bpa0(i).Kmax) / bpa0(i).rest);
    Xq = bpa0(i).P .* ones(size(Yq));
    Vq = interp2(Xp, Yp', z, Xq, Yq, 'spline');
    Fold = Vq .* bpa0(i).unitD;
    Fq = hypot(Fold(:,1), Fold(:,2));
    Mo = -bpa0(i).mA .* Fq;
    Mo(Yq > 1) = 0;
    Mo(bpa0(i).strain < 0) = NaN;
    Mold{i} = Mo;
    idx = round(interp1(bpa0(i).Ak, 1:numel(bpa0(i).Ak), bpa0(i).Aexp, 'nearest', 'extrap'));
    Mm = Mo(idx); Me = bpa0(i).Mexp;
    v = ~isnan(Mm) & isfinite(Me);
    r = Mm(v) - Me(v);
    gofH(i,:) = [sqrt(mean(r.^2)), sum(r.^2)/sum((Me(v)-mean(Me(v))).^2), max(abs(r))];
end
fprintf('=== GoF Hunt equation [RMSE FVU Max] ===\n');
for i = 1:5, fprintf('%-22s %8.4f %8.4f %8.4f\n', labels5{i}, gofH(i,:)); end
fprintf('MEANS training: %.4f/%.4f/%.4f  validation: %.4f/%.4f/%.4f\n', ...
    mean(gofH(tr,:),1), mean(gofH(va,:),1));

%% ---- FIG P1: pre-optimized (Hunt, Unoptimized, Measured, Hybrid) ----
sty = {'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 15],'XTick',[-120 -90 -60 -30 0],'GridLineStyle','none'};
p1order = [1 5 2 4];   % 48.5 T, 41.0 T, 45.7 V, 40.6t V
p1sub = {'Training','Training','Validation','Validation'};
p1ttl = {'48 cm','41 cm','46 cm','40 cm Tendon'};
fig1 = figure('Name','FlxPin_preopt','Color','w'); fig1.Position = [80 80 950 700];
t1 = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
allY1 = [];
for k = 1:4
    i = p1order(k);
    nexttile(k); hold on
    plot(bpa0(i).Ak, Mold{i}, '--', 'Color', cgray, 'LineWidth', 2, 'DisplayName','Hunt equation');
    plot(bpa0(i).Ak, bpa0(i).M, '-.', 'Color', c3, 'LineWidth', 2.5, 'DisplayName','Unoptimized');
    scatter(bpa0(i).Aexp, bpa0(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
    scatter(bpa0(i).A_h, bpa0(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName','Hybrid');
    title(['\bf ' p1ttl{k}]); subtitle(['\bf ' p1sub{k}], 'FontName','Arial');
    set(gca, sty{:});
    allY1 = [allY1; Mold{i}(~isnan(Mold{i})); bpa0(i).M; bpa0(i).Mexp; bpa0(i).M_h]; %#ok<AGROW>
    if k == 4, lg = legend('Location','southwest'); lg.FontSize = 8; end
end
ylabel(t1, '\bf Torque, N\cdotm'); xlabel(t1, '\bf \theta_{k} , \circ');
set(t1, 'TileSpacing','loose');
ann1 = {'(A)','(B)','(C)','(D)'}; xp1 = [0.001 0.48 0.001 0.48]; yp1 = [0.935 0.935 0.45 0.45];
for j = 1:4
    annotation(fig1,'textbox',[xp1(j) yp1(j) 0.05 0.05],'String',['\bf ' ann1{j}], ...
        'FontSize',12,'FontName','Arial','EdgeColor','none');
end
for ax = reshape(findall(fig1, 'Type', 'axes'), 1, [])
    ax.YLim = [floor(min(allY1)/5)*5, 0];
end
exportgraphics(fig1, fullfile(outdir,'FlxPin_preopt.eps'), 'ContentType','vector');
exportgraphics(fig1, fullfile(qcdir,'QC_FlxPin_preopt.png'), 'Resolution',150);
fprintf('WROTE FlxPin_preopt.eps\n'); close(fig1);

%% ---- FIG P2: optimized (Unoptimized dashed, Optimized, Measured filled) ----
fig2 = figure('Name','FlxPin_optimized','Color','w'); fig2.Position = [80 80 950 700];
t2 = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
allY2 = [];
for k = 1:4
    i = p1order(k);
    nexttile(k); hold on
    plot(bpa1(i).Ak, bpa1(i).M, '--', 'Color', c3, 'LineWidth', 2, 'DisplayName','Unoptimized');
    plot(bpa1(i).Ak, bpa1(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Optimized');
    scatter(bpa1(i).Aexp, bpa1(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
    title(['\bf ' p1ttl{k}]); subtitle(['\bf ' p1sub{k}], 'FontName','Arial');
    set(gca, sty{:});
    allY2 = [allY2; bpa1(i).M; bpa1(i).M_p(:,3); bpa1(i).Mexp]; %#ok<AGROW>
    if k == 4, lg = legend('Location','southwest'); lg.FontSize = 8; end
end
ylabel(t2, '\bf Torque, N\cdotm'); xlabel(t2, '\bf \theta_{k} , \circ');
for j = 1:4
    annotation(fig2,'textbox',[xp1(j) yp1(j) 0.05 0.05],'String',['\bf ' ann1{j}], ...
        'FontSize',12,'FontName','Arial','EdgeColor','none');
end
for ax = reshape(findall(fig2, 'Type', 'axes'), 1, [])
    ax.YLim = [floor(min(allY2)/5)*5, 0];
end
exportgraphics(fig2, fullfile(outdir,'FlxPin_optimized.eps'), 'ContentType','vector');
exportgraphics(fig2, fullfile(qcdir,'QC_FlxPin_optimized.png'), 'Resolution',150);
fprintf('WROTE FlxPin_optimized.eps\n'); close(fig2);

%% ---- FIG P3: replacement for online fig 5 (3x2, all five tests + legend cell) ----
p3order = [1 5 2 4 3];  % 48.5 T, 41.0 T, 45.7 V, 40.6t V, 47.9 N
p3sub = {'Training','Training','Validation','Validation','Not used'};
p3ttl = {'48 cm','41 cm','46 cm','40 cm Tendon','47 cm'};
fig3 = figure('Name','FlxPin_group','Color','w'); fig3.Position = [60 60 1150 780];
t3 = tiledlayout(2,3,'TileSpacing','loose','Padding','loose');
allY3 = [];
for k = 1:5
    i = p3order(k);
    nexttile(k); hold on
    scatter(bpa1(i).A_h, bpa1(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName','Hybrid');
    scatter(bpa1(i).Aexp, bpa1(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
    plot(bpa1(i).Ak, bpa1(i).M, '--', 'Color', c3, 'LineWidth', 2, 'DisplayName','Unoptimized');
    plot(bpa1(i).Ak, bpa1(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Optimized');
    title(['\bf ' p3ttl{k}]); subtitle(['\bf ' p3sub{k}], 'FontName','Arial');
    set(gca, sty{:});
    allY3 = [allY3; bpa1(i).M; bpa1(i).M_p(:,3); bpa1(i).Mexp; bpa1(i).M_h]; %#ok<AGROW>
end
% tile 6: shared legend + fold key
nexttile(6); hold on
hL(1) = scatter(NaN,NaN,sz,'filled','MarkerFaceColor',c1,'MarkerEdgeColor','none'); %#ok<SAGROW>
hL(2) = scatter(NaN,NaN,sz,'filled','MarkerFaceColor',c7,'MarkerEdgeColor','none'); %#ok<SAGROW>
hL(3) = plot(NaN,NaN,'--','Color',c3,'LineWidth',2); %#ok<SAGROW>
hL(4) = plot(NaN,NaN,'-','Color',c5,'LineWidth',2.5); %#ok<SAGROW>
lgL = legend(hL, {'Hybrid','Measured','Unoptimized','Optimized'}, 'Location','northwest');
lgL.FontSize = 11;
text(0.05, 0.35, {'Fold of record:', 'Training: 48.5, 41.0 cm', ...
    'Validation: 45.7, 40.6 cm tendon', 'Not used: 47.9 cm', ...
    '', 'Negative values denote', 'knee flexion.'}, ...
    'Units','normalized', 'FontSize', 11, 'FontWeight','bold', 'FontName','Arial', ...
    'VerticalAlignment','top');
axis off;
set(gca, sty{:});
for ax = reshape(findall(fig3, 'Type', 'axes'), 1, [])
    ax.YLim = [floor(min(allY3)/5)*5, 0];
end
ylabel(t3, '\bf Torque, N\cdotm'); xlabel(t3, '\bf \theta_{k} , \circ');
ann3 = {'(A)','(B)','(C)','(D)','(E)'};
xp3 = {[0.01 0.93],[0.345 0.93],[0.68 0.93],[0.01 0.45],[0.345 0.45]};
for j = 1:5
    annotation(fig3,'textbox',[xp3{j}(1) xp3{j}(2) 0.04 0.05],'String',['\bf ' ann3{j}], ...
        'FontSize',12,'FontName','Arial','EdgeColor','none');
end
exportgraphics(fig3, fullfile(outdir,'FlxPin_group.eps'), 'ContentType','vector');
exportgraphics(fig3, fullfile(qcdir,'QC_FlxPin_group_new.png'), 'Resolution',150);
fprintf('WROTE FlxPin_group.eps (replacement)\n'); close(fig3);

%% ---- FIG P4: FlxGrp replacement (biomimetic, adopted Xi terms) ----
[~, ~, ~, bpab] = minimizeFlx(g(1), g(2), g(3));
fprintf('\n=== biomimetic BPA values (non-stiffness-aware dump at adopted terms) ===\n');
for i = 1:numel(bpab)
    fprintf('bpa(%d): rest=%.4f m  Kmax=%.4f  tendon=%.4f m  P=%.1f kPa\n', ...
        i, bpab(i).rest, bpab(i).Kmax, bpab(i).ten, mean(bpab(i).P));
end
Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)'; Bifemsh_T = Tab(:,4)';

fig4 = figure('Name','FlxGrp','Color','w'); fig4.Position = [90 90 950 480];
t4 = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');
% Tile A: 10 mm
nexttile(1); hold on
i = 1;
plot(bpab(i).Ak, bpab(i).M, '--', 'Color', c3, 'LineWidth', 2, 'DisplayName','Unoptimized');
plot(bpab(i).Ak, bpab(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Optimized');
scatter(bpab(i).Aexp, bpab(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
title('\phi 10 mm'); set(gca, sty{:});
lgA = legend('Location','southwest'); lgA.FontSize = 8;
% Tile B: 20 mm at 620 kPa + human + reduced-pressure points
nexttile(2); hold on
i = 2;
plot(knee_angle_rT, Bifemsh_T, ':', 'Color', c8, 'LineWidth', 4, 'DisplayName','Human model (max)');
plot(bpab(i).Ak, bpab(i).M, '--', 'Color', c3, 'LineWidth', 2, 'DisplayName','Unoptimized');
plot(bpab(i).Ak, bpab(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Optimized');
scatter(bpab(i).Aexp, bpab(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured, 620 kPa');
% exactly four reduced-pressure measured points (recipe from minimizeFlx10mm.m)
angle_data = readmatrix('Results_table_FullSize.xlsx','Sheet','FlxTest20mm_42cm (3)','Range','C7:X7');
pressure_data = readmatrix('Results_table_FullSize.xlsx','Sheet','FlxTest20mm_42cm (3)','Range','C8:X8');
torque_data = readmatrix('Results_table_FullSize.xlsx','Sheet','FlxTest20mm_42cm (3)','Range','C17:X17');
angle281 = angle_data(pressure_data == 281);
mask = ismember(pressure_data, [560, 421, 281]);
near325 = abs(pressure_data - 325) <= 5;
idx325 = find(near325 & angle_data < angle281, 1, 'first');
if ~isempty(idx325), mask(idx325) = true; end
scatter(angle_data(mask), torque_data(mask), sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c4, 'MarkerEdgeColor', 'none', 'DisplayName','Measured < 620 kPa');
title('\phi 20 mm'); set(gca, sty{:});
lgB = legend('Location','southwest'); lgB.FontSize = 8;
ylabel(t4, '\bf Torque, N\cdotm'); xlabel(t4, '\bf \theta_{k} , \circ');
annotation(fig4,'textbox',[0.01 0.9 0.04 0.05],'String','\bf (A)','FontSize',12,'FontName','Arial','EdgeColor','none');
annotation(fig4,'textbox',[0.46 0.9 0.04 0.05],'String','\bf (B)','FontSize',12,'FontName','Arial','EdgeColor','none');
exportgraphics(fig4, fullfile(outdir,'FlxGrp.eps'), 'ContentType','vector');
exportgraphics(fig4, fullfile(qcdir,'QC_FlxGrp_new.png'), 'Resolution',150);
fprintf('WROTE FlxGrp.eps (replacement)\n'); close(fig4);

disp('DONE_FlexorCampaign');
