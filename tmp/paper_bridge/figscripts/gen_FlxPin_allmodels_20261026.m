%gen_FlxPin_allmodels_20261026.m -- P1-b rebuild of the merged flexor
% pinned-knee figure. Identical layout/series to gen_FlxPin_allmodels_
% 20261023.m EXCEPT the Hunt-equation series is now a DIRECT call of the
% equation of record (HuntEq.m: epsilon = per-point relative contraction,
% epsilon_max = 0.171 the equation's own fit basis, S = 0 isometric,
% pressure = the test pressure), replacing the legacy ForceStrainForFit z
% lookup (NaN->0 + spline, relative strain on the per-actuator KMAX basis,
% mislabeled pressure axis). Audit: lookup vs direct RMSE 5.779->5.223,
% 6.633->6.249, 5.315->6.048, 3.945->3.956 (Hunt_series_audit_20261026.csv).
% Torque axis label now M (uppercase, knee-axis subscript bf k) per the
% notation ruling. Output: Figures/FlxPin_allmodels_20261026.eps + QC png.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');

c1 = [1 0.843 0];         % gold        - Hybrid
c3 = [0.980 0.529 0.459]; % light orange - Unoptimized
c5 = [0.804 0.204 0.710]; % magenta     - Corrected (adopted terms)
c7 = [0 0 1];             % indigo      - Measured
cgray = [0.4 0.4 0.4];    % gray dashed - Hunt equation
sz = 60;

%% adopted flexor terms of record (row 77, ledger-verified)
S = load('minimizeFlxPin10_results.mat');
if isfield(S, 'filtered_results') && isfield(S, 'xCols')
    g = S.filtered_results(77, S.xCols);
else
    error('Cannot locate adopted flexor terms in minimizeFlxPin10_results.mat');
end
ledger = [3.937370e-03, 3.997873e+04, 1.473375e+04];
if any(abs(g(1:3) - ledger) > 1e-6 * ledger)
    error('Row 77 mismatch vs ledger: got [%.6e %.6e %.6e]', g(1), g(2), g(3));
end
fprintf('ADOPTED row77: X0=%.6e m  X1=%.6e N/m  X2=%.6e N/m\n', g(1), g(2), g(3));

%% evaluator calls over all five tests, index the four needed
tests5 = 1:5;   % family: 1=48.5, 2=45.7, 3=47.9, 4=40.6t, 5=41.0
[f0, bpa0] = minimizeFlxPin(0, Inf, Inf, tests5);      % zero-correction
[f1, bpa1] = minimizeFlxPin(g(1), g(2), g(3), tests5); % adopted terms
ordr = [1 5 2 4];      % family indices in panel order 48.5, 41.0, 45.7, 40.6t

%% Hunt-equation series: DIRECT call of HuntEq.m (equation of record)
Mold = cell(4, 1);
for k = 1:4
    i = ordr(k);
    b = bpa0(i);
    eps = b.strain(:);
    evalc('Bh = HuntEq(eps, 0.171, b.P, 0);');
    Fq = Bh(:, 2);
    Fq(Fq < 0) = 0;
    Mo = -b.mA .* Fq;
    Mo(b.strain > 0.171) = 0;
    Mo(b.strain < 0) = NaN;
    Mold{k} = Mo;
end

%% figure: P1 layout + corrected line
sty = {'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 15],'XTick',[-120 -90 -60 -30 0],'GridLineStyle','none'};
p1order = [1 4 2 3];   % 48.5 T, 41.0 T, 45.7 V, 40.6t V
p1sub = {'Training','Training','Validation','Validation'};
p1ttl = {'48 cm','41 cm','46 cm','40 cm Tendon'};
fig1 = figure('Name','FlxPin_allmodels_20261026','Color','w'); fig1.Position = [80 80 950 700];
t1 = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');
allY1 = []; axP1 = gobjects(1,4);
for k = 1:4
    i = p1order(k);
    nexttile(k); hold on
    scatter(bpa0(i).A_h, bpa0(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName','Hybrid');
    scatter(bpa0(i).Aexp, bpa0(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
    plot(bpa0(i).Ak, Mold{k}, '--', 'Color', cgray, 'LineWidth', 2, 'DisplayName','Hunt equation');
    plot(bpa0(i).Ak, bpa0(i).M, '-.', 'Color', c3, 'LineWidth', 2.5, 'DisplayName','Unoptimized');
    plot(bpa1(i).Ak, bpa1(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, ...
        'DisplayName','Corrected model');
    title(['\bf ' p1ttl{k}]); subtitle(['\bf ' p1sub{k}], 'FontName','Arial');
    set(gca, sty{:});
    axP1(k) = gca;
    allY1 = [allY1; Mold{k}(~isnan(Mold{k})); bpa0(i).M; bpa0(i).Mexp; bpa0(i).M_h; ...
        bpa1(i).M_p(:,3)]; %#ok<AGROW>
    if k == 4, lg = legend('Location','southwest'); lg.FontSize = 8; end
end
ylabel(t1, '\bf M_{\bf k}, N\cdotm'); xlabel(t1, '\bf \theta_{\bf k} , \circ');
set(t1, 'TileSpacing', 'loose');
ann1 = {'(A)','(B)','(C)','(D)'}; xp1 = [0.001 0.48 0.001 0.48];
drawnow;
for j = 1:4
    axj = axP1(j); p4 = axj.Position;
    tcY = 0.935;
    try
        Tj = axj.Title; Tj.Units = 'normalized';
        te = Tj.Extent;
        cand = p4(2) + (te(2) + te(4)/2) * p4(4);
        if isfinite(cand) && cand > 0.3, tcY = cand; end
    catch
    end
    a1 = annotation(fig1,'textbox',[xp1(j) tcY-0.025 0.05 0.05],'String',['\bf ' ann1{j}], ...
        'FontSize',12,'FontName','Arial','FontWeight','bold','EdgeColor','none');
    try
        drawnow; ae = a1.Extent;
        a1.Position(2) = a1.Position(2) + (tcY - (ae(2) + ae(4)/2));
    catch
    end
    fprintf('AUDIT tile %d: titleCenter=%.4f\n', j, tcY);
end
for ax = reshape(findall(fig1, 'Type', 'axes'), 1, [])
    ax.YLim = [floor(min(allY1)/5)*5, 0];
end
exportgraphics(fig1, fullfile(outdir, 'FlxPin_allmodels_20261026.eps'), 'ContentType','vector');
exportgraphics(fig1, fullfile(qcdir, 'QC_FlxPin_allmodels_20261026.png'), 'Resolution',150);
fprintf('WROTE FlxPin_allmodels_20261026.eps\n'); close(fig1);
disp('DONE_FlxPin_allmodels_20261026');
