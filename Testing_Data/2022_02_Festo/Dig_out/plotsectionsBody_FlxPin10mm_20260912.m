%% --- Define color scheme and labels ---
c = cell(8,1);
c{1} = '#FFD700'; % gold -> Hybrid
c{2} = '#FFB14E'; % orange
c{3} = '#FA8775'; % light orange
c{4} = '#EA5F94'; % pink
c{5} = '#CD34B5'; % magenta -> Predicted
c{6} = '#9D02D7'; % magenta 2
c{7} = '#0000FF'; % indigo -> Measured
c{8} = '#000000'; % black

% allBPA = allBPA;    %Plot only the train and validation BPAs
allBPA = [1, 2, 3, 4, 5]; %Plot all tests
numBPA = numel(allBPA); %recalculate if allBPA has changed

% Auto panel letters (A), (B), ... one per plotted test, journal caption style
tileLabels = arrayfun(@(k) sprintf('(%c)', 'A' + k - 1), 1:numBPA, 'UniformOutput', false);
% Annotation positions [x, y] in normalized figure units
xAnn = [0, 0.48, 0, 0.48];
yAnn = [0.94, 0.94, 0.45, 0.45];
sz = 60;

%% Plot torque curves, pre- and post-Optimized
load ForceStrainForFit.mat z
z(isnan(z)) = 0;
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

figTpre = figure('Name','Torque, Pre-Optimized','Color','w');
figTpre.Position = [100 100 950 700];
tTpre = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

titles = ["\bf 48.5 cm", "\bf 45.7 cm","\bf 47.9 cm", "\bf 40.6 cm", "\bf 41.7 cm"]; %indexed by test number j (allBPA may plot a subset)
% titles = validLabels;
subtitles = ["\bf Pre-optimized","\bf Pre-optimized","\bf Pre-optimized","\bf Pre-optimized","\bf Optimized","\bf Optimized","\bf Optimized","\bf Optimized"];

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    hold on
    % Pre-optimization: Mold calculation
    Yq = bpa(j).strain./((bpa(j).rest - bpa(j).Kmax)/bpa(j).rest);
    Xq = bpa(j).P.*ones(size(Yq));
    Vq = interp2(X, Y', z, Xq, Yq,'spline');
    Fold = Vq.*bpa(j).unitD;     %Old force calc              
    Fq = hypot(Fold(:,1), Fold(:,2)); %Old force magnitude on the XY plane
    Mold = -bpa(j).mA.*Fq;       %Torque with old force calc
    Mold(Yq>1) = 0;
    Mold(bpa(j).strain < 0) = NaN;
    scatter(bpa(j).A_h, bpa(j).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{1}, 'DisplayName', 'Hybrid');
    scatter(bpa(j).Aexp, bpa(j).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Mold, '--', 'Color', c{8}, 'LineWidth', 2, 'DisplayName', 'Original');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');

    clear Vq Fold Fq Mold
    title(titles(j), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    text(ax, 0.002, 1.03, tileLabels{k}, 'Units', 'normalized', 'Clipping', 'off', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'VerticalAlignment', 'bottom');
    
    % ylabel('\bf Torque, N \cdot m', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    % xlabel('\bf \theta_{k}, \circ', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    % subtitle(subtitles(k), 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlim([-120 20]); ylim([-25 0]);
end

ylabel(tTpre,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpre,'\bf \theta_{k} , \circ','Interpreter','tex');


figTpost = figure('Name','Torque, Post-Optimized','Color','w');
figTpost.Position = [100 100 950 700];
tTpost = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    hold on
    
    scatter(bpa(j).Aexp, bpa(j).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).M_p(:,3), '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    
    title(titles(j), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    text(ax, 0.002, 1.03, tileLabels{k}, 'Units', 'normalized', 'Clipping', 'off', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'VerticalAlignment', 'bottom');

    % ylabel('\bf Torque, N \cdot m', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    % xlabel('\bf \theta_{k}, \circ', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    % subtitle(subtitles(k), 'FontSize', 10, 'FontName', 'Arial', 'FontWeight', 'bold');
    xlim([-120 20]); ylim([-25 0]);
end

% Shared labels
ylabel(tTpost,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpost,'\bf \theta_{k} , \circ','Interpreter','tex');

% (A)-(D) annotations
% for j = 1:4
%     annotation(gcf, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end

% Column titles
% annotation(gcf, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(gcf, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% Legend placement: even number of tests -> legend inside tile (1,2);
% odd -> the first empty tile of the (ceil(n/2),2) grid (e.g. tile 6 for 5 tests)
if mod(numBPA,2) == 0
    lg = legend(tTpre.Children(numBPA-1)); %tile (1,2) series
else
    lg = legend(tTpre.Children(end));      %tile 1 series, moved to the empty tile
    lg.Layout.Tile = 2*ceil(numBPA/2);
end
lg.Location = 'best';
lg.FontSize = 8;

if mod(numBPA,2) == 0
    lg2 = legend(tTpost.Children(numBPA-1)); %tile (1,2) series
else
    lg2 = legend(tTpost.Children(end));      %tile 1 series, moved to the empty tile
    lg2.Layout.Tile = 2*ceil(numBPA/2);
end
lg2.Location = 'best';
lg2.FontSize = 8;

%% Plot muscle length, optimization and validation
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    text(ax, 0.002, 1.03, tileLabels{k}, 'Units', 'normalized', 'Clipping', 'off', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'VerticalAlignment', 'bottom');
    hold on
    
    % Calculate predicted
    Lm_p = bpa(j).Lmt_p - 2 * bpa(j).fitn - bpa(j).ten;
    Lm   = bpa(j).Lmt   - 2 * bpa(j).fitn - bpa(j).ten;
    % Measured
    scatter(bpa(j).A_h, bpa(j).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    %Old prediction
    plot(bpa(j).Ak, Lm, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Original prediction');
    % New prediction
    plot(bpa(j).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');


    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tL,'\bf Muscle length, m','Interpreter','tex');
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex');


% for j = 1:4
%     annotation(figL, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end
% 
% annotation(figL, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(figL, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

if mod(numBPA,2) == 0
    lg = legend(tL.Children(numBPA-1)); %tile (1,2) series
else
    lg = legend(tL.Children(end));      %tile 1 series, moved to the empty tile
    lg.Layout.Tile = 2*ceil(numBPA/2);
end
lg.Location = 'best';
lg.FontSize = 8;

% lg2 = legend(tL.Children(end-(numBPA-1)));
% lg2.Location = 'best';
% lg2.FontSize = 8;

%% Plot moment arm, optimization and validation
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');


for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    text(ax, 0.002, 1.03, tileLabels{k}, 'Units', 'normalized', 'Clipping', 'off', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'VerticalAlignment', 'bottom');
    hold on
    
    G_p = hypot(bpa(j).mA_p(:,1), bpa(j).mA_p(:,2));
    
    scatter(bpa(j).A_h, bpa(j).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).mA, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tMA,'\bf Moment arm, m','Interpreter','tex');
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex');
% for j = 1:4
%     annotation(figMA, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end
% annotation(figMA, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(figMA, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
if mod(numBPA,2) == 0
    lg = legend(tMA.Children(numBPA-1)); %tile (1,2) series
else
    lg = legend(tMA.Children(end));      %tile 1 series, moved to the empty tile
    lg.Layout.Tile = 2*ceil(numBPA/2);
end
lg.Location = 'best';
lg.FontSize = 8;
% legend(tMA.Children(end-(numBPA-1)),'Location','best','FontSize',8);

%% Plot relative strain, optimization and validation
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');

for k = 1:numBPA
    ax = nexttile(k);
    j = allBPA(k);
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    text(ax, 0.002, 1.03, tileLabels{k}, 'Units', 'normalized', 'Clipping', 'off', ...
        'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial', 'VerticalAlignment', 'bottom');
    hold on
    
    strain_h = (bpa(j).rest - bpa(j).Lm_h)/bpa(j).rest;
    kmax = (bpa(j).rest - bpa(j).Kmax)/bpa(j).rest;
    scatter(bpa(j).A_h, strain_h/kmax, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).strain/kmax, '-.', 'Color', c{3}, 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).strain_p/kmax, '-', 'Color', c{5}, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

    % if j == 1
    %     ylabel('\bf Optimized', 'Interpreter', 'tex', ...
    %         'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    % end
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end

ylabel(tS,'\bf Relative strain','Interpreter','tex');
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex');
% for j = 1:4
%     annotation(figS, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
%         'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% end
% annotation(figS, 'textbox', [0.2, 0.95, 0.1, 0.05], 'String', '\bf Optimization', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% annotation(figS, 'textbox', [0.7, 0.95, 0.1, 0.05], 'String', '\bf Validation', ...
%     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
if mod(numBPA,2) == 0
    lg = legend(tS.Children(numBPA-1)); %tile (1,2) series
else
    lg = legend(tS.Children(end));      %tile 1 series, moved to the empty tile
    lg.Layout.Tile = 2*ceil(numBPA/2);
end
lg.Location = 'best';
lg.FontSize = 8;
% legend(tS.Children(end-(numBPA+1)),'Location','best','FontSize',8);


