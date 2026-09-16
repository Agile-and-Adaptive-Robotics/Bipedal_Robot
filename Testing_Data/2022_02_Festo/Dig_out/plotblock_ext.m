%% --- Define color scheme and labels ---
c = cell(8,1);
c{1} = '#FFD700'; % gold → Hybrid
c{2} = '#FFB14E'; % orange
c{3} = '#FA8775'; % light orange
c{4} = '#EA5F94'; % pink
c{5} = '#CD34B5'; % magenta → Predicted
c{6} = '#9D02D7'; % magenta 2
c{7} = '#0000FF'; % indigo → Measured
c{8} = '#000000'; % black

% tileIdxs = [1, 5, 8, 12, 15, 19, 22, 26, 29]; 
tileIdxs = [1, 5, 8, 12, 15];  % A, B, C, D
% tileIdxs = [1, 5, 10];  % A, B, C
el = numel(tileIdxs);
tileSpans = [1 3];      % Span: [rows cols]
% tileSpans = [1 1];      % Span: [rows cols]
% tileOrder = [1, 2, 3, 4, 5, 6, 7, 8, 9];
tileOrder = allBPA;
% tileLabels = {'(A)', '(B)', '(C)', '(D)'};
% Annotation positions [x, y] in normalized figure units
% xAnn = [0.035, 0.51, 0.035, 0.51];  % (A), (B), (C), (D)
% yAnn = [0.89, 0.89, 0.41, 0.41];    % (A), (B), (C), (D)
% xAnn = [0.035, 0.51, 0.265];  % (A), (B), (C)
% yAnn = [0.89, 0.89, 0.41];    % (A), (B), (C)
sz = 60;

%for plotting 

%% --- Torque Figure with tiles ---
figT = figure('Name','Torque','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    % Plot order: hybrid (gold), experimental (indigo), original (gray dash), predicted (magenta)
    scatter(bpa(i).A_h, bpa(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{1}, 'DisplayName', 'Hybrid');
    scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
    plot(bpa(i).Ak, bpa(i).M, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, ...
        'DisplayName', 'Original');
    plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c{5}, 'LineWidth', 2.5, ...
        'DisplayName', 'Predicted');

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % ylabel('\bf Torque, N\cdotm','Interpreter','tex')
    % xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    % annotation(figT, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'YLim', [0 15], ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');
end

%shared axes labels
ylabel(tT,'\bf Torque, N\cdotm','Interpreter','tex')
xlabel(tT,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tT.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Muscle Length Figure with tiles---
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    % Predicted
    Lm_p = bpa(i).Lmt_p - 2 * bpa(i).fitn - bpa(i).ten;
    Lm   = bpa(i).Lmt   - 2 * bpa(i).fitn - bpa(i).ten;

    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured'); % Hybrid (gold)
    plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');      % Original
    plot(bpa(i).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');            % Predicted
    % ylabel('\bf Length','Interpreter','tex')
    % xlabel('\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % annotation(figL, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');

end

%shared axes labels
ylabel(tL,'\bf Length','Interpreter','tex')
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tL.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Moment Arm Figure with tiles ---
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));

    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured');  % Hybrid
    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');      % Original
    plot(bpa(i).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');                   % Predicted

    % ylabel(ax,'\bf Moment arm, m','Interpreter','tex')
    % xlabel(ax,'\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % % annotation(figMA, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    % %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');
end

%shared axes labels
ylabel(tMA,'\bf Moment arm, m','Interpreter','tex')
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tMA.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;

%% --- Strain Figure with tiles ---
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = tileOrder(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on
    
    strain_h = (bpa(i).rest - bpa(i).Lm_h)/bpa(i).rest;
    kmax = (bpa(i).rest - bpa(i).Kmax)/bpa(i).rest;
    scatter(bpa(i).A_h, strain_h/kmax, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7},'DisplayName', 'Measured');
    plot(bpa(i).Ak, bpa(i).strain/kmax, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');
    plot(bpa(i).Ak, bpa(i).strain_p/kmax, '-', 'Color', '#CD34B5', 'LineWidth', 2.5,'DisplayName', 'Predicted');
    % ylabel(tS,'\bf \epsilon^*','Interpreter','tex')
    % xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    % annotation(figS, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
    %     'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

    % Axis config
    set(gca, ...
        'FontSize', 12, ...
        'FontWeight', 'bold', ...
        'FontName', 'Arial', ...
        'LineWidth', 2, ...
        'XMinorTick', 'on', ...
        'YMinorTick', 'on', ...
        'TickLength', [0.025 0.05], ...
        'GridLineStyle','none');
end

%shared axes labels
ylabel(tS,'\bf \epsilon^*','Interpreter','tex')
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex')

% Legend in top-right tile only
lg = legend(tS.Children(1));
lg.Location = 'northeast';
lg.FontSize = 8;
