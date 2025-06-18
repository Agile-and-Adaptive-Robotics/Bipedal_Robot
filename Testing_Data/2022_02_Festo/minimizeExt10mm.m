%% MinimizeExtPinX3.m
%Minimization scheme

clear; clc; close all

%% Initial Baseline Evaluation (for constraint bounds)
[a0, ~] = minimizeExt(0, Inf, Inf, 0, 1);  % All for baseline
baselineScores = a0;  % RMSE, FVU, Max Residual
fprintf('Baseline: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n', mean(baselineScores(:,1)),mean(baselineScores(:,2)),mean(baselineScores(:,3)));

load minimizeExtPin10_results_20250617.mat filtered_results
pick = 12;
sol_actual = filtered_results(pick, 2:5);
sol_actual1 = sol_actual;
[f, bpa] = minimizeExt(sol_actual1(1), sol_actual1(2), sol_actual1(3), sol_actual1(4)*.5831, 1);   % Use solution from Flexor bracket, and compare results
% clear sol_actual
baselineScores1 = f;  % RMSE, FVU, Max Residual
fprintf('Baseline using previous opt: RMSE %.4f, FVU %.4f, Max. Residual %.4f\n', mean(baselineScores1(:,1)),mean(baselineScores1(:,2)),mean(baselineScores1(:,3)));


%% Plot results

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

labels = ["52cm", "52cm", "52cm", "52cm"];
% tileIdxs = [1, 5, 8, 12];  % A, B, C, D
tileIdxs = [1, 5, 10];  % A, B, C
tileSpans = [1 1];      % Span: [rows cols]
tileOrder = [1, 2, 3, 4];
tileLabels = {'(A)', '(B)', '(C)', '(D)'};
% Annotation positions [x, y] in normalized figure units
% xAnn = [0.035, 0.51, 0.035, 0.51];  % (A), (B), (C), (D)
% yAnn = [0.89, 0.89, 0.41, 0.41];    % (A), (B), (C), (D)
xAnn = [0.035, 0.51, 0.265];  % (A), (B), (C)
yAnn = [0.89, 0.89, 0.41];    % (A), (B), (C)
sz = 60;

%for plotting 
%% --- Muscle Length Figure with tiles---
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');

for j = 1
    i = 1;
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    % Predicted
    Lm_p = bpa(i).Lmt_p - 2 * bpa(i).fitn - bpa(i).ten;
    Lm   = bpa(i).Lmt   - 2 * bpa(i).fitn - bpa(i).ten;

    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured'); % Hybrid (gold)
    plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');      % Original
    plot(bpa(i).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');            % Predicted
    ylabel(tL,'\bf Length, m','Interpreter','tex')
    xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex')

    % Tile-specific title and annotation label
    title('\bf 52 cm', 'Interpreter','tex');
    annotation(figL, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

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
    lg = legend;
    lg.Location = 'best';
    lg.FontSize = 8;

end

%% --- Moment Arm Figure with tiles ---
figMA = figure('Name','Moment Arm - 2x2','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(1,1,'TileSpacing','loose','Padding','loose');

for j = 1
    i = 1;
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));

    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured');  % Hybrid
    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');      % Original
    plot(bpa(i).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');                   % Predicted

    title('\bf 52 cm', 'Interpreter', 'tex');

    % Tile-specific title and annotation label
    title(['\bf ' labels(i)], 'Interpreter','tex');
    annotation(figMA, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

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
    
        lg = legend;
        lg.Location = 'best';
        lg.FontSize = 8;
end

%shared axes labels
ylabel(tMA,'\bf Moment arm, m','Interpreter','tex')
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex')

%% --- Torque and relative strain Figure with tiles ---
figT = figure('Name','Torque and relative strain','Color','w');
figT.Position = [100 100 950 400];
tT = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

for j = 1
    i = 1;
    ax = nexttile(tileIdxs(j));
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
    ylabel('\bf Torque, N\cdotm','Interpreter','tex')
    xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    % Tile-specific title and annotation label
    title('Torque', 'Interpreter','tex');
    annotation(figT, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

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
    
        lg = legend;
        lg.Location = 'best';
        lg.FontSize = 8;
end

for j = 2
    i = 1;
    ax = nexttile(j);
    hold on
    
    strain_h = (bpa(i).rest - bpa(i).Lm_h)/bpa(i).rest;
    kmax = (bpa(i).rest - bpa(i).Kmax)/bpa(i).rest;
    scatter(bpa(i).A_h, strain_h/kmax, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7},'DisplayName', 'Measured');
    plot(bpa(i).Ak, bpa(i).strain/kmax, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');
    plot(bpa(i).Ak, bpa(i).strain_p/kmax, '-', 'Color', '#CD34B5', 'LineWidth', 2.5,'DisplayName', 'Predicted');
    
    ylabel('\bf \epsilon^*','Interpreter','tex')
    xlabel('\bf \theta_{k} , \circ','Interpreter','tex')
    title('\bf Relative Strain', 'Interpreter', 'tex');

   
    annotation(figT, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment','center');

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
    
        lg = legend;
        lg.Location = 'best';
        lg.FontSize = 8;          
end



%% --- Strain Figure with tiles ---
% figS = figure('Name','Relative Strain','Color','w');
% figS.Position = [100 100 950 700];
% tS = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');



%% Helper functions
function ff = min1(x, trainIdx)
    if numel(x) == 4 && size(x,1) == 1
        % OK
    else
        error('min1: Input x must be a 1x4 vector');
    end
    Xi0 = x(1) / 100;
    Xi1 = 10^x(2);
    Xi2 = 10^x(3);
    Xi3 = x(4);
    try
        [f_all, ~] = minimizeExt(Xi0, Xi1, Xi2, Xi3, trainIdx); % Nx3 matrix (e.g., 3x3 if 3 training BPAs)
        ff = mean(f_all(trainIdx, :), 1, 'omitnan');              % Return 1x3: [mean RMSE, mean FVU, mean MaxResidual]
        if ~isnumeric(ff) || numel(ff) ~= 3
            ff = [Inf, Inf, Inf];  % Defensive return if shape is wrong
        end
    catch
        ff = [Inf, Inf, Inf];      % Defensive return if minimizeExt throws
    end
end



