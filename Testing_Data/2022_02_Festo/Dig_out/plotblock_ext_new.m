%% Plot results (2026-09-16: two figures per metric -- (1) the BPAs trained &
% validated (allBPA) with per-tile subtitles naming each tile's training and
% validation folds, and (2) the tests left out of allBPA)
c = cell(8,1);
c{1} = '#FFD700'; % gold → Hybrid
c{2} = '#FFB14E'; % orange
c{3} = '#FA8775'; % light orange
c{4} = '#EA5F94'; % pink
c{5} = '#CD34B5'; % magenta → Predicted
c{6} = '#9D02D7'; % magenta 2
c{7} = '#0000FF'; % indigo → Measured
c{8} = '#000000'; % black
sz = 60;

% --- figure groups: {allBPA trained/validated, tests left out of allBPA} ---
grpTests = {allBPA, setdiff(1:numel(labels), allBPA)};
grpNames = {'Trained \& Validated', 'Left Out of CV'};

% --- fold structure for the training/validation subtitles ---
foldList = nchoosek(allBPA, numHold);
nFolds = size(foldList, 1);
trainFolds = cell(numel(labels),1); valFolds = cell(numel(labels),1);
for ii = 1:numel(labels)
    if any(allBPA == ii)
        valFolds{ii} = find(any(foldList == ii, 2));
        trainFolds{ii} = setdiff(1:nFolds, valFolds{ii});
    end
end
tileIdxsFor = @(n) arrayfun(@(j) (ceil(j/2)-1)*7 + 1 + 4*mod(j-1,2), 1:n);

for g = 1:2
idxs = grpTests{g};
el = numel(idxs);
tileIdxs = tileIdxsFor(el);
tileSpans = [1 3];      % Span: [rows cols]

%% --- Torque ---
figT = figure('Name', ['Torque - ' grpNames{g}], 'Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = idxs(j);
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

    % Tile-specific title and train/validation subtitle
    title(['\bf ' labels(i)], 'Interpreter','tex');
    if g == 1
        subT = subtitle(['Training: ' strjoin(string(trainFolds{i}), ' ') ...
            '  |  Validation: ' strjoin(string(valFolds{i}), ' ')]);
    else
        subT = subtitle('Left out of allBPA');
    end
    subT.FontSize = 9;

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

%% --- Muscle Length ---
figL = figure('Name', ['Muscle Length - ' grpNames{g}], 'Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = idxs(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    % Predicted
    Lm_p = bpa(i).Lmt_p - 2 * bpa(i).fitn - bpa(i).ten;
    Lm   = bpa(i).Lmt   - 2 * bpa(i).fitn - bpa(i).ten;

    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured');
    plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');
    plot(bpa(i).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');

    % Tile-specific title and train/validation subtitle
    title(['\bf ' labels(i)], 'Interpreter','tex');
    if g == 1
        subT = subtitle(['Training: ' strjoin(string(trainFolds{i}), ' ') ...
            '  |  Validation: ' strjoin(string(valFolds{i}), ' ')]);
    else
        subT = subtitle('Left out of allBPA');
    end
    subT.FontSize = 9;

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

%% --- Moment Arm ---
figMA = figure('Name', ['Moment Arm - ' grpNames{g}], 'Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = idxs(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));

    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c{7},'DisplayName', 'Measured');  % Hybrid
    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');
    plot(bpa(i).Ak, G_p, '-', 'Color', c{5}, 'LineWidth', 2.5,'DisplayName', 'Predicted');

    % Tile-specific title and train/validation subtitle
    title(['\bf ' labels(i)], 'Interpreter','tex');
    if g == 1
        subT = subtitle(['Training: ' strjoin(string(trainFolds{i}), ' ') ...
            '  |  Validation: ' strjoin(string(valFolds{i}), ' ')]);
    else
        subT = subtitle('Left out of allBPA');
    end
    subT.FontSize = 9;

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

%% --- Strain ---
figS = figure('Name', ['Strain - ' grpNames{g}], 'Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(el/2),7,'TileSpacing','tight','Padding','tight');

for j = 1:el
    i = idxs(j);
    ax = nexttile(tileIdxs(j), tileSpans);
    hold on

    strain_h = (bpa(i).rest - bpa(i).Lm_h)/bpa(i).rest;
    kmax = (bpa(i).rest - bpa(i).Kmax)/bpa(i).rest;
    scatter(bpa(i).A_h, strain_h/kmax, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', c{7},'DisplayName', 'Measured');
    plot(bpa(i).Ak, bpa(i).strain/kmax, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2,'DisplayName', 'Original');
    plot(bpa(i).Ak, bpa(i).strain_p/kmax, '-', 'Color', '#CD34B5', 'LineWidth', 2.5,'DisplayName', 'Predicted');

    % Tile-specific title and train/validation subtitle
    title(['\bf ' labels(i)], 'Interpreter','tex');
    if g == 1
        subT = subtitle(['Training: ' strjoin(string(trainFolds{i}), ' ') ...
            '  |  Validation: ' strjoin(string(valFolds{i}), ' ')]);
    else
        subT = subtitle('Left out of allBPA');
    end
    subT.FontSize = 9;

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

end
