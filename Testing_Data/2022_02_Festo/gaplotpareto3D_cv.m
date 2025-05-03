function stop = gaplotpareto3D_cv(options, state, flag, holdoutIdx, Color_var, labels)
% Custom 3D Pareto plot for cross-validation — shows holdout info and live Pareto front

persistent ax hScatter

stop = false;

% Reset persistent vars on 'init'
if strcmp(flag, 'init')
    hScatter = [];
    ax = [];
    return;
end

% Only plot on 'iter' or 'done'
if ~ismember(flag, {'iter','done'})
    return;
end

scores = state.Score;
if isempty(scores) || size(scores,2) < 3 || all(isnan(scores), 'all')
    return;
end

% Validate labels
if holdoutIdx > numel(labels)
    warning('Invalid holdoutIdx index for label.');
    return;
end

% Create figure and axes if needed
if isempty(ax) || ~isgraphics(ax, 'axes')
    fig = figure(100 + holdoutIdx);
    clf(fig);
    ax = axes('Parent', fig);
    hold(ax, 'on'); grid(ax, 'on');
    xlabel(ax, 'RMSE'); ylabel(ax, 'FVU'); zlabel(ax, 'Max Residual');
    view(ax, [135 25]);
    title(ax, sprintf('Pareto Front – CV Holdout: %s', labels(holdoutIdx)), 'FontWeight','bold');
    ax.FontWeight = 'bold'; ax.FontSize = 11;
end

% Get color from training set mapping
trainSet = Color_var(holdoutIdx, :);
cmap = lines(numel(labels));
thisColor = mean(cmap(trainSet, :), 1);  % mean color of training set

% Clear previous plot
if ~isempty(hScatter) && isgraphics(hScatter)
    delete(hScatter);
end

% Plot new 3D scatter
hScatter = scatter3(ax, scores(:,1), scores(:,2), scores(:,3), ...
    36, repmat(thisColor, size(scores,1), 1), 'filled', ...
    'MarkerEdgeColor', 'k', ...
    'DisplayName', ['Holdout: ' char(labels(holdoutIdx))]);

legend(ax, 'show', 'Location', 'northeast');
drawnow;

% Don't delete on final generation
if strcmp(flag, 'done')
    hScatter = [];
end

end
