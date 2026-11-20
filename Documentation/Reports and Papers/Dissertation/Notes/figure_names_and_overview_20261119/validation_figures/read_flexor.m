function read_flexor(figure_file, output_json)
% Read native arrays directly from the accepted MATLAB publication figure.
f = openfig(figure_file, 'invisible');
close_guard = onCleanup(@() close(f));
ax = findall(f, 'Type', 'axes');
assert(isscalar(ax), 'One torque axis is required.');
data = struct('xlim', ax.XLim, 'ylim', ax.YLim, 'x_ticks', ax.XTick, ...
    'y_ticks', ax.YTick, 'series', []);
objects = findall(ax,'-property','XData','-property','YData');
for index = 1:numel(objects)
    h = objects(index);
    % Three direct children are one-NaN legend keys. The actual curves and
    % measured trials are children of the copied plotting group.
    if numel(h.XData) < 15 || ~all(isfinite(h.YData))
        continue
    end
    if isa(h,'matlab.graphics.chart.primitive.Scatter')
        if contains(h.DisplayName,'Prediction')
            name = 'prediction';
        else
            name = 'measured';
        end
    elseif strcmp(h.LineStyle,':')
        name = 'human';
    else
        name = 'design';
    end
    record = struct('name',name, 'type',char(h.Type), ...
        'x',h.XData(:)', 'y',h.YData(:)', 'area',[], 'alpha',1);
    if isa(h,'matlab.graphics.chart.primitive.Scatter')
        record.area = h.SizeData(:)';
        record.alpha = h.MarkerFaceAlpha;
    end
    data.series = [data.series record];
end
assert(numel(data.series) == 4, 'Expected two curves and two trial series.');
fid = fopen(output_json,'w','n','UTF-8');
assert(fid >= 0, 'Unable to write the source data.');
file_guard = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n', jsonencode(data));
end
