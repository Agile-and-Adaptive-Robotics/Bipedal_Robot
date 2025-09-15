clear;
clc;
close all;

% Input data
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
%Distance from femur origin
% knee_x_Pam = [23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]'/1000;
% knee_y_Pam = [-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]'/1000;
%Distance from theta1
t1_ICR_x = [29.66	28.54	27.86	27.4	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]'/1000;
t1_ICR_y = [25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]'/1000;


% Create cubic spline fits for robot data, theta1 to ICR
fcn_t1x = fit(knee_angle, t1_ICR_x, 'cubicspline');
fcn_t1y = fit(knee_angle, t1_ICR_y, 'cubicspline');

% Generate dense angle data for parametric curves
theta_dense = linspace(min(knee_angle), max(knee_angle), 200);
% x_dense = fcn_x(theta_dense);
% y_dense = fcn_y(theta_dense);

% "Human knee" parametric curve (use common theta range for simplicity)
x_t1 = fcn_t1x(theta_dense);
y_t1 = fcn_t1y(theta_dense);

% Plot parametric curve in magenta color
% figICR = figure('Name','Knee ICR','Color','w');
t = tiledlayout(1, 2);
t.Parent.Color = 'w';
t.Parent.Name = 'Knee ICR';
t.Parent.Position = [78 387 1250 552];
tileLabels = {'(A)', '(B)'};
% Annotation positions [x, y] in normalized figure units
xAnn = [0, 0.54];
yAnn = [0.86, 0.86];
    for j = 1:2
    annotation(gcf, 'textbox', [xAnn(j), yAnn(j), 0.05, 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 10, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    end

% Plot distance from t1 parametric line in light orange with display name
nexttile
plot(x_t1, y_t1, '--','Color', '#FA8775', 'LineWidth', 2, 'DisplayName', '\theta_{1}');
hold on;

% Plot white circles at t1 to ICR data points
q_points = plot(t1_ICR_x, t1_ICR_y, 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
q_points.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Round knee angles to nearest whole number in degrees
knee_angle_deg_rounded = round(rad2deg(knee_angle));

% Create custom datatip row for degree angles
customRow = dataTipTextRow(' ', knee_angle_deg_rounded, 'degrees');

% Replace default datatip rows
% h_points.DataTipTemplate.DataTipRows = customRow;
q_points.DataTipTemplate.DataTipRows = customRow;

% Adjust axes properties
ax = gca;
ax.FontSize = 10;
ax.FontWeight = 'bold';
ax.FontName = 'Arial';
ax.XLabel.String = 'X (m)';
ax.YLabel.String = 'Y (m)';
ax.XLabel.FontSize = 10;
ax.YLabel.FontSize = 10;
ax.XLabel.FontWeight = 'bold';
ax.YLabel.FontWeight = 'bold';
ax.XLabel.FontName = 'Arial';
ax.YLabel.FontName = 'Arial';
ax.TickLength = [0.025 0.05];
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.GridLineStyle = 'none';
ax.LineWidth = 2;

% Remove title
title('');

% Add 10% margin around data to prevent touching the edges or tick marks

y_min = min([y_t1]);
y_max = max([y_t1]);
y_margin = 0.30 * (y_max - y_min);

ylim([y_min - y_margin, y_max + y_margin]);

x_min = min([x_t1]);
x_max = max([x_t1]);
x_margin = 0.30 * (x_max - x_min);

xlim([x_min - x_margin, x_max + x_margin]);

lgd = legend('Location','best');
title(lgd,'Distance from')

nexttile
B0 = imread('KneeLinkage.PNG');
B = imresize(B0,2);
image(B)
axis image
axis off

boxLabels = ["\theta_{3}", "\theta_{2}", "\theta_{1}", "\theta_{4}"];
xbox = [0.6488,0.7212; 0.8176,0.7576; 0.6176,0.6908; 0.8192,0.76720];
ybox = [0.6431,0.6417; 0.5797,0.57645; 0.4438,0.4417; 0.4438,0.4453];
    for j = 1:4
    annotation(gcf, 'textarrow', xbox(j,:), ybox(j,:), 'String', ['\bf ' boxLabels(j)], ...
        'FontSize', 10, 'FontName', 'Arial');
    end
