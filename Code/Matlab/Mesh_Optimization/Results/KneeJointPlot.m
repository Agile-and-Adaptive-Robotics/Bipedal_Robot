% Input data
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
%Distance from femur origin
knee_x_Pam = [23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]'/1000;
knee_y_Pam = [-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]'/1000;
%Distance from theta1
t1_ICR_x = [29.66	28.54	27.86	27.4	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]'/1000;
t1_ICR_y = [25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	21.32	20.99	21.2	22.33]'/1000;

% Create cubic spline fits for robot data, femur origin to ICR
fcn_x = fit(knee_angle, knee_x_Pam, 'cubicspline');
fcn_y = fit(knee_angle, knee_y_Pam, 'cubicspline');

% Create cubic spline fits for robot data, femur origin to ICR
fcn_t1x = fit(knee_angle, t1_ICR_x, 'cubicspline');
fcn_t1y = fit(knee_angle, t1_ICR_y, 'cubicspline');

% Create cubic spline fits for "Human knee"
% knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
% knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
% fcn1 = fit(knee_angle_x, knee_x, 'cubicspline');
% knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
% knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
% fcn2 = fit(knee_angle_y, knee_y, 'cubicspline');

% Generate dense angle data for parametric curves
theta_dense = linspace(min(knee_angle), max(knee_angle), 200);
x_dense = fcn_x(theta_dense);
y_dense = fcn_y(theta_dense);

% "Human knee" parametric curve (use common theta range for simplicity)
x_t1 = fcn_t1x(theta_dense);
y_t1 = fcn_t1y(theta_dense);

% % Compute Human knee data points for same angles
% knee_x_Bio = fcn1(knee_angle);
% knee_y_Bio = fcn2(knee_angle);

% Plot parametric curve in magenta color
figICR = figure('Name','Knee ICR','Color','w');
% Plot distance from t1 parametric line in light orange with display name
yyaxis left
plot(x_t1, y_t1, '--','Color', '#FA8775', 'LineWidth', 2, 'DisplayName', '\theta_{1}');
hold on;

% Plot white circles at t1 to ICR data points
q_points = plot(t1_ICR_x, t1_ICR_y, 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
q_points.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Add 10% margin around data to prevent touching the edges or tick marks
y_min = min([y_t1]);
y_max = max([y_t1]);

y_margin = 0.30 * (y_max - y_min);

ylim([y_min - y_margin, y_max + y_margin]);

% Plot femur to ICR values
yyaxis right
plot(x_dense, y_dense, 'Color', '#CD34B5', 'LineWidth', 2, 'DisplayName','Femur origin');

% Plot white circles at data points (Robot knee)
h_points = plot(knee_x_Pam, knee_y_Pam, 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
h_points.Annotation.LegendInformation.IconDisplayStyle = 'off';

% Round knee angles to nearest whole number in degrees
knee_angle_deg_rounded = round(rad2deg(knee_angle));

% Create custom datatip row for degree angles
customRow = dataTipTextRow(' ', knee_angle_deg_rounded, 'degrees');

% Replace default datatip rows
h_points.DataTipTemplate.DataTipRows = customRow;
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
x_min = min([x_t1; x_dense]);
x_max = max([x_t1; x_dense]);
y_min = min([y_dense]);
y_max = max([y_dense]);

x_margin = 0.20 * (x_max - x_min);
y_margin = 0.30 * (y_max - y_min);

xlim([x_min - x_margin, x_max + x_margin]);
ylim([y_min - y_margin, y_max + y_margin]);

lgd = legend('Location','best');
title(lgd,'Distance from')