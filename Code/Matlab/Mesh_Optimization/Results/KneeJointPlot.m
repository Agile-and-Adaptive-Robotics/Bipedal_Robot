% Input data
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam = [0.0065; 0.0083; 0.0094; 0.0101; 0.0120; 0.0140; 0.0161; 0.0220; 0.0269; 0.0302; 0.0311; 0.0295; 0.0253; 0.0189; 0.0109; 0.0021];
knee_y_Pam = [-0.3981; -0.3968; -0.3961; -0.3957; -0.3949; -0.3943; -0.3941; -0.3950; -0.3982; -0.4034; -0.4098; -0.4165; -0.4227; -0.4273; -0.4297; -0.4289];

% Create cubic spline fits for original data
fcn_x = fit(knee_angle, knee_x_Pam, 'cubicspline');
fcn_y = fit(knee_angle, knee_y_Pam, 'cubicspline');

% Create cubic spline fits for "Human knee"
knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x, knee_x, 'cubicspline');

knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y, knee_y, 'cubicspline');

% Generate dense angle data for parametric curves
theta_dense = linspace(min(knee_angle), max(knee_angle), 200);
x_dense = fcn_x(theta_dense);
y_dense = fcn_y(theta_dense);

% "Human knee" parametric curve (use common theta range for simplicity)
x_dense_human = fcn1(theta_dense);
y_dense_human = fcn2(theta_dense);

% Compute Human knee data points for same angles
knee_x_Bio = fcn1(knee_angle);
knee_y_Bio = fcn2(knee_angle);

% Plot parametric curve in magenta color
figICR = figure('Name','Knee ICR','Color','w');
% Plot "Human knee" parametric line in light orange with display name
plot(x_dense_human, y_dense_human, '--','Color', '#FA8775', 'LineWidth', 2, 'DisplayName', 'Human knee');
hold on;

% Plot Robot knee values
plot(x_dense, y_dense, 'Color', '#CD34B5', 'LineWidth', 2, 'DisplayName','RobotKnee');

% Plot white circles at data points (Robot knee)
h_points = plot(knee_x_Pam, knee_y_Pam, 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
h_points.Annotation.LegendInformation.IconDisplayStyle = 'off';
% Plot white circles at "Human knee" data points
q_points = plot(knee_x_Bio, knee_y_Bio, 'wo', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
q_points.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
x_min = min([x_dense; x_dense_human]);
x_max = max([x_dense; x_dense_human]);
y_min = min([y_dense; y_dense_human]);
y_max = max([y_dense; y_dense_human]);

x_margin = 0.10 * (x_max - x_min);
y_margin = 0.10 * (y_max - y_min);

xlim([x_min - x_margin, x_max + x_margin]);
ylim([y_min - y_margin, y_max + y_margin]);

lgd = legend('Location','best');
