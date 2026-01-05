%% Visualize Curvature vs. Knee Angle
clear; clc;

% --- Parameters ---
theta_deg = -120:1:60;                 % Knee angle [°] (flexion negative)
angleRad = deg2rad(theta_deg);         % Convert to radians

% --- Call curvature function ---
kappa = Curvature(angleRad);           % Curvature [1/m]
R = 1 ./ kappa;                        % Radius [m]

% --- Plot curvature ---
figure('Color','w'); hold on;
yyaxis left
plot(theta_deg, kappa, 'LineWidth', 2, 'Color', [0.2 0.2 1]);
ylabel('\bf Curvature, \kappa (1/m)','Interpreter','tex')
ylim([0 40])

yyaxis right
plot(theta_deg, R*1000, '--', 'LineWidth', 2, 'Color', [1 0.3 0.3]);
ylabel('\bf Radius of Curvature, R (mm)','Interpreter','tex')
ylim([30 Inf])

xlabel('\bf Knee Angle, \theta_k (°)', 'Interpreter','tex')
title('\bf Estimated BPA Curvature vs. Knee Angle')
grid on;
set(gca, 'FontSize',12, 'FontWeight','bold', 'LineWidth', 1.5, ...
    'XMinorTick','on', 'YMinorTick','on', 'FontName','Arial')
legend({'Curvature (\kappa)', 'Radius (R)'}, 'Location','NorthWest');

%% Curvature in bending BPA estimate
function kappa = Curvature(angleRad)
    % Parameters
    R_min = 0.03;               % [m] min radius of curvature (30 mm)
    kappa_max = 1 / R_min;      % Max curvature
    angleMid = deg2rad(-60);    % Midpoint of transition (steepest change)
    sharpness = 12;             % Larger value => steeper transition

    % Logistic-style transition: 0 → kappa_max
    kappa = kappa_max ./ (1 + exp(sharpness * (angleRad - angleMid)));

    % Optional clamp (shouldn't be necessary, but safe)
    kappa = min(kappa, kappa_max);
end