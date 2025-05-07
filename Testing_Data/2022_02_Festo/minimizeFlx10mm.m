
%Minimization scheme

clear; clc; close all

[a, b, d, ~] = minimizeFlx(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness

%% Use solution from optimizer and check validity on biomimetic knee

load minimizeFlxPin10_results_20250507.mat results_sort_actual

pick = 1;
g = results_sort_actual(pick,2:4);
[u, v, w, bpa] = minimizeFlx(g(1),g(2),g(3));           % Now pull bpa structures out

%% Plot setup;
% Define a colorblind-friendly palette
labels = ["original, ","experiment, ","predicted, "];
labels2 = ["10 mm", "620 kPa", "325 kPa"];
labelz = labels+labels2';

% Color palette (accessible)
c1 = '#FFD700'; % gold
c2 = '#FFB14E'; % orange
c3 = '#FA8775'; % light orange
c4 = '#EA5F94'; % pink
c5 = '#CD34B5'; % magenta
c6 = '#9D02D7'; % violet
c7 = '#0000FF'; % indigo (for Mexp)
c8 = '#000000'; % black (for OpenSim)
c = {c1; c2; c3; c4; c5; c6; c7; c8};

sz = 60; % marker size

%% Plot torque curves, Optimized and validation 
% Plot for 10 mm
i = 1;
figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
plot(bpa(i).Ak, bpa(i).M, '--', 'Color', c2, 'LineWidth', 1.25, 'DisplayName', 'Original');
scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled', 'MarkerFaceColor', c7, 'DisplayName', 'Experiment');
plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2, 'DisplayName', 'Predicted');
title('Torque - 10 mm');
xlabel('\theta_k (°)'); ylabel('Torque, N /cdot m');
legend('Location', 'best'); ylim padded;

% Plot for 20 mm
Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)';           %Angle values directly from O
Bifemsh_T = Tab(:,4)';              %Torque values directly from OpenSim

figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
plot(knee_angle_rT, Bifemsh_T, ':', 'Color', c8, 'LineWidth', 2, 'DisplayName', 'Human Model');
for i = 2:3
    plot(bpa(i).Ak, bpa(i).M, '--', 'Color', c{i}, 'LineWidth', 1.25, 'DisplayName', labelz(i,1));
    scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled', 'MarkerFaceColor', c{9-i}, 'DisplayName',labelz(i,2));
    plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c{7-i}, 'LineWidth', 2, 'DisplayName', labelz(i,3));
end
title('Torque Prediction - 20 mm');
xlabel('\theta_k, °'); ylabel('Torque, N /cdot m)');
legend('Location', 'best'); ylim padded;

%% Plot muscle length, optimization and validation
% 10 mm
figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
i = 1;
Lm = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;
plot(bpa(i).Ak, Lm_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).A_h, bpa(i).Lm_h, sz, hex2rgb(c{7}), 'filled', 'DisplayName', 'Measured');
plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
title('Muscle Length - 10 mm');
xlabel('\theta_k, °'); ylabel('Length, m');
legend('Location', 'best'); ylim padded;

% 20 mm
figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
for i = 2:3
    Lm = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
    Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;
    plot(bpa(i).Ak, Lm_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', labels(i));
    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, hex2rgb(c{7}), 'filled', 'DisplayName',labels(i));
end
title('Muscle Length - 20 mm, Multiple Pressures');
xlabel('\theta_k (°)'); ylabel('Length (m)');
legend('Location', 'best'); ylim padded;


%% Plot moment arm, optimization and validation
% 10 mm
figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
i = 1;
G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));
plot(bpa(i).Ak, G_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).A_h, bpa(i).mA_h, sz, hex2rgb(c{7}), 'filled', 'DisplayName', 'Measured');
plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
title('Moment Arm - 10 mm');
xlabel('\theta_k (°)'); ylabel('Moment Arm (m)');
legend('Location', 'best'); ylim padded;

% ✅ Moment Arm: 20 mm
TabMA = readmatrix('OpenSim_Bifem_MomentArm.txt');
knee_angle_rMA = TabMA(:,2)';
Bifemsh_MA = -TabMA(:,3)';

figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';
plot(knee_angle_rMA, Bifemsh_MA, '--', 'Color', c{8}, 'LineWidth', 2, 'DisplayName', 'Human');
for i = 2:3
    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));
    plot(bpa(i).Ak, G_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', labels(i));
    scatter(bpa(i).A_h, bpa(i).mA_h, sz, hex2rgb(c{7}), 'filled', 'DisplayName',labels(i));
end
title('Moment Arm - 20 mm, Multiple Pressures');
xlabel('\theta_k (°)'); ylabel('Moment Arm (m)');
legend('Location', 'best'); ylim padded;

%% Plot normalized strain, optimization and validation
% 10 mm
figure; hold on;
ax = gca; 
ax.FontWeight = 'bold'; 
ax.FontSize = 11; 
ax.LineWidth = 2; 
ax.Box = 'on'; 
ax.XMinorTick = 'on'; 
ax.YMinorTick = 'on';

i = 1;
kmax = (bpa(i).rest - bpa(i).Kmax) / bpa(i).rest;
Lm   = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;

E    = (bpa(i).rest - Lm) ./ bpa(i).rest / kmax;
E_p  = (bpa(i).rest - Lm_p) ./ bpa(i).rest / kmax;
E_h  = (bpa(i).rest - bpa(i).Lm_h) ./ bpa(i).rest / kmax;

plot(bpa(i).Ak, E_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2.5, 'DisplayName', 'Predicted');
scatter(bpa(i).A_h, E_h, sz, 'filled', 'MarkerFaceColor', hex2rgb(c{5}), 'DisplayName', 'Measured');
plot(bpa(i).Ak, E, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
title('Normalized Strain - 10 mm');
xlabel('\theta_k, °'); ylabel('Normalized Strain');
legend('Location', 'best'); ylim padded;

% 20 mm
figure; hold on;
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 11;
ax.LineWidth = 2;
ax.Box = 'on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';

for i = 2:3
    kmax = (bpa(i).rest - bpa(i).Kmax) / bpa(i).rest;
    Lm   = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
    Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;

    E    = (bpa(i).rest - Lm) ./ bpa(i).rest / kmax;
    E_p  = (bpa(i).rest - Lm_p) ./ bpa(i).rest / kmax;
    E_h  = (bpa(i).rest - bpa(i).Lm_h) ./ bpa(i).rest / kmax;

    plot(bpa(i).Ak, E_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2.5, 'DisplayName', labels(i));
    scatter(bpa(i).A_h, E_h, sz, 'filled', 'MarkerFaceColor', hex2rgb(c{5}), 'DisplayName', labels(i));
end
title('Normalized Strain - 20 mm, Multiple Pressures');
xlabel('\theta_k, °'); ylabel('Normalized Strain');
legend('Location', 'best'); ylim padded;



function rgb = hex2rgb(hex)
% Converts hex color string to an RGB vector
    hex = char(hex);
    if hex(1) == '#'
        hex = hex(2:end);
    end
    if numel(hex) ~= 6
        error('Hex color must be 6 characters long');
    end
    r = hex2dec(hex(1:2)) / 255;
    g = hex2dec(hex(3:4)) / 255;
    b = hex2dec(hex(5:6)) / 255;
    rgb = [r, g, b];
end
