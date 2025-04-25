
%Minimization scheme

clear; clc; close all

[a, b, d, ~] = minimizeFlx(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness

%% Use solution from optimizer and check validity on biomimetic knee

load minimizeFlxPin10_results.mat sol_actual

g = sol_actual(1,1:3);
[u, v, w, bpa] = minimizeFlx(g(1),g(2),g(3));           % Now pull bpa structures out

%% Plot setup;
% Define a colorblind-friendly palette
labels = ["10 mm", "20 mm, 625 kPa", "20 mm, 325 kPa"];

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
figure; hold on; grid on;
i = 1;
plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).Aexp, bpa(i).Mexp, sz, hex2rgb(c{7}), 'filled', 'DisplayName', 'Experiment');
plot(bpa(i).Ak, bpa(i).M, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
title('Torque Prediction - 10 mm');
xlabel('\theta_k (°)'); ylabel('Torque (Nm)');
legend('Location', 'best'); ylim padded;

% Plot for 20 mm
Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)';           %Angle values directly from O
Bifemsh_T = Tab(:,4)';              %Torque values directly from OpenSim

figure; hold on; grid on;
plot(knee_angle_rT, Bifemsh_T, '--', 'Color', hex2rgb(c{8}), 'LineWidth', 2, 'DisplayName', 'Human');
for i = 2:3
    plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', labels(i));
    scatter(bpa(i).Aexp, bpa(i).Mexp, sz, hex2rgb(c{7}), 'filled', 'DisplayName',labels(i));
end
title('Torque Prediction - 20 mm, Multiple Pressures');
xlabel('\theta_k (°)'); ylabel('Torque (Nm)');
legend('Location', 'best'); ylim padded;

%% Plot muscle length, optimization and validation
% 10 mm
figure; hold on; grid on;
i = 1;
Lm = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;
plot(bpa(i).Ak, Lm_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).A_h, bpa(i).Lm_h, sz, hex2rgb(c{7}), 'filled', 'DisplayName', 'Measured');
plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
title('Muscle Length - 10 mm');
xlabel('\theta_k (°)'); ylabel('Length (m)');
legend('Location', 'best'); ylim padded;

% 20 mm
figure; hold on; grid on;
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
figure; hold on; grid on;
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

figure; hold on; grid on;
plot(knee_angle_rMA, Bifemsh_MA, '--', 'Color', c{8}, 'LineWidth', 2, 'DisplayName', 'Human');
for i = 2:3
    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));
    plot(bpa(i).Ak, G_p, '-', 'Color', hex2rgb(c{i}), 'LineWidth', 2, 'DisplayName', labels(i));
    scatter(bpa(i).A_h, bpa(i).mA_h, sz, hex2rgb(c{7}), 'filled', 'DisplayName',labels(i));
end
title('Moment Arm - 20 mm, Multiple Pressures');
xlabel('\theta_k (°)'); ylabel('Moment Arm (m)');
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
