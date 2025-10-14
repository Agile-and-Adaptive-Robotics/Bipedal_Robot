
%Minimization scheme

clear; clc; close all

[a, b, d, ~] = minimizeFlx(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness

%% Use solution from optimizer and check validity on biomimetic knee

load minimizeFlxPin10_results_20251014_1transform.mat results_sort_actual

pick = 13;
g = results_sort_actual(pick,2:4);
[u, v, w, bpa] = minimizeFlx(0,g(2),g(3));           % Now pull bpa structures out

%% Plot setup;
% Define a colorblind-friendly palette
labels = ["Original, ","Measured, ","Predicted, "];
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
%% === 1x2 Tile Plot for Torque Curves (10 mm and 20 mm) ===
Tab = readmatrix('OpenSim_Bifem_Results.txt');
knee_angle_rT = Tab(:,2)';    % Angle values from OpenSim
Bifemsh_T = Tab(:,4)';        % Torque values from OpenSim

figT = figure('Name', 'Torque Comparison', 'Color', 'w');
figT.Position = [100 100 950 450];
tT = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

% Tile 1: 10 mm
ax1 = nexttile(1);
hold on

i = 1;
plot(bpa(i).Ak, bpa(i).M, '--', 'Color', c2, 'LineWidth', 2, 'DisplayName', 'Original');
plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled', 'MarkerFaceColor', c7, 'DisplayName', 'Measured');

title('\phi 10 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Torque, N \cdot m', 'Interpreter', 'tex', ...
       'FontWeight', 'bold', 'FontSize', 11, 'FontName', 'Arial');
xlabel('\bf \theta_{k} , \circ', 'Interpreter', 'tex', ...
       'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);

legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

% Tile 2: 20 mm
ax2 = nexttile(2);
hold on

% plot(knee_angle_rT, Bifemsh_T, ':', 'Color', c8, 'LineWidth', 4, 'DisplayName', 'Human Model');
% 
% for i = 2:3
%     plot(bpa(i).Ak, bpa(i).M, '--', 'Color', c{i}, 'LineWidth', 2, 'DisplayName', labelz(i,1));
%     scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled', 'MarkerFaceColor', c{9-i}, 'DisplayName', labelz(i,2));
%     plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c{7-i}, 'LineWidth', 2, 'DisplayName', labelz(i,3));
% end

% 1) human model
plot(knee_angle_rT, Bifemsh_T, ':', 'Color', c8, 'LineWidth', 4, 'DisplayName', 'Human Model');

% 2) only the 20 mm @ 620 kPa curves (bpa(2))
i = 2;
plot(bpa(i).Ak, bpa(i).M, '--', 'Color', c{i}, 'LineWidth', 2, 'DisplayName', 'Original, 620 kPa');
plot(bpa(i).Ak, bpa(i).M_p(:,3), '-', 'Color', c{7-i}, 'LineWidth', 2, 'DisplayName', 'Predicted, 620 kPa');
scatter(bpa(i).Aexp, bpa(i).Mexp, sz, 'filled','MarkerFaceColor', c7, 'DisplayName', 'Measured, 620 kPa');

% 3) overlay exactly four measured points < 620 kPa as one scatter
%  - 560 kPa, 421 kPa, nearest ~325 kPa with angle < 281 kPa, and 281 kPa
angle_data    = readmatrix('Results_table_FullSize.xlsx', ...
    'Sheet','FlxTest20mm_42cm (3)', 'Range','C7:X7');
pressure_data = readmatrix('Results_table_FullSize.xlsx', ...
    'Sheet','FlxTest20mm_42cm (3)', 'Range','C8:X8');
torque_data   = readmatrix('Results_table_FullSize.xlsx', ...
    'Sheet','FlxTest20mm_42cm (3)', 'Range','C17:X17');

% find index for 281 kPa
angle281 = angle_data(pressure_data == 281);
% mask for 560, 421, 281 kPa
mask = ismember(pressure_data, [560, 421, 281]);
% find one ~325 kPa point with angle < angle281
near325 = abs(pressure_data - 325) <= 5;
idx325 = find(near325 & angle_data < angle281, 1, 'first');
if ~isempty(idx325)
    mask(idx325) = true;
end
% plot all four as one scatter
scatter(angle_data(mask), torque_data(mask), sz, 'filled', ...
    'MarkerFaceColor', c4, 'DisplayName', 'Measured < 620 kPa');

title('\phi 20 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Torque, N \cdot m', 'Interpreter', 'tex', ...
       'FontWeight', 'bold', 'FontSize', 11, 'FontName', 'Arial');
xlabel('\bf \theta_{k} , \circ', 'Interpreter', 'tex', ...
       'FontSize', 12, 'FontName', 'Arial');

set(gca, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on','TickLength', [0.025 0.05]);

legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

% Shared x-label

% Textboxes (A) and (B)
annotation(figT, 'textbox', [0.01, 0.9, 0.05, 0.05], 'String', '\bf (A)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figT, 'textbox', [0.46, 0.9, 0.05, 0.05], 'String', '\bf (B)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');


%% Plot muscle length, optimization and validation
figL = figure('Name', 'Muscle Length', 'Color', 'w');
figL.Position = [100 100 950 450];
tL = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

% Tile 1: 10 mm
ax1 = nexttile(1);
hold on
i = 1;
Lm = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;

plot(bpa(i).Ak, Lm, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2, 'DisplayName', 'Original');
scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceColor', c{7}, 'DisplayName', 'Measured');
plot(bpa(i).Ak, Lm_p, '-', 'Color', c{5}, 'LineWidth', 2, 'DisplayName', 'Predicted');

title('\phi 10 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Length, m', 'FontWeight', 'bold', 'FontSize', 11, 'FontName', 'Arial');

set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

% Tile 2: 20 mm
ax2 = nexttile(2);
hold on
for i = 2:3
    Lm = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
    Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;

    plot(bpa(i).Ak, Lm, '--', 'Color', c{i}, 'LineWidth', 2, 'DisplayName', labelz(i,1));
    scatter(bpa(i).A_h, bpa(i).Lm_h, sz, 'filled', 'MarkerFaceColor', c{9-i}, 'DisplayName', labelz(i,2));
    plot(bpa(i).Ak, Lm_p, '-', 'Color', c{7-i}, 'LineWidth', 2, 'DisplayName', labelz(i,3));
end

title('\phi 20 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontWeight', 'bold', 'FontSize', 12, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

xlabel(tL, '\bf \theta_{k} , \circ', 'Interpreter', 'tex', ...
       'FontSize', 12, 'FontName', 'Arial');

annotation(figL, 'textbox', [0.01, 0.9, 0.05, 0.05], 'String', '\bf (A)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figL, 'textbox', [0.46, 0.9, 0.05, 0.05], 'String', '\bf (B)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');


%% Plot moment arm, optimization and validation
figMA = figure('Name', 'Moment Arm', 'Color', 'w');
figMA.Position = [100 100 950 450];
tMA = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

TabMA = readmatrix('OpenSim_Bifem_MomentArm.txt');
knee_angle_rMA = TabMA(:,2)';
Bifemsh_MA = -TabMA(:,3)';

% Tile 1: 10 mm
ax1 = nexttile(1);
hold on
i = 1;
G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));

plot(bpa(i).Ak, G_p, '-', 'Color', c2, 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceColor', c7, 'DisplayName', 'Measured');
plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', c5, 'LineWidth', 2, 'DisplayName', 'Original');

title('\phi 10 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Moment Arm (m)', 'FontWeight', 'bold', 'FontSize', 11, 'FontName', 'Arial');

set(gca, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

% Tile 2: 20 mm
ax2 = nexttile(2);
hold on
plot(knee_angle_rMA, Bifemsh_MA, '--', 'Color', c{8}, 'LineWidth', 2, 'DisplayName', 'Human');

for i = 2:3
    G_p = hypot(bpa(i).mA_p(:,1), bpa(i).mA_p(:,2));
    plot(bpa(i).Ak, bpa(i).mA, '--', 'Color', c{i}, 'LineWidth', 2, 'DisplayName', labelz(i,1));
    scatter(bpa(i).A_h, bpa(i).mA_h, sz, 'filled', 'MarkerFaceColor', c{9-i}, 'DisplayName', labelz(i,2));
    plot(bpa(i).Ak, G_p, '-', 'Color', c{7-i}, 'LineWidth', 2, 'DisplayName', labelz(i,3));
end

title('\phi 20 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

xlabel(tMA, '\bf \theta_{k} , \circ', 'Interpreter', 'tex', ...
       'FontSize', 12, 'FontName', 'Arial');

annotation(figMA, 'textbox', [0.01, 0.9, 0.05, 0.05], 'String', '\bf (A)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figMA, 'textbox', [0.46, 0.9, 0.05, 0.05], 'String', '\bf (B)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Plot normalized strain, optimization and validation
figS = figure('Name', 'Normalized Strain', 'Color', 'w');
figS.Position = [100 100 950 450];
tS = tiledlayout(1,2,'TileSpacing','loose','Padding','loose');

% Tile 1: 10 mm
ax1 = nexttile(1);
hold on
i = 1;
kmax = (bpa(i).rest - bpa(i).Kmax) / bpa(i).rest;
Lm   = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;

E    = (bpa(i).rest - Lm) ./ bpa(i).rest / kmax;
E_p  = (bpa(i).rest - Lm_p) ./ bpa(i).rest / kmax;
E_h  = (bpa(i).rest - bpa(i).Lm_h) ./ bpa(i).rest / kmax;

plot(bpa(i).Ak, E_p, '-', 'Color', c2, 'LineWidth', 2, 'DisplayName', 'Predicted');
scatter(bpa(i).A_h, E_h, sz, 'filled', 'MarkerFaceColor', c7, 'DisplayName', 'Measured');
plot(bpa(i).Ak, E, '--', 'Color', c5, 'LineWidth', 2, 'DisplayName', 'Original');

title('\phi 10 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
ylabel('Normalized Strain', 'FontWeight', 'bold', 'FontSize', 11, 'FontName', 'Arial');

set(gca, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

% Tile 2: 20 mm
ax2 = nexttile(2);
hold on
for i = 2:3
    kmax = (bpa(i).rest - bpa(i).Kmax) / bpa(i).rest;
    Lm   = bpa(i).Lmt - 2*bpa(i).fitn - bpa(i).ten;
    Lm_p = bpa(i).Lmt_p - 2*bpa(i).fitn - bpa(i).ten;

    E    = (bpa(i).rest - Lm) ./ bpa(i).rest / kmax;
    E_p  = (bpa(i).rest - Lm_p) ./ bpa(i).rest / kmax;
    E_h  = (bpa(i).rest - bpa(i).Lm_h) ./ bpa(i).rest / kmax;

    plot(bpa(i).Ak, E, '--', 'Color', c{i}, 'LineWidth', 2, 'DisplayName', labelz(i,1));
    scatter(bpa(i).A_h, E_h, sz, 'filled', 'MarkerFaceColor', c{9-i}, 'DisplayName', labelz(i,2));
    plot(bpa(i).Ak, E_p, '-', 'Color', c{7-i}, 'LineWidth', 2, 'DisplayName', labelz(i,3));
end

title('\phi 20 mm', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial');
set(gca, 'FontWeight', 'bold', 'FontSize', 11, 'LineWidth', 2, ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
legend('Location', 'best', 'FontSize', 8);
xlim([-120 10]);

xlabel(tS, '\bf \theta_{k} , \circ', 'Interpreter', 'tex', ...
       'FontSize', 12, 'FontName', 'Arial');

annotation(figS, 'textbox', [0.01, 0.9, 0.05, 0.05], 'String', '\bf (A)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
annotation(figS, 'textbox', [0.46, 0.9, 0.05, 0.05], 'String', '\bf (B)', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');


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
