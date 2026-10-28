function vasmed_rp7_fig_torque_20261028()
% vasmed_rp7_fig_torque_20261028 -- torque-angle validation figure
% (2026-10-28). Panels:
%   (A) Robot vastus-medialis-extension torque vs knee angle on the dense
%       131-pose verification grid: human single-muscle target (vas_med_r,
%       OpenSim Gait2392, pinned source), the direct-drive saved design
%       (1,1, regraded), the best fixed-route sizes-only cell (BEFORE
%       placement optimization), and the attachment-optimized winner
%       (AFTER). Y-label reads EXACTLY 'Torque, N m'.
%   (B) Margin (tau_robot - tau_human) for the same designs; zero line
%       (plain target) and the 5% route-design margin curve.
% Style: Arial >= 10 pt at final size, bold margin panel letters (A)/(B)
% left of each tile straddling the plot-box top (FlxPin_group
% convention), Tol bright CVD-safe palette, distinct markers AND line
% styles (never color alone). Exports: FIG/PDF/EPS/PNG in this folder.

here = fileparts(mfilename('fullpath'));
T = readtable(fullfile(here, 'vasmed_rp6_dense_curves.csv'), ...
    'VariableNamingRule', 'preserve', 'Delimiter', ',', 'TextType', 'char');
dcol = string(T.design);
tags = unique(dcol, 'stable');
deg = numel(tags);
sayT = cell(deg, 1);
for k = 1:deg
    sayT{k} = strrep(tags{k}, '_', '-');
end
% Map display names
dispNames = {};
for k = 1:deg
    t = tags{k};
    tok = regexp(t, 'nB(\d+)_G(\d+)', 'tokens', 'once');
    if isempty(tok), tok = {'1', '1'}; end
    if startsWith(t, 'saved')
        dn = sprintf('direct-drive saved design (n=1, G=1)');
    elseif startsWith(t, 'fixedroute')
        dn = sprintf('fixed route, sizes re-solved (n=%s, G=%s)', tok{1}, tok{2});
    else
        dn = sprintf('placement-optimized (n=%s, G=%s)', tok{1}, tok{2});
    end
    dispNames{end + 1, 1} = dn; %#ok<AGROW>
end

% Tol bright palette, CVD-safe
col = [214, 88, 108; 0, 110, 191; 247, 143, 21; 0, 158, 137; ...
    153, 0, 114; 0, 0, 0] / 255;
lst = {'-', '--', '-.', ':', '-', '--'};
mk = {'o', 's', '^', 'd', 'v', 'x'};

fig = figure('Units', 'centimeters', 'Position', [2, 2, 17.8, 7.6], ...
    'Color', 'w', 'PaperPositionMode', 'auto');

axA = subplot(1, 2, 1); hold(axA, 'on');
axB = subplot(1, 2, 2); hold(axB, 'on');

tauH = [];
for k = 1:deg
    m = dcol == tags{k};
    if ~any(m), continue; end
    ang = T.knee_angle_deg(m); tb = T.tau_robot_b_Nm(m);
    th = T.tau_h_Nm(m);
    if isempty(tauH), tauH = th(1); tauHang = ang; end
    ii = 1:6:numel(ang);
    plot(axA, ang, tb, lst{k}, 'Color', col(k, :), 'LineWidth', 1.6, ...
        'DisplayName', dispNames{k});
    plot(axA, ang(ii), tb(ii), mk{k}, 'MarkerEdgeColor', col(k, :), ...
        'MarkerFaceColor', 'w', 'MarkerSize', 4.5, 'LineStyle', 'none', ...
        'HandleVisibility', 'off');
    plot(axB, ang, tb - th, lst{k}, 'Color', col(k, :), 'LineWidth', 1.6, ...
        'HandleVisibility', 'off');
end
plot(axA, tauHang, tauH, 'k-', 'LineWidth', 2.2, ...
    'DisplayName', 'human vas medialis (OpenSim Gait2392)');
plot(axB, tauHang, 0.05 * tauH, 'k:', 'LineWidth', 1.6);
text(axB, tauHang(end) - 4, 0.05 * tauH(end) + 2, '+5\% design margin', ...
    'FontSize', 10, 'FontName', 'Arial', 'HorizontalAlignment', 'right');
plot(axB, xlim(axB), [0 0], 'k-', 'LineWidth', 0.8);

xlabel(axA, 'Knee angle, degrees (flexion negative)', 'FontSize', 10, ...
    'FontName', 'Arial');
ylabel(axA, 'Torque, N m', 'FontSize', 10, 'FontName', 'Arial');
xlabel(axB, 'Knee angle, degrees (flexion negative)', 'FontSize', 10, ...
    'FontName', 'Arial');
ylabel(axB, 'Torque margin, N m', 'FontSize', 10, 'FontName', 'Arial');
set([axA, axB], 'FontSize', 10, 'FontName', 'Arial', 'Box', 'on', ...
    'TickDir', 'out', 'Layer', 'bottom');
legend(axA, 'Location', 'northoutside', 'Orientation', 'vertical', ...
    'FontSize', 8.5, 'FontName', 'Arial', 'NumColumns', 1);
grid(axA, 'on'); grid(axB, 'on');
axA.Title.String = 'Torque about the biomimetic knee (cable arm, moving ICR)';
axB.Title.String = 'Margin vs human target over the full domain';

% Panel letters (FlxPin_group convention): margin labels LEFT of each
% tile, straddling the plot-box top.
drawnow;
aP = axA.Position; bP = axB.Position;
annotation(fig, 'textbox', [aP(1) - 0.065, aP(2) + aP(4) - 0.028, 0.05, 0.07], ...
    'String', '(A)', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', ...
    'Arial', 'EdgeColor', 'none', 'VerticalAlignment', 'bottom');
annotation(fig, 'textbox', [bP(1) - 0.055, bP(2) + bP(4) - 0.028, 0.05, 0.07], ...
    'String', '(B)', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', ...
    'Arial', 'EdgeColor', 'none', 'VerticalAlignment', 'bottom');

set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
saveas(fig, fullfile(here, 'vasmed_reverse_pulley_torque_20261028.fig'));
exportgraphics(fig, fullfile(here, 'vasmed_reverse_pulley_torque_20261028.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(here, 'vasmed_reverse_pulley_torque_20261028.png'), ...
    'Resolution', 300);
print(fig, fullfile(here, 'vasmed_reverse_pulley_torque_20261028.eps'), ...
    '-depsc2', '-vector');
fprintf('RP7 torque figure DONE\n');
end
