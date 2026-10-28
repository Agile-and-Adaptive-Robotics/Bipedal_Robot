function vasmed_rp5_fig_schematic_20261028()
% vasmed_rp5_fig_schematic_20261028 -- vector reverse-pulley schematic
% (2026-10-28). Based on the ACTUAL system: the MonoPam_pulley /
% BiPam_pulley2 hardware contract (Ben's 2026-09-24 block-and-tackle rig:
% n parallel 20 mm BPAs + tackle mounted on the proximal body, one tendon
% line crossing the joint to the tibial insertion) and the vas-med
% distal-ring route of record (femur-side BPA wrap contacts, tibia-side
% cable ring contacts). NOTE: the file named in the task prompt,
% "Documentation\Posters\Bolen et al. - 2024 - Self-balancing Robot Leg
% .pdf", does not exist in this repo (checked 2026-10-28); the closest
% poster on disk is 4008_Bolen_Ben_PDF.pdf ("Torque Values About
% Artificial Knee Using Braided Pneumatic Actuators"). The schematic
% therefore follows the class hardware contract, and says so.
%
% Panel (A): tackle free-body diagram -- n BPAs, G line parts, fixed and
% moving anchors, forces: per-BPA force F_BPA, internal/reaction force on
% the fixed block R = n*F_BPA*u_bpa + F_out*u_t (NOT assumed equal to the
% output tension), transmitted tendon force F_out = n*F_BPA/G (ideal),
% travel relations (ideal: F_out = n*F_BPA/G, tendon travel = G x BPA
% contraction), with the assumptions stated on the figure.
% Panel (B): the actual vas-med route of record at -120 deg and at
% +10 deg (from the route builder): BPA-side wrap (femur) vs cable-side
% ring (tibia), showing WHERE the BPA bends (Xi3) vs the cable bends.
%
% Style: Arial >= 10 pt at final size, bold margin panel letters (A)/(B)
% left of each tile straddling the plot-box top (FlxPin_group
% convention), CVD-safe (Tol bright palette), no color-only encoding.
% Export: FIG (editable source), PDF, EPS, PNG into this folder.

here = fileparts(mfilename('fullpath'));
R = 'D:\GitHub\Bipedal_Robot';
addpath(genpath(fullfile(R, 'Code', 'Matlab')));
addpath(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization'));

% Tol bright palette (Colors.m order), CVD-safe + distinct markers
col = [214, 88, 108; 0, 110, 191; 247, 143, 21; 0, 158, 137; 66, 74, 93] / 255;
colBPA = col(1, :); colCab = col(2, :); colGray = col(5, :);

fig = figure('Units', 'centimeters', 'Position', [2, 2, 17.8, 8.6], ...
    'Color', 'w', 'PaperPositionMode', 'auto');
axA = axes('Parent', fig, 'Position', [0.075, 0.16, 0.36, 0.72]);
hold(axA, 'on'); axis(axA, 'off'); axis(axA, 'equal');
xlim(axA, [0, 10]); ylim(axA, [0, 10]);

%% ---------------- Panel A: tackle free-body diagram -------------------
% Fixed block (proximal mount, top); moving block (bottom) tied to the
% bladder end; G = 2 line parts drawn (2:1); n = 2 BPAs in parallel.
xb = [3.4, 6.6]; yb = 8.9;                       % fixed block
xm = [4.3, 6.9]; ym = 3.1;                       % moving block
rectangle(axA, 'Position', [xb(1), yb - 0.45, xb(2) - xb(1), 0.9], ...
    'Curvature', 0.4, 'EdgeColor', colGray, 'LineWidth', 1.4, 'FaceColor', [1 1 1]);
rectangle(axA, 'Position', [xm(1), ym - 0.4, xm(2) - xm(1), 0.8], ...
    'Curvature', 0.4, 'EdgeColor', colGray, 'LineWidth', 1.4, 'FaceColor', [1 1 1]);
text(axA, mean(xb), yb + 0.75, 'fixed block (proximal mount)', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Arial');
text(axA, mean(xm), ym - 0.95, 'moving block (bladder end)', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontName', 'Arial');

% BPAs (n = 2, parallel, from origin bar to moving block)
for k = 1:2
    x0 = 0.6 + (k - 1) * 0.95;
    x1 = 5.05 + (k - 1) * 0.55;
    rectangle(axA, 'Position', [x0, 4.4, x1 - x0, 0.85], 'Curvature', 0.5, ...
        'EdgeColor', colBPA, 'LineWidth', 1.4, 'FaceColor', [1 0.93 0.93]);
    text(axA, (x0 + x1) / 2, 4.82, sprintf('BPA_%d', k), 'HorizontalAlignment', ...
        'center', 'FontSize', 10, 'FontName', 'Arial', 'Color', colBPA);
    plot(axA, [x0 - 0.7, x0], [4.825, 4.825], '-', 'Color', colGray, 'LineWidth', 1.2);
    plot(axA, [x1, x1], [5.25, ym], '-', 'Color', colBPA, 'LineWidth', 1.6);
end
plot(axA, [0.55, 1.6], [4.05, 4.05], '-', 'Color', colGray, 'LineWidth', 2);
text(axA, 0.35, 3.55, 'origin bracket (X_i1, X_i2)', 'FontSize', 10, ...
    'FontName', 'Arial', 'Color', colGray);

% Tackle rope: G = 2 line parts (down, around moving sheave, back up)
xrope1 = 5.6; xrope2 = 6.4;
plot(axA, [xrope1, xrope1], [yb - 0.45, ym + 0.4], '-', ...
    'Color', colCab, 'LineWidth', 1.6);
plot(axA, [xrope2, xrope2], [yb - 0.45, ym + 0.4], '-', ...
    'Color', colCab, 'LineWidth', 1.6);
viscirclesA(axA, [xrope1, ym + 0.4], 0.28, colCab);
viscirclesA(axA, [xrope2, yb - 0.45], 0.28, colCab);
text(axA, 7.6, 6.4, {'G = 2 line parts', '(travel gain G)'}, 'FontSize', 10, ...
    'FontName', 'Arial', 'Color', colCab);

% Exit tendon to the insertion (crosses the joint)
plot(axA, [xm(2) + 0.0, 9.3], [ym + 0.0, 1.6], '-', 'Color', colCab, 'LineWidth', 2);
plot(axA, 9.3, 1.6, 'v', 'MarkerFaceColor', colCab, 'MarkerEdgeColor', colCab, ...
    'MarkerSize', 7);
text(axA, 9.35, 0.85, 'tibial insertion', 'FontSize', 10, 'FontName', 'Arial', ...
    'HorizontalAlignment', 'right');

% Force arrows
quiverA(axA, 1.05, 3.9, 0, 0.0, []);  % placeholder no-op
arrowA(axA, [5.05, 4.0], [5.05, 3.35], colBPA);   % bladder pull -> moving block
text(axA, 4.85, 3.62, 'n\cdotF_BPA', 'FontSize', 10, 'FontName', 'Arial', ...
    'Color', colBPA, 'HorizontalAlignment', 'right', 'Interpreter', 'tex');
arrowA(axA, [xrope1 - 0.15, yb + 0.6], [xrope1 - 0.15, yb + 1.5], colGray);
text(axA, xrope1 - 0.3, yb + 1.75, {'R = n\cdotF_BPA\cdotu_bpa + F_t\cdotu_t', ...
    '(mount reaction; NOT = F_t)'}, 'FontSize', 10, 'FontName', 'Arial', ...
    'Color', colGray, 'Interpreter', 'tex', 'HorizontalAlignment', 'left');
arrowA(axA, [8.0, 2.45], [8.9, 1.95], colCab);
text(axA, 8.05, 2.85, 'F_t = n\cdotF_BPA/G', 'FontSize', 10, 'FontName', ...
    'Arial', 'Color', colCab, 'Interpreter', 'tex');

% Travel relations box
text(axA, 0.1, 9.6, {'Ideal relations (lossless rope, rigid blocks):', ...
    'F_t = n\cdotF_BPA / G;  tendon travel = G \times BPA contraction;', ...
    'bracket share cb\cdotF per BPA (own bracket); \delta_t = F_t/kSpr.'}, ...
    'FontSize', 10, 'FontName', 'Arial', 'Color', 'k', 'Interpreter', 'tex', ...
    'VerticalAlignment', 'top');

%% ---------------- Panel B: the actual vas-med route -------------------
axB = axes('Parent', fig, 'Position', [0.575, 0.16, 0.385, 0.72]);
hold(axB, 'on'); axis(axB, 'off');
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'xFullBest', 'ctxCamBest');
x = C.xFullBest(:).'; ctx = C.ctxCamBest;
[~, ~, infoFlex] = buildDistalRingLocation20mm(x(1:3), x(4:6), x(8), ctx);
rawFlex = infoFlex.raw(:, :, 1);                 % -120 deg (full flexion)
iExt = size(infoFlex.raw, 3);
rawExt = infoFlex.raw(:, :, iExt);               % +10 deg (extension)

% Route in the femur frame at each pose
Pf = zeros(9, 3, 2);
for jj = 1:2
    raw = infoFlex.raw(:, :, (jj - 1) * (iExt - 1) + 1);
    for r = 1:9
        if r <= 5
            Pf(r, :, jj) = raw(r, :);
        else
            T2 = ctx.T_Pam(:, :, (jj - 1) * (iExt - 1) + 1) * ...
                ctx.T_ICR_t1(:, :, (jj - 1) * (iExt - 1) + 1);
            Pf(r, :, jj) = RowVecTrans(T2, raw(r, :));
        end
    end
end
mk = {'o', 's'};
lg = {'-120 deg (full flexion)', '+10 deg (full extension)'};
for jj = 1:2
    ii = (jj - 1) * (iExt - 1) + 1;
    plot(axB, Pf(1:5, 1, jj), Pf(1:5, 2, jj), '-', 'Color', col(jj, :), ...
        'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(axB, Pf(1, 1, jj), Pf(1, 2, jj), mk{jj}, 'MarkerFaceColor', ...
        col(jj, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
        'LineStyle', 'none', 'DisplayName', lg{jj});
end
text(axB, Pf(1, 1, 1), Pf(1, 2, 1) + 0.015, 'p1 origin (BPA)', 'FontSize', 10, ...
    'FontName', 'Arial');
text(axB, -0.05, -0.335, 'BPA bends (femur cylinder + condyle/clip): \chi_3 wrap loss', ...
    'FontSize', 10, 'FontName', 'Arial', 'Color', colBPA);
text(axB, -0.05, -0.395, 'Cable bends (tibia ring p6-p8): no \chi_3', ...
    'FontSize', 10, 'FontName', 'Arial', 'Color', colCab);
axis(axB, 'equal');
xl = [min(Pf(:, 1, :), [], 'all') - 0.05, max(Pf(:, 1, :), [], 'all') + 0.05];
yl = [min(Pf(:, 2, :), [], 'all') - 0.14, max(Pf(:, 2, :), [], 'all') + 0.06];
xlim(axB, xl); ylim(axB, yl);
xlabel(axB, 'x, m (femur frame)', 'FontSize', 10, 'FontName', 'Arial');
ylabel(axB, 'y, m', 'FontSize', 10, 'FontName', 'Arial');
legend(axB, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'FontSize', 9, 'FontName', 'Arial');

%% ---------------- panel letters (FlxPin_group convention) -------------
drawnow;
apos = axA.Position;
tcY = apos(2) + apos(4);
annotation(fig, 'textbox', [0.012, tcY - 0.03, 0.05, 0.07], 'String', '(A)', ...
    'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial', ...
    'EdgeColor', 'none', 'VerticalAlignment', 'bottom');
bpos = axB.Position;
tcY2 = bpos(2) + bpos(4);
annotation(fig, 'textbox', [bpos(1) - 0.055, tcY2 - 0.03, 0.05, 0.07], ...
    'String', '(B)', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial', ...
    'EdgeColor', 'none', 'VerticalAlignment', 'bottom');

set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig, '-property', 'FontSize'), 'FontSize', 10);

saveas(fig, fullfile(here, 'vasmed_reverse_pulley_schematic_20261028.fig'));
exportgraphics(fig, fullfile(here, 'vasmed_reverse_pulley_schematic_20261028.pdf'), ...
    'ContentType', 'vector');
exportgraphics(fig, fullfile(here, 'vasmed_reverse_pulley_schematic_20261028.png'), ...
    'Resolution', 300);
% EPS via print (exportgraphics has no EPS)
print(fig, fullfile(here, 'vasmed_reverse_pulley_schematic_20261028.eps'), ...
    '-depsc2', '-vector');
fprintf('RP5 schematic DONE\n');
end

function viscirclesA(ax, c, r, colr)
th = linspace(0, 2 * pi, 40);
plot(ax, c(1) + r * cos(th), c(2) + r * sin(th), '-', 'Color', colr, ...
    'LineWidth', 1.4);
plot(ax, c(1), c(2), '.', 'Color', colr, 'MarkerSize', 10);
end

function arrowA(ax, p0, p1, colr)
quiver(ax, p0(1), p0(2), p1(1) - p0(1), p1(2) - p0(2), 0, ...
    'Color', colr, 'LineWidth', 1.8, 'MaxHeadSize', 0.6);
end

function quiverA(ax, x, y, dx, dy, ~) %#ok<STOUT,INUSD>
if abs(dx) + abs(dy) > 0
    quiver(ax, x, y, dx, dy, 0, 'LineWidth', 1.4);
end
end
