% gen_RouteVariants_20261020.m -- EDIT_SPEC section-5 figure 3: the four
% route variants of Table 5 (tab:Redesign) as curves instead of a table.
%   (A) extensor: direct route (Vas_Pam_20mm_Result_20260925) vs cam design
%       of record (Vas_Pam_20mm_Result_cam_20261014_0416) vs the vasti
%       human requirement (+5% envelope);
%   (B) flexor: direct route (Bifemsh_20mm_Result_20260921_2037) vs the
%       mono-pulley insertion-side torque (Bifemsh_20mm_Result_pulley_
%       20260926_1652, TorqueInsZ) vs the bifemsh short-head requirement
%       (+5% envelope).
% Campaign style per EDIT_SPEC section 0. EPS -> ProofFinal/Figures,
% QC PNG (150 dpi) -> figscripts.
root = 'D:/GitHub/Bipedal_Robot';
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');

c1 = [1 0.843 0];          % gold        - mono-pulley flexor variant
c3 = [0.980 0.529 0.459];  % light orange - direct routes (superseded)
c5 = [0.804 0.204 0.710];  % magenta     - cam extensor (design of record)
c8 = [0 0 0];
cgray = [0.4 0.4 0.4];

sty = {'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
    'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', [0.025 0.05], 'GridLineStyle', 'none', ...
    'XLim', [-120 30], 'XTick', [-120 -90 -60 -30 0 30]};

V = load(fullfile(resDir, 'Vas_Pam_20mm_Result_20260925.mat'), 'predBest', 'ctx');
C = load(fullfile(resDir, 'Vas_Pam_20mm_Result_cam_20261014_0416.mat'), 'predBest', 'ctx');
B = load(fullfile(resDir, 'Bifemsh_20mm_Result_20260921_2037.mat'), 'predBest', 'ctx');
M = load(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_20260926_1652.mat'), 'predBest', 'ctx');

fig = figure('Name', 'RouteVariants', 'Color', 'w');
fig.Position = [80 60 1120 470];
tl = tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'loose');
ax = gobjects(1, 2);

nexttile(1); hold on
plot(V.ctx.phiD, V.predBest.TorqueZ, '-', 'Color', c3, 'LineWidth', 2, ...
    'DisplayName', 'Direct extensor route');
plot(C.ctx.phiD, C.predBest.TorqueZ, '--', 'Color', c5, 'LineWidth', 2.5, ...
    'DisplayName', 'Cam extensor (design of record)');
plot(V.ctx.phiD, 1.05 * V.ctx.humanTorque, ':', 'Color', cgray, 'LineWidth', 2.5, ...
    'DisplayName', 'Required envelope (+5%)');
plot(V.ctx.phiD, V.ctx.humanTorque, '-.', 'Color', c8, 'LineWidth', 2, ...
    'DisplayName', 'Human target (vasti)');
title('\bf Extensor route variants'); set(gca, sty{:});
ax(1) = gca;
lg1 = legend('Location', 'southwest'); lg1.FontSize = 8;

nexttile(2); hold on
plot(B.ctx.phiD, B.predBest.TorqueZ, '-', 'Color', c3, 'LineWidth', 2, ...
    'DisplayName', 'Direct flexor route');
plot(M.ctx.phiD, M.predBest.TorqueInsZ, '--', 'Color', c1, 'LineWidth', 2.5, ...
    'DisplayName', 'Mono-pulley flexor (insertion side)');
plot(B.ctx.phiD, 1.05 * B.ctx.humanTorque, ':', 'Color', cgray, 'LineWidth', 2.5, ...
    'DisplayName', 'Required envelope (+5%)');
plot(B.ctx.phiD, B.ctx.humanTorque, '-.', 'Color', c8, 'LineWidth', 2, ...
    'DisplayName', 'Human target (bifemsh)');
title('\bf Flexor route variants'); set(gca, sty{:});
ax(2) = gca;
lg2 = legend('Location', 'northeast'); lg2.FontSize = 8;

% Y limits: ceil-to-5 per section 0
y1 = [min([V.predBest.TorqueZ; C.predBest.TorqueZ; 1.05 * V.ctx.humanTorque]), ...
    max([V.predBest.TorqueZ; C.predBest.TorqueZ])];
ax(1).YLim = [floor(y1(1) / 5) * 5, ceil(y1(2) / 5) * 5];
y2 = [min([B.predBest.TorqueZ; M.predBest.TorqueInsZ]), 5];
ax(2).YLim = [floor(y2(1) / 5) * 5, ceil(y2(2) / 5) * 5];

ylabel(tl, '\bf Torque, N\cdotm'); xlabel(tl, '\bf \theta_{k} , \circ');
set(tl, 'TileSpacing', 'loose');
drawnow;
xp2 = [0.005 0.505];
for j = 1:2
    axj = ax(j); p4 = axj.Position;
    tcY = 0.935;
    try
        Tj = axj.Title; Tj.Units = 'normalized';
        te = Tj.Extent;
        cand = p4(2) + (te(2) + te(4) / 2) * p4(4);
        if isfinite(cand) && cand > 0.3, tcY = cand; end
    catch
    end
    a1 = annotation(fig, 'textbox', [xp2(j) tcY - 0.025 0.05 0.05], ...
        'String', sprintf('\\bf (%c)', 'A' + j - 1), 'FontSize', 12, ...
        'FontName', 'Arial', 'FontWeight', 'bold', 'EdgeColor', 'none');
    try
        drawnow; ae = a1.Extent;
        a1.Position(2) = a1.Position(2) + (tcY - (ae(2) + ae(4) / 2));
    catch
    end
end
exportgraphics(fig, fullfile(outdir, 'RouteVariants_20261020.eps'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(qcdir, 'QC_RouteVariants_20261020.png'), 'Resolution', 150);
fprintf('WROTE RouteVariants_20261020.eps\n'); close(fig);
