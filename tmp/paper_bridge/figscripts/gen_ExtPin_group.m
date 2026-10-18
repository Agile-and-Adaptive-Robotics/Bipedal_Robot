%gen_ExtPin_group.m — regenerate ExtPin_group.eps for the Frontiers manuscript.
%3 panels (A)(B)(C) in the PlotPrompt ext staggered layout (tiledlayout(2,7),
%tiles [1 5 10] spanning [1 3]); panels show pinned-extensor tests
%ke(9)=48.0 cm, ke(7)=45.7 cm, ke(3)=41.5 cm (TRUE resting lengths, ledger §4 /
%chapter fig:ExtPin_group caption).
%Series: Hybrid gold, Measured indigo, Original gray dashed, Predicted magenta.
%Model: minimizeExtX3 at the INITIAL extensor terms (chi0 = -10.9 mm,
%chi3 = 0.123, stiffness pair shared with the flexor identification; plan §2 3.2).
%YLim [0 15]. Writes ONLY the deliverable EPS (+ scratch QC PNG).
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('EXTX3_TRANS', '');   % default transMode (2trans)

c1 = [1 0.843 0];        % #FFD700 gold    - Hybrid
c5 = [0.804 0.204 0.710];% #CD34B5 magenta - Predicted
c7 = [0 0 1];            % #0000FF indigo  - Measured
cgray = [0.4 0.4 0.4];   % Original dashed
sz = 60;

pmap = [9 7 3];                       % panels A/B/C -> current ke indices
[fE, bpaE] = minimizeExtX3(-0.0109, 1.02e5, 1.06e4, 0.123, pmap); %#ok<ASGLU>
for k = 1:3
    fprintf('ke(%d) rest=%.3f: GoF %.4f/%.4f/%.4f\n', pmap(k), bpaE(pmap(k)).rest, ...
        fE(k,1), fE(k,2), fE(k,3));
end

tileIdxs = [1, 5, 10];     % PlotPrompt ext layout
tileSpans = [1 3];
titles = {'\bf 48.0 cm', '\bf 45.7 cm', '\bf 41.5 cm'};
xAnn = [0.035, 0.51, 0.265];
yAnn = [0.89, 0.89, 0.41];
tileLabels = {'(A)', '(B)', '(C)'};

figT = figure('Name','ExtPin_group','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(2,7,'TileSpacing','loose','Padding','loose');

axPanels = gobjects(1,3);
for j = 1:3
    i = pmap(j);
    axPanels(j) = nexttile(tileIdxs(j), tileSpans);
    hold on
    scatter(bpaE(i).A_h, bpaE(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName', 'Hybrid');
    scatter(bpaE(i).Aexp, bpaE(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
        'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName', 'Measured');
    plot(bpaE(i).Ak, bpaE(i).M, '--', 'Color', cgray, 'LineWidth', 2, ...
        'DisplayName', 'Original');
    plot(bpaE(i).Ak, bpaE(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, ...
        'DisplayName', 'Predicted');
    title(titles{j}, 'Interpreter', 'tex');
    annotation(figT, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'XTick', [-120 -90 -60 -30 0], 'YLim', [0 15], ...
        'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
end

% shared axes labels
ylabel(tT, '\bf Torque, N\cdotm', 'Interpreter', 'tex');
xlabel(tT, '\bf \theta_{k} , \circ', 'Interpreter', 'tex');

% legend on panel (B), northeast, 8 pt (PlotPrompt ext convention)
lg = legend(axPanels(2)); lg.Location = 'northeast'; lg.FontSize = 8;

exportgraphics(figT, fullfile(outdir, 'ExtPin_group.eps'), 'ContentType', 'vector');
exportgraphics(figT, fullfile(qcdir, 'QC_ExtPin_group.png'), 'Resolution', 150);
fprintf('WROTE %s\n', fullfile(outdir, 'ExtPin_group.eps'));
close(figT);
disp('DONE_ExtPin_group');
