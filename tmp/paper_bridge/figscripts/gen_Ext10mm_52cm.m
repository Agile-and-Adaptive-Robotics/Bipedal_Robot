%gen_Ext10mm_52cm.m — regenerate Ext10mm_52cm.eps for the Frontiers manuscript.
%Single panel, biomimetic-knee extensor, l0 = 51.8 cm (filename keeps the legacy
%52 cm name per plan §6). Four series (plan spec / chapter fig:Ext10mm_52cm):
%  Measured indigo scatter, Hybrid gold scatter, Original gray dashed,
%  Optimized prediction magenta 2.5.
%Producer path (located): minimizeExt.m — the bio-extensor evaluator that loads
%KneeExt_10mm_52cm.mat + Plot_KneeExt10mm_52cm.mat; corrected series evaluated at
%the INITIAL extensor terms (chi0 = -10.9 mm, chi3 = 0.123 on the shared
%stiffness pair 1.02e5 / 1.06e4 N/m).
%Writes ONLY the deliverable EPS (+ scratch QC PNG).
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('EXTX3_TRANS', '');

c1 = [1 0.843 0];        % #FFD700 gold    - Hybrid
c5 = [0.804 0.204 0.710];% #CD34B5 magenta - Optimized prediction
c7 = [0 0 1];            % #0000FF indigo  - Measured
cgray = [0.4 0.4 0.4];   % Original dashed
sz = 60;

[fB, bpaB] = minimizeExt(-0.0109, 1.02e5, 1.06e4, 0.123, 1); %#ok<ASGLU>
fprintf('bio ext GoF (current evaluator, initial coeffs): %.4f/%.4f/%.4f  (chapter initial-era: 0.8920/0.3396/1.4188)\n', ...
    fB(1,1), fB(1,2), fB(1,3));
b = bpaB(1);

figT = figure('Name','Ext10mm_52cm','Color','w');
figT.Position = [100 100 760 560];
ax = axes(figT); hold(ax, 'on');

scatter(b.A_h, b.M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName', 'Hybrid');
scatter(b.Aexp, b.Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName', 'Measured');
plot(b.Ak, b.M, '--', 'Color', cgray, 'LineWidth', 2, 'DisplayName', 'Original');
plot(b.Ak, b.M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');

title('\bf l_{0} = 51.8 cm', 'Interpreter', 'tex');
xlabel('\bf \theta_{k} , \circ', 'Interpreter', 'tex');
ylabel('\bf Torque, N\cdotm', 'Interpreter', 'tex');
set(ax, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
    'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'XLim', [-125 15], 'XTick', [-120 -90 -60 -30 0], 'YLim', [0 20], ...
    'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
lg = legend(ax); lg.Location = 'northeast'; lg.FontSize = 8;

exportgraphics(figT, fullfile(outdir, 'Ext10mm_52cm.eps'), 'ContentType', 'vector');
exportgraphics(figT, fullfile(qcdir, 'QC_Ext10mm_52cm.png'), 'Resolution', 150);
fprintf('WROTE %s\n', fullfile(outdir, 'Ext10mm_52cm.eps'));
close(figT);
disp('DONE_Ext10mm_52cm');
