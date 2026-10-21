%gen_ExtBiomimetic10mm.m — rebuild Ext10mm_52cm.eps from minimizeExt (the
%biomimetic extensor evaluator, Ben 2026-10-19 ruling: minimizeFlx feeds the
%flexor biomimetic plots, minimizeExt feeds this one; the X3 bend model with
%R1=0.022/28deg + R2=0.176 lives inside its Contraction). One pipeline:
%zero-correction baseline, corrected at adopted row-32 terms, corrected at
%initial terms (for the prose), measured, hybrid — all from the evaluator.
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');

[f0, b0] = minimizeExt(0, Inf, Inf, 0);                            % zero correction
[f1, b1] = minimizeExt(-6.378022e-3, 3.997873e4, 1.473375e4, 1.582887e-1);  % adopted row 32
[f2, b2] = minimizeExt(-0.0109, 1.02e5, 1.06e4, 0.123);            % initial terms
fprintf('AUDIT minimizeExt GoF [RMSE FVU Max]:\n');
fprintf('  zero-correction : %8.4f %8.4f %8.4f\n', f0(1,:));
fprintf('  adopted (row32) : %8.4f %8.4f %8.4f\n', f1(1,:));
fprintf('  initial terms   : %8.4f %8.4f %8.4f\n', f2(1,:));
fprintf('  hybrid (M_h)    : %8.4f %8.4f %8.4f\n', ...
    sqrt(mean((b1(1).M_h(:)-b1(1).Mexp(:)).^2)), ...
    sum((b1(1).M_h(:)-b1(1).Mexp(:)).^2)/sum((b1(1).Mexp(:)-mean(b1(1).Mexp(:))).^2), ...
    max(abs(b1(1).M_h(:)-b1(1).Mexp(:))));
fprintf('  series check max|b0.M - b1.M| = %.4f (expect ~0, terms-independent)\n', ...
    max(abs(b0(1).M(:) - b1(1).M(:))));

c1 = [1 0.843 0]; c5 = [0.804 0.204 0.710]; c7 = [0 0 1]; cgray = [0.4 0.4 0.4];
sz = 60;
fig = figure('Name','Ext10mm_52cm','Color','w'); fig.Position = [90 90 700 480];
hold on
% series order per minimizeExt10mm.m figT: Hybrid, Measured, Original, Predicted
scatter(b1(1).A_h, b1(1).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName','Hybrid');
scatter(b1(1).Aexp, b1(1).Mexp, sz, 'd', 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
plot(b1(1).Ak, b1(1).M, '--', 'Color', cgray, 'LineWidth', 2, 'DisplayName','Original');
plot(b1(1).Ak, b1(1).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Predicted');
set(gca, 'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 30],'XTick',[-120 -90 -60 -30 0 30], 'YLim',[0 15], ...
    'GridLineStyle','none');
xlabel('\bf \theta_{k} , \circ', 'Interpreter','tex');
ylabel('\bf Torque, N\cdotm', 'Interpreter','tex');
lg = legend('Location','southwest'); lg.FontSize = 8;
exportgraphics(fig, fullfile(outdir, 'Ext10mm_52cm.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(qcdir, 'QC_ExtBiomimetic10mm.png'), 'Resolution',150);
fprintf('WROTE Ext10mm_52cm.eps (minimizeExt pipeline)\n');
close(fig);
disp('DONE_genExtBio');
