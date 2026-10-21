%rerun_ExtBiomimetic.m — refresh the 10mm/52cm biomimetic EXTENSOR chain
%per Ben 2026-10-19: rerun Knee_Extensor_10mm (now at the adopted row-32
%X3 terms; routing unchanged), rewrite KneeExt_10mm_52cm.mat (old mat
%archived), rerun Plot_KneeExt10mm_52cm, print GoF triplets, and rebuild
%the paper's single-panel Ext10mm_52cm figure in campaign style.
root = 'D:/GitHub/Bipedal_Robot';
resdir = fullfile(root, 'Code/Matlab/Mesh_Optimization/Results');
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
addpath(fullfile(root, 'Code/Matlab/Human_Data'));
addpath(genpath(fullfile(root, 'Code/Matlab/Mesh_Optimization')));
warning('off', 'all');

%% 1. data script -> fresh workspace
cd(resdir);
matpath = fullfile(festodir, 'KneeExt_10mm_52cm.mat');
archive = fullfile(festodir, 'KneeExt_10mm_52cm_stale_20261019.mat');
if exist(matpath, 'file') && ~exist(archive, 'file')
    copyfile(matpath, archive);
    fprintf('ARCHIVED old mat -> %s\n', archive);
end
run('Knee_Extensor_10mm.m');
% NOTE: the data script starts with `clear`, wiping this driver's variables;
% re-derive every path from literals before using them again.
matpath = fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo', 'KneeExt_10mm_52cm.mat');
outdir = fullfile('D:/GitHub/Bipedal_Robot', 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile('D:/GitHub/Bipedal_Robot', 'tmp/paper_bridge/figscripts');
fprintf('AUDIT Xi used: Xi0=%.6g Xi1=%.6g Xi2=%.6g Xi3=%.6g wraps=%d bpaCount=%d\n', ...
    Xi0, Xi1, Xi2, Xi3, wraps, BPAcount);
save(matpath);
fprintf('WROTE %s\n', matpath);
close all;

%% 2. Ben's plot function (regenerates TorqueHand from ExtTest stats)
cd(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo'));
addpath(fullfile(pwd, 'Ext_10mm'));   % transient: ExtTest stats mats
run('Plot_KneeExt10mm_52cm.m');
rmpath(fullfile('D:/GitHub/Bipedal_Robot', 'Testing_Data/2022_02_Festo/Ext_10mm'));
close all;
% the Plot script also clears the workspace: re-derive export paths
outdir = fullfile('D:/GitHub/Bipedal_Robot', 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile('D:/GitHub/Bipedal_Robot', 'tmp/paper_bridge/figscripts');

%% 3. GoF triplets at the measured angles
ph = phiD(:); th = Theoretical(:); tc = Vas_Pam_adj.Torque_p(:,3);
th(th ~= th) = NaN; tc(tc ~= tc) = NaN;
G = griddedInterpolant(ph, th); Mc = G(Angle(:));
G2 = griddedInterpolant(ph, tc); Cc = G2(Angle(:));
Me = Torque(:); %#ok<NODEF> workspace from Plot script
gof = @(p, m) deal(sqrt(mean((p(:)-m(:)).^2)), ...
    sum((p(:)-m(:)).^2)/sum((m(:)-mean(m(:))).^2), max(abs(p(:)-m(:))));
[oR, oF, oM] = gof(Mc, Me);
[cR, cF, cM] = gof(Cc, Me);
[hR, hF, hM] = gof(TorqueHand(:), Me);
fprintf('AUDIT GoF ext-10mm-52cm vs measured [RMSE FVU Max]:\n');
fprintf('  Original (rigid):        %8.4f %8.4f %8.4f\n', oR, oF, oM);
fprintf('  Corrected (row-32 Xi):   %8.4f %8.4f %8.4f\n', cR, cF, cM);
fprintf('  Hybrid (back-calc):      %8.4f %8.4f %8.4f\n', hR, hF, hM);

%% 4. campaign-style single-panel figure
c1 = [1 0.843 0]; c5 = [0.804 0.204 0.710]; c7 = [0 0 1]; cgray = [0.4 0.4 0.4];
sz = 60;
fig = figure('Name','Ext10mm_52cm','Color','w'); fig.Position = [90 90 700 480];
hold on
scatter(Angle, Torque, sz, 'd', 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName','Measured');
scatter(Angle, TorqueHand, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
    'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName','Hybrid');
plot(ph, th, '--', 'Color', cgray, 'LineWidth', 2, 'DisplayName','Original');
plot(ph, tc, '-', 'Color', c5, 'LineWidth', 2.5, 'DisplayName','Corrected');
set(gca, 'FontSize',12,'FontWeight','bold','FontName','Arial','LineWidth',2, ...
    'XMinorTick','on','YMinorTick','on','TickLength',[0.025 0.05], ...
    'XLim',[-120 30],'XTick',[-120 -90 -60 -30 0 30], ...
    'YLim',[0 ceil(max([Me(:); TorqueHand(:); th(~isnan(th)); tc(~isnan(tc))])/5)*5], ...
    'GridLineStyle','none');
xlabel('\bf \theta_{k} , \circ', 'Interpreter','tex');
ylabel('\bf Torque, N\cdotm', 'Interpreter','tex');
lg = legend('Location','southwest'); lg.FontSize = 8;
exportgraphics(fig, fullfile(outdir, 'Ext10mm_52cm.eps'), 'ContentType','vector');
exportgraphics(fig, fullfile(qcdir, 'QC_ExtBiomimetic10mm.png'), 'Resolution',150);
fprintf('WROTE Ext10mm_52cm.eps (adopted row-32 terms)\n');
close(fig);
disp('DONE_ExtBiomimetic');
