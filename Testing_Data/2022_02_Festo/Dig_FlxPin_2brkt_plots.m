%Dig_FlxPin_2brkt_plots.m — standard 2brk plots for ANY pick, from the saved mat.
%Set PICK (row in filtered_results, sorted by validation distance), Run, figures appear.
%In interactive MATLAB the figures open live; in -batch they save to Dig_out\ as .png/.fig.
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');

MATNAME = 'minimizeFlxPin10_results_20260907_2brkt_2trans.mat';
PICK    = 168;

S = load(MATNAME, 'filtered_results', 'xCols', 'ALLBPA', 'labels');
F = S.filtered_results; xC = S.xCols;
valCols = 10:12; distCol = 13;   %[ind, hold(2), x(3), train(3), val(3), dist]
ALLBPA = S.ALLBPA; labels = string(S.labels); numBPA = numel(ALLBPA);
p = min(max(PICK, 1), size(F,1));
k1 = F(p,xC(1)); k2 = F(p,xC(2)); k3 = F(p,xC(3));
fprintf('pick %d: Xi0=%.4f m  Xi1=%.3g  Xi2=%.3g  | valRMSE=%.3f valFVU=%.3f dist=%.3f\n', ...
    p, k1, k2, k3, F(p,valCols(1)), F(p,valCols(2)), F(p,distCol));

[f, bpa] = minimizeFlxPin2brk(k1, k2, k3, [], true, '2trans');
disp(array2table(f, 'VariableNames', {'RMSE','FVU','MaxResidual'}, 'RowNames', cellstr(labels')));

%% Plot torque curves, pre- and post-Optimized
load ForceStrainForFit.mat z
z(isnan(z)) = 0;
X = linspace(0,620,20); X = X(2:20);
Y = linspace(0,1,30);

figTpre = figure('Name','Torque, Pre-Optimized','Color','w');
figTpre.Position = [100 100 950 700];
tTpre = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');
titles = ["\bf 48.5 cm", "\bf 45.7 cm","\bf 47.9 cm", "\bf 40.6 cm", "\bf 41.7 cm"];

for k = 1:numBPA
    nexttile(k); j = ALLBPA(k); hold on
    Yq = bpa(j).strain./((bpa(j).rest - bpa(j).Kmax)/bpa(j).rest);
    Xq = bpa(j).P.*ones(size(Yq));
    Vq = interp2(X, Y', z, Xq, Yq,'spline');
    Fold = Vq.*bpa(j).unitD;
    Fq = hypot(Fold(:,1), Fold(:,2));
    Mold = -bpa(j).mA.*Fq;
    Mold(Yq>1) = 0;
    Mold(bpa(j).strain < 0) = NaN;
    scatter(bpa(j).A_h, bpa(j).M_h, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#FFD700', 'DisplayName', 'Hybrid');
    scatter(bpa(j).Aexp, bpa(j).Mexp, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Mold, '--', 'Color', '#000000', 'LineWidth', 2, 'DisplayName', 'Original');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    clear Vq Fold Fq Mold
    title(titles(k), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    xlim([-120 20]); ylim([-25 0]);
end
ylabel(tTpre,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpre,'\bf \theta_{k} , \circ','Interpreter','tex');
lg = legend(tTpre.Children(end-1)); lg.Location = 'best'; lg.FontSize = 8;

figTpost = figure('Name','Torque, Post-Optimized','Color','w');
figTpost.Position = [100 100 950 700];
tTpost = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');
for k = 1:numBPA
    nexttile(k); j = ALLBPA(k); hold on
    scatter(bpa(j).Aexp, bpa(j).Mexp, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).M, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).M_p(:,3), '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    title(titles(k), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
    xlim([-120 20]); ylim([-25 0]);
end
ylabel(tTpost,'\bf Torque, N\cdotm','Interpreter','tex');
xlabel(tTpost,'\bf \theta_{k} , \circ','Interpreter','tex');
lg2 = legend(tTpost.Children(end-1)); lg2.Location = 'best'; lg2.FontSize = 8;

%% Muscle length
figL = figure('Name','Muscle Length','Color','w');
figL.Position = [100 100 950 700];
tL = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');
for k = 1:numBPA
    nexttile(k); j = ALLBPA(k); hold on
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    Lm_p = bpa(j).Lmt_p - 2*bpa(j).fitn - bpa(j).ten;
    Lm   = bpa(j).Lmt   - 2*bpa(j).fitn - bpa(j).ten;
    scatter(bpa(j).A_h, bpa(j).Lm_h, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, Lm, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Original prediction');
    plot(bpa(j).Ak, Lm_p, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end
ylabel(tL,'\bf Muscle length, m','Interpreter','tex');
xlabel(tL,'\bf \theta_{k} , \circ','Interpreter','tex');
lg = legend(tL.Children(end-1)); lg.Location = 'best'; lg.FontSize = 8;

%% Moment arm
figMA = figure('Name','Moment Arm','Color','w');
figMA.Position = [100 100 950 700];
tMA = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');
for k = 1:numBPA
    nexttile(k); j = ALLBPA(k); hold on
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    G_p = hypot(bpa(j).mA_p(:,1), bpa(j).mA_p(:,2));
    scatter(bpa(j).A_h, bpa(j).mA_h, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).mA, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, G_p, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end
ylabel(tMA,'\bf Moment arm, m','Interpreter','tex');
xlabel(tMA,'\bf \theta_{k} , \circ','Interpreter','tex');
legend(tMA.Children(end-1),'Location','best','FontSize',8);

%% Relative strain
figS = figure('Name','Relative Strain','Color','w');
figS.Position = [100 100 950 700];
tS = tiledlayout(ceil(numBPA/2),2,'TileSpacing','loose','Padding','loose');
for k = 1:numBPA
    nexttile(k); j = ALLBPA(k); hold on
    if bpa(j).ten > 0
        title(sprintf('\\bf l_0 = %0.1f cm, %0.0f mm tendon',bpa(j).rest*100, bpa(j).ten*10^3),'Interpreter','tex')
    else
        title(sprintf('\\bf l_0 = %0.1f cm',bpa(j).rest*100),'Interpreter','tex')
    end
    strain_h = (bpa(j).rest - bpa(j).Lm_h)/bpa(j).rest;
    kmax = (bpa(j).rest - bpa(j).Kmax)/bpa(j).rest;
    scatter(bpa(j).A_h, strain_h/kmax, 60, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerFaceColor', '#0000FF', 'DisplayName', 'Measured');
    plot(bpa(j).Ak, bpa(j).strain/kmax, '-.', 'Color', '#FA8775', 'LineWidth', 2.5, 'DisplayName', 'Improved BPA model');
    plot(bpa(j).Ak, bpa(j).strain_p/kmax, '-', 'Color', '#CD34B5', 'LineWidth', 2.5, 'DisplayName', 'Optimized prediction');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
        'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', 'TickLength', [0.025 0.05]);
end
ylabel(tS,'\bf Relative strain','Interpreter','tex');
xlabel(tS,'\bf \theta_{k} , \circ','Interpreter','tex');
legend(tS.Children(end-1),'Location','best','FontSize',8);

%% Save if batch, otherwise leave figures open
if batchStartupOptionUsed
    out = 'Dig_out';
    figs = [figTpre, figTpost, figL, figMA, figS];
    names = {'torquePre','torquePost','muscleLength','momentArm','relStrain'};
    for q = 1:numel(figs)
        exportgraphics(figs(q), fullfile(out, sprintf('Dig_FlxPin_2brkt_pick%d_%s.png', p, names{q})), 'Resolution', 150);
        savefig(figs(q), fullfile(out, sprintf('Dig_FlxPin_2brkt_pick%d_%s.fig', p, names{q})));
        close(figs(q));
    end
    fprintf('Figures saved to Dig_out\\Dig_FlxPin_2brkt_pick%d_*.png/.fig\n', p);
end
