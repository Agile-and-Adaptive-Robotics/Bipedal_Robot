%gen_FlxPin_group.m — regenerate FlxPin_group.eps for the Frontiers manuscript.
%2x2 tiles: columns Optimization (48.5 cm) / Validation (45.7 cm);
%rows Pre-optimization / Optimized; legends in tiles 2 and 4.
%Series (plan §6 / PlotPrompt / Dig_FlxPin_2brkt palette):
%  row 1: Hybrid gold, Measured indigo, Original gray dashed (old force surface),
%         Improved BPA model light orange dash-dot
%  row 2: Measured indigo, Improved BPA model light orange dash-dot,
%         Optimized prediction magenta 2.5
%XLim [-120 20], YLim [-35 0], XTick [-120 -90 -60 -30 0] (PlotPrompt final para).
%Writes ONLY the deliverable EPS (+ a scratch PNG for QC in tmp/paper_bridge).
root = 'D:/GitHub/Bipedal_Robot';
festodir = fullfile(root, 'Testing_Data/2022_02_Festo');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');
addpath(fullfile(root, 'Code/Matlab/Functions'));
addpath(fullfile(root, 'Code/Matlab/Functions/ModernRobotics'));
addpath(fullfile(root, 'Code/Matlab/Robot_Data'));
cd(festodir);
warning('off', 'all');
setenv('FLX_T5_YMM', '');   % ensure off (default)

% palette (PlotPrompt c{1..8})
c1 = [1 0.843 0];        % #FFD700 gold    - Hybrid
c3 = [0.980 0.529 0.459];% #FA8775 lt org  - Improved BPA model
c5 = [0.804 0.204 0.710];% #CD34B5 magenta - Optimized prediction
c7 = [0 0 1];            % #0000FF indigo  - Measured
cgray = [0.4 0.4 0.4];   % Original dashed
sz = 60;

% Model of record: two-bracket flexor evaluator at the INITIAL identification
% coefficients (chi0 = 11.9 mm, chi1 = 102 kN/m, chi2 = 10.6 kN/m; plan §2 3.1)
[fF, bpaF] = minimizeFlxPin(0.0119, 1.02e5, 1.06e4, [1 2]); %#ok<ASGLU>
fprintf('GoF (current evaluator, initial coeffs): test1 %.4f/%.4f/%.4f  test2 %.4f/%.4f/%.4f\n', ...
    fF(1,1), fF(1,2), fF(1,3), fF(2,1), fF(2,2), fF(2,3));

% "Original" series: old force-surface prediction (ForceStrainForFit z),
% per minimizeFlxPin10mm.m / Dig_FlxPin_2brkt_plots.m
load('ForceStrainForFit.mat', 'z');
z(isnan(z)) = 0;
Xp = linspace(0,620,20); Xp = Xp(2:20);
Yp = linspace(0,1,30);
Mold = cell(2,1);
for i = 1:2
    Yq = bpaF(i).strain./((bpaF(i).rest - bpaF(i).Kmax)/bpaF(i).rest);
    Xq = bpaF(i).P.*ones(size(Yq));
    Vq = interp2(Xp, Yp', z, Xq, Yq, 'spline');
    Fold = Vq.*bpaF(i).unitD;
    Fq = hypot(Fold(:,1), Fold(:,2));
    Mo = -bpaF(i).mA.*Fq;
    Mo(Yq>1) = 0;
    Mo(bpaF(i).strain < 0) = NaN;
    Mold{i} = Mo;
end

figT = figure('Name','FlxPin_group','Color','w');
figT.Position = [100 100 950 700];
tT = tiledlayout(2,2,'TileSpacing','loose','Padding','loose');

colTitle = {'\bf Optimization, l_{0} = 48.5 cm', '\bf Validation, l_{0} = 45.7 cm'};
axTiles = gobjects(2,2);
for row = 1:2
    i = row;   % test 1 = Optimization column, test 2 = Validation column
    for col = 1:2
        axTiles(row,col) = nexttile((row-1)*2 + col); hold on
        if row == 1
            scatter(bpaF(i).A_h, bpaF(i).M_h, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
                'MarkerFaceColor', c1, 'MarkerEdgeColor', 'none', 'DisplayName', 'Hybrid');
            scatter(bpaF(i).Aexp, bpaF(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
                'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName', 'Measured');
            plot(bpaF(i).Ak, Mold{i}, '--', 'Color', cgray, 'LineWidth', 2, ...
                'DisplayName', 'Original');
            plot(bpaF(i).Ak, bpaF(i).M, '-.', 'Color', c3, 'LineWidth', 2.5, ...
                'DisplayName', 'Improved BPA model');
        else
            scatter(bpaF(i).Aexp, bpaF(i).Mexp, sz, 'filled', 'MarkerFaceAlpha', 0.75, ...
                'MarkerFaceColor', c7, 'MarkerEdgeColor', 'none', 'DisplayName', 'Measured');
            plot(bpaF(i).Ak, bpaF(i).M, '-.', 'Color', c3, 'LineWidth', 2.5, ...
                'DisplayName', 'Improved BPA model');
            plot(bpaF(i).Ak, bpaF(i).M_p(:,3), '-', 'Color', c5, 'LineWidth', 2.5, ...
                'DisplayName', 'Optimized prediction');
        end
        if row == 1, title(colTitle{col}, 'Interpreter', 'tex'); end
        set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
            'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
            'XLim', [-120 20], 'XTick', [-120 -90 -60 -30 0], 'YLim', [-35 0], ...
            'TickLength', [0.025 0.05], 'GridLineStyle', 'none');
    end
end

% shared axes labels
ylabel(tT, '\bf Torque, N\cdotm', 'Interpreter', 'tex');
xlabel(tT, '\bf \theta_{k} , \circ', 'Interpreter', 'tex');

% legends in tiles 2 and 4 (right column; SouthWest, 8 pt — master convention)
lg2 = legend(axTiles(1,2), 'Location', 'southwest'); lg2.FontSize = 8;
lg4 = legend(axTiles(2,2), 'Location', 'southwest'); lg4.FontSize = 8;

% panel letters (exemplar positions)
tileLabels = {'(A)', '(B)', '(C)', '(D)'};
xAnn = [0.001, 0.48, 0.001, 0.48];
yAnn = [0.935, 0.935, 0.45, 0.45];
for j = 1:4
    annotation(figT, 'textbox', [xAnn(j) yAnn(j) 0.05 0.05], 'String', ['\bf ' tileLabels{j}], ...
        'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'HorizontalAlignment', 'left');
end

% row labels (PlotPrompt: y-label titles on tiles 1 and 3)
annotation(figT, 'textbox', [0.128 0.815 0.012 0.09], 'String', '\bf Pre-optimization', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation(figT, 'textbox', [0.128 0.33 0.012 0.09], 'String', '\bf Optimized', ...
    'FontSize', 12, 'FontName', 'Arial', 'EdgeColor', 'none', 'Rotation', 90, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% export EPS deliverable (same filename) + scratch PNG for QC
exportgraphics(figT, fullfile(outdir, 'FlxPin_group.eps'), 'ContentType', 'vector');
exportgraphics(figT, fullfile(qcdir, 'QC_FlxPin_group.png'), 'Resolution', 150);
fprintf('WROTE %s\n', fullfile(outdir, 'FlxPin_group.eps'));
close(figT);
disp('DONE_FlxPin_group');
