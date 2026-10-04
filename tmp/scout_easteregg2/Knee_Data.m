

%% Reverse-Pulley Knee Flexor Data (MonoPam_pulley spin-off)
% Knee_Data.m -- journal-style results caller for the reverse-pulley
% (block-and-tackle) knee-flexor optimization (Ben, 2026-09-24). In the
% image of Knee_Flexor_Data_20mm.m: self-locates the repo root, loads the
% newest (or a named) Bifemsh_20mm_Result_pulley_<stamp>.mat, rebuilds the
% pulley prediction from xBest + routeCtx + Xi3 + ctx with MonoPam_pulley
% objects, and plots torque vs human target, muscle length, strain, moment
% arm, torque margin, plus the pulley quantities (insertion force,
% reaction force, and the transmission config -- nPulleyBPA,
% tackleLineParts, routing mode -- as an annotated panel and in the
% figure title) so the transmission trade is visible. Plotting convention
% as in minimizeFlxPin10mm's color-scheme sections: Colors.m palette,
% Arial sizes 10/12/8, the Ben 3-line allBPA override block, auto
% tileLabels (A),(B),... and the legend rule (even panel count -> legend
% in tile (1,2); odd -> first empty tile).

%% Freshen up the workspace
clc
clear
close all

set(groot, ...
    'defaultAxesFontName','Arial', ...
    'defaultAxesFontSize',10, ...
    'defaultAxesFontWeight','bold', ...
    'defaultAxesLabelFontSizeMultiplier',1, ...
    'defaultAxesTitleFontSizeMultiplier',1.2, ...
    'defaultAxesLineWidth',2, ...
    'defaultAxesBox','off', ...
    'defaultAxesXMinorTick','on', ...
    'defaultAxesYMinorTick','on', ...
    'defaultAxesTickLength',[0.025 0.05], ...
    'defaultAxesXGrid','off', ...
    'defaultAxesYGrid','off', ...
    'defaultLineLineWidth',2, ...
    'defaultLegendFontName','Arial', ...
    'defaultLegendFontSize',8, ...
    'defaultLegendFontWeight','bold', ...
    'defaultLegendBox','off')

%% Add paths to the muscle and pam calculators
% Repo root and path setup (same block as Knee_Flexor_Data_20mm.m):
% resolves the muscle classes, Mesh_Optimization, Colors.m, and the
% OpenSim_Bifem txt files in Testing_Data\2022_02_Festo, regardless of cwd.
scriptDir = fileparts(mfilename('fullpath'));
root = scriptDir;
for k = 1:8
    [parent, name] = fileparts(root);
    if strcmpi(name, 'Bipedal_Robot')
        break
    end
    if strcmp(parent, root)
        error('Could not locate the Bipedal_Robot repo root from %s', scriptDir)
    end
    root = parent;
end
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
% Mesh_Optimization must win any shadowing contest against data subfolders.
addpath(fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization'));
% Append (do not prepend) so Code\Matlab keeps winning name collisions.
addpath(fullfile(root, 'Testing_Data', '2022_02_Festo'), '-end');

%% Load the newest (or a named) pulley result mat
% Set resultFile to a full path to pin one run; empty = newest match.
resultFile = '';

resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
if isempty(resultFile)
    dResults = dir(fullfile(resDir, 'Bifemsh_20mm_Result_pulley_*.mat'));
    if isempty(dResults)
        error('Knee_Data:NoResultMat', ...
            ['No Bifemsh_20mm_Result_pulley_*.mat in %s. Run ' ...
             'Opt_run_pulley.m (a real run, not SMOKE) first.'], resDir)
    end
    [~, idxNewest] = max([dResults.datenum]);
    resultFile = fullfile(dResults(idxNewest).folder, dResults(idxNewest).name);
end
fprintf('Loading %s\n', resultFile)
resultData = load(resultFile);

requiredResultFields = {'xBest', 'routeCtx', 'Xi3', 'ctx'};
for iField = 1:numel(requiredResultFields)
    if ~isfield(resultData, requiredResultFields{iField})
        error('Knee_Data:MissingResultField', ...
            ['%s does not contain %s. Rerun the updated ' ...
             'Opt_run_pulley.m.'], resultFile, requiredResultFields{iField})
    end
end

xBest = resultData.xBest;
routeCtx = resultData.routeCtx;
Xi3 = resultData.Xi3;
ctx = resultData.ctx;
if isfield(resultData, 'fBest')
    fBest = resultData.fBest;
    fprintf('fBest                        = %.9g\n', fBest)
end

% Backfill pulley knobs for mats from before the spin-off knobs existed;
% keep values in sync with Opt_run_pulley.
if ~isfield(ctx, 'pulleyBPACount'), ctx.pulleyBPACount = 2; end
if ~isfield(ctx, 'tackleLineParts'), ctx.tackleLineParts = 1; end
if ~isfield(ctx, 'pulleyRoutingMode'), ctx.pulleyRoutingMode = 'moving_via'; end
if ~isfield(ctx, 'pulleyExitIndex'), ctx.pulleyExitIndex = 1; end
if ~isfield(ctx, 'bowdenBossDia'), ctx.bowdenBossDia = 0.007; end
if ~isfield(ctx, 'bowdenRunClearance'), ctx.bowdenRunClearance = 0.005; end
if ~isfield(ctx, 'optimizePulleyGain'), ctx.optimizePulleyGain = false; end
if ~isfield(routeCtx, 'pulleyBPACount'), routeCtx.pulleyBPACount = ctx.pulleyBPACount; end
if ~isfield(routeCtx, 'tackleLineParts'), routeCtx.tackleLineParts = ctx.tackleLineParts; end
if ~isfield(routeCtx, 'pulleyRoutingMode'), routeCtx.pulleyRoutingMode = ctx.pulleyRoutingMode; end
if ~isfield(routeCtx, 'pulleyExitIndex'), routeCtx.pulleyExitIndex = ctx.pulleyExitIndex; end
if ~isfield(routeCtx, 'bowdenBossDia'), routeCtx.bowdenBossDia = ctx.bowdenBossDia; end
if ~isfield(routeCtx, 'bowdenRunClearance'), routeCtx.bowdenRunClearance = ctx.bowdenRunClearance; end

%% Rebuild the pulley prediction from xBest + routeCtx + Xi3 + ctx
p1 = xBest(1:3); %Origin, femur frame
p2 = xBest(4:6); %End/insertion, t1 frame
rest   = xBest(7);
tendon = xBest(8);
KMAX = ctx.KMAX;
kmax = (1 - KMAX)*rest;
fitting = ctx.fitting;

% The gain: the ninth variable when it was optimized, else the run's
% tackle line parts (tackleLineParts = travel gain G).
if ctx.optimizePulleyGain && numel(xBest) > 8
    pulleyGain = xBest(9);
    fprintf('Pulley gain OPTIMIZED: G = %.6f\n', pulleyGain)
else
    pulleyGain = routeCtx.tackleLineParts;
    fprintf('Pulley gain fixed: G = %.6f\n', pulleyGain)
end
fprintf('Xi3 (saved by Opt_run_pulley, unused by MonoPam_pulley) = %.6g\n', Xi3)

% Run-level transmission config for the MonoPam_pulley constructor.
% Shared-tackle aggregation (same rule as the driver's predictor): with
% the two-route pipeline the physical parallel BPAs ARE the two routes, so
% each route object carries ONE BPA line and the route SUM is the rig's
% total pull F_t = (F1 + F2)/G. The class's nPulleyBPA > 1 path is used
% only for single-corridor designs.
if ctx.BPAcount == 1
    perRouteN = routeCtx.pulleyBPACount;
else
    perRouteN = 1;
end
pulleyConfig = struct( ...
    'nPulleyBPA', perRouteN, ...
    'tackleLineParts', routeCtx.tackleLineParts, ...
    'gain', pulleyGain, ...
    'routingMode', routeCtx.pulleyRoutingMode, ...
    'pulleyExitIndex', routeCtx.pulleyExitIndex, ...
    'bowdenBossDia', routeCtx.bowdenBossDia, ...
    'bowdenRunClearance', routeCtx.bowdenRunClearance);
configNote = sprintf('n_{BPA} = %d, G = %g, %s mode', ...
    routeCtx.pulleyBPACount, pulleyGain, pulleyConfig.routingMode);

% Use the saved kinematic/routing context with the CURRENT geometry
% definitions (same hand-tune policy as Knee_Flexor_Data_20mm).
routeCtx.geo = buildGeoExclusion();

[Location, ~, routeInfo] = ...
    buildKneeFlexorRoute20mm(p1, p2, tendon, routeCtx);

% BPA 2 (Ben, 2026-09-21 asymmetric routing): pEnd{2} keeps the mirrored
% distal attachment and p1{2} shares BPA 1's side of the knee.
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[LocationB, ~, routeInfoB] = ...
    buildKneeFlexorRoute20mm(p1B, p2B, tendon, routeCtx);

% Optimized design: BOTH routes through MonoPam_pulley, as in Opt_run_pulley.
Pam_opt_1 = MonoPam_pulley(ctx.Name, Location, ctx.CrossPoint, ctx.Dia, ...
    ctx.T_Pam, rest, kmax, tendon, fitting, ctx.targetPressure, ...
    ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.wraps, pulleyConfig);
Pam_opt_2 = MonoPam_pulley(ctx.Name, LocationB, ctx.CrossPoint, ctx.Dia, ...
    ctx.T_Pam, rest, kmax, tendon, fitting, ctx.targetPressure, ...
    ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.wraps, pulleyConfig);

% Original two-point route, no wrap, and NO tackle (the original rig had a
% straight tendon): same MonoPam_pulley class at G = 1. Its rest length is
% sized INTO the strain window for this short route (ctx.originalRest =
% 0.415 m sits out of the model's range on the straight two-point route;
% MonoPam_pulley flags that infeasible instead of silently zeroing force,
% so the comparator is length-matched to stay comparable).
Location0 = zeros(2, 3, ctx.N);
routeLen0 = zeros(ctx.N, 1);
for i = 1:ctx.N
    p20ICR = RowVecTrans(ctx.T_ICR_t1(:,:,i), ctx.originalPEnd);
    Location0(:,:,i) = [ctx.originalP1; p20ICR];
    routeLen0(i) = norm(ctx.originalP1 - RowVecTrans(ctx.T_Pam(:,:,i), p20ICR));
end
rest0 = max(routeLen0) - ctx.Xi0 - ctx.originalTendon - 2*fitting;
kmax0 = (1 - KMAX)*rest0;
cfgStraight = struct('nPulleyBPA', 1, 'tackleLineParts', 1, ...
    'pulleyExitIndex', 1);
Pam_orig = MonoPam_pulley(ctx.Name, Location0, ctx.CrossPoint, ctx.Dia, ...
    ctx.T_Pam, rest0, kmax0, ctx.originalTendon, fitting, ...
    ctx.targetPressure, ctx.Xi0, ctx.Xi1, ctx.Xi2, ctx.originalWraps, ...
    cfgStraight);

% Route-pair sums (the same aggregation the driver's predictor uses).
TorqueZ = Pam_opt_1.Torque_p(:,3) + Pam_opt_2.Torque_p(:,3);
TorqueInsZ = Pam_opt_1.Torque_ins(:,3) + Pam_opt_2.Torque_ins(:,3);
TorqueX = Pam_opt_1.Torque_p(:,1) + Pam_opt_2.Torque_p(:,1);
TorqueY = Pam_opt_1.Torque_p(:,2) + Pam_opt_2.Torque_p(:,2);
activeLength = rest .* (1 - Pam_opt_1.strain_p(:));
strain_p = Pam_opt_1.strain_p(:);
Contraction = Pam_opt_1.Contraction(:);
momentArm = hypot(Pam_opt_1.mA_p(:,1), Pam_opt_1.mA_p(:,2));
momentArm0 = hypot(Pam_orig.mA_p(:,1), Pam_orig.mA_p(:,2));
FinsMag = vecnorm(Pam_opt_1.F_ins, 2, 2) + vecnorm(Pam_opt_2.F_ins, 2, 2);
ReactionMag = Pam_opt_1.ReactionFmag + Pam_opt_2.ReactionFmag;
activeLength0 = rest0 .* (1 - Pam_orig.strain_p(:));
TorqueZ0 = Pam_orig.Torque_p(:,3);
infeasibleCount = nnz(Pam_opt_1.PulleyInfeasible | Pam_opt_2.PulleyInfeasible);
slackCount = nnz(Pam_opt_1.PulleySlack | Pam_opt_2.PulleySlack);

%% Human target
H = readmatrix('OpenSim_Bifem_Results.txt', ...
    'FileType', 'text', ...
    'NumHeaderLines', 7);
humanAngle = H(:,2);
TorqueHz = H(:,4);
humanTorqueAbs = abs(TorqueHz);

%% Plotting convention (Ben, 2026-09-12)
% Use the established accessible project palette.  Colors.m supplies c for
% line colors and d for RGB scatter-marker colors; these figures use c.
run('Colors.m')
originalColor = [0.4 0.4 0.4];
optimizedColor = c{5};
humanColor = '#000000';
limitColor = c{7};

fontName = 'Arial';
axesFontSize = 10;
titleFontSize = 12;
legendFontSize = 8;
axesLineWidth = 2;
originalLineWidth = 2;
optimizedLineWidth = 2.5;
humanLineWidth = 4;
tickLength = [0.025 0.05];
xLimits = [min(ctx.phiD), max(ctx.phiD)];

% Ben's 3-line allBPA override block (kept in the established form). Here
% allBPA selects the plotted BPA ROUTES (1 = optimized route, 2 = derived
% BPA 2); the PANEL count plays numBPA's role from minimizeFlxPin10mm for
% tileLabels and the legend rule.
% allBPA = allBPA;
allBPA = [1, 2];
% allBPA = [1, 2, 3, 4, 5];
numBPA = numel(allBPA);

nPanels = 7;  % torque, length, strain, moment arm, margin, forces, config
tileLabels = arrayfun(@(k) sprintf('(%c)', 'A' + k - 1), 1:nPanels, ...
    'UniformOutput', false);

%% Journal figure: seven panels
figure('Name', 'Reverse-Pulley Flexor Data', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2 2 18 21])
tLayout = tiledlayout(ceil(nPanels/2), 2, ...
    'TileSpacing', 'loose', 'Padding', 'loose');
title(tLayout, ['Reverse-Pulley Flexor (' configNote ')' ], ...
    'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')

%% (A) Torque vs human target
axA = nexttile(1);
hold(axA, 'on')
plot(axA, ctx.phiD, TorqueZ0, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA (no tackle)')
plot(axA, ctx.phiD, TorqueZ, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Route-bookkeeping torque')
plot(axA, ctx.phiD, TorqueInsZ, '-.', ...
    'Color', limitColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Insertion torque')
plot(axA, humanAngle, -humanTorqueAbs, ':', ...
    'Color', humanColor, 'LineWidth', humanLineWidth, ...
    'DisplayName', 'Human target')
title(axA, 'Flexor Torque', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
stylePanel(axA, fontName, axesFontSize, axesLineWidth, tickLength, ...
    xLimits, 'Torque, N\cdotm')

%% (B) Muscle length
axB = nexttile(2);
hold(axB, 'on')
plot(axB, ctx.phiD, activeLength0, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(axB, ctx.phiD, activeLength, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
title(axB, 'Muscle Length, L_m', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', titleFontSize, 'FontWeight', 'bold')
stylePanel(axB, fontName, axesFontSize, axesLineWidth, tickLength, ...
    xLimits, 'Muscle Length, m')

%% (C) Strain definitions
axC = nexttile(3);
hold(axC, 'on')
plot(axC, ctx.phiD, strain_p, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'strain\_p')
plot(axC, ctx.phiD, Contraction, '-.', ...
    'Color', c{6}, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Contraction')
minLine = yline(axC, 0, ':', 'Minimum strain', ...
    'Color', humanColor, 'LineWidth', originalLineWidth, ...
    'HandleVisibility', 'off');
maxLine = yline(axC, KMAX, ':', 'KMAX', ...
    'Color', limitColor, 'LineWidth', originalLineWidth, ...
    'HandleVisibility', 'off');
set([minLine, maxLine], 'FontName', fontName, ...
    'FontSize', legendFontSize, 'FontWeight', 'bold')
title(axC, 'Strain Definitions', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
stylePanel(axC, fontName, axesFontSize, axesLineWidth, tickLength, ...
    xLimits, 'Strain')

%% (D) Moment arm
axD = nexttile(4);
hold(axD, 'on')
plot(axD, ctx.phiD, momentArm0, '--', ...
    'Color', originalColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Original BPA')
plot(axD, ctx.phiD, momentArm, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Optimized BPA')
title(axD, 'Moment Arm', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
stylePanel(axD, fontName, axesFontSize, axesLineWidth, tickLength, ...
    xLimits, 'Moment Arm, m')

%% (E) Torque margin fraction (BPA side and after the tackle)
humanAbsAtRobotAngles = interp1(humanAngle, humanTorqueAbs, ctx.phiD, ...
    'pchip', 'extrap');
validHumanTorque = humanAbsAtRobotAngles > 100*eps;
torqueMarginFraction = nan(size(TorqueZ));
torqueMarginFraction(validHumanTorque) = ...
    abs(TorqueZ(validHumanTorque)) ./ ...
    humanAbsAtRobotAngles(validHumanTorque) - 1;
insertionMarginFraction = nan(size(TorqueInsZ));
insertionMarginFraction(validHumanTorque) = ...
    abs(TorqueInsZ(validHumanTorque)) ./ ...
    humanAbsAtRobotAngles(validHumanTorque) - 1;

axE = nexttile(5);
hold(axE, 'on')
plot(axE, ctx.phiD, 100*torqueMarginFraction, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', 'Route-bookkeeping margin')
plot(axE, ctx.phiD, 100*insertionMarginFraction, '--', ...
    'Color', limitColor, 'LineWidth', originalLineWidth, ...
    'DisplayName', 'Insertion margin')
zeroLine = yline(axE, 0, ':', 'Human target', ...
    'Color', humanColor, 'LineWidth', originalLineWidth, ...
    'HandleVisibility', 'off');
zeroLine.FontName = fontName;
zeroLine.FontSize = legendFontSize;
zeroLine.FontWeight = 'bold';
title(axE, 'Torque Margin vs Human', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
stylePanel(axE, fontName, axesFontSize, axesLineWidth, tickLength, ...
    xLimits, 'Torque Margin, %')

%% (F) Pulley transmission trade: insertion force and mount reaction
axF = nexttile(6);
hold(axF, 'on')
plot(axF, ctx.phiD, FinsMag, '-', ...
    'Color', optimizedColor, 'LineWidth', optimizedLineWidth, ...
    'DisplayName', '|F\_ins|')
plot(axF, ctx.phiD, ReactionMag, '-.', ...
    'Color', c{6}, 'LineWidth', originalLineWidth, ...
    'DisplayName', '||ReactionF||')
title(axF, 'Insertion and Reaction Forces', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')
stylePanel(axF, fontName, axesFontSize, axesLineWidth, tickLength, ...
    xLimits, 'Force, N')

%% (G) Transmission config annotation (nPulleyBPA, G, routing mode)
axG = nexttile(7);
hold(axG, 'on')
axis(axG, 'off')
configLines = { ...
    sprintf('nPulleyBPA (parallel BPAs)  = %d', routeCtx.pulleyBPACount), ...
    sprintf(['(each of the %d route corridors carries 1 BPA line'], ...
        ctx.BPAcount), ...
    sprintf(' into the shared tackle; sum = rig total)'), ...
    sprintf('tackleLineParts = gain G    = %d', ...
        pulleyConfig.tackleLineParts), ...
    sprintf('routing mode                = %s', ...
        pulleyConfig.routingMode), ...
    sprintf('tackle exit row index       = %d', ...
        pulleyConfig.pulleyExitIndex), ...
    sprintf(['Bowden envelope: boss %.1f mm, run clearance %.1f mm'], ...
        1000*pulleyConfig.bowdenBossDia, ...
        1000*pulleyConfig.bowdenRunClearance), ...
    sprintf('infeasible / slack frames   = %d / %d of %d', ...
        infeasibleCount, slackCount, ctx.N)};
text(axG, 0.02, 0.90, configLines, 'Units', 'normalized', ...
    'FontName', fontName, 'FontSize', legendFontSize, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
    'Interpreter', 'none', 'Clipping', 'off')
title(axG, 'Transmission Config', 'FontName', fontName, ...
    'FontSize', titleFontSize, 'FontWeight', 'bold')

%% Auto tile labels (A)-(G), bold, top-left of each tile (Ben convention)
for k = 1:nPanels
    axK = tLayout.Children(end - k + 1);   % Children are creation-reversed
    text(axK, 0.002, 1.03, tileLabels{k}, 'Units', 'normalized', ...
        'Clipping', 'off', 'FontSize', 14, 'FontWeight', 'bold', ...
        'FontName', 'Arial', 'VerticalAlignment', 'bottom')
end

%% Legend rule (as in minimizeFlxPin10mm): even number of panels -> legend
% inside tile (1,2); odd -> the first empty tile of the grid.
if mod(nPanels, 2) == 0
    lg = legend(tLayout.Children(nPanels-1));  %tile (1,2) series
else
    lg = legend(tLayout.Children(end));        %tile 1 series, moved to the empty tile
    lg.Layout.Tile = 2*ceil(nPanels/2);
end
lg.Location = 'best';
lg.FontName = fontName;
lg.FontSize = legendFontSize;
lg.FontWeight = 'bold';
lg.Box = 'off';

ylabel(tLayout, 'Torque / Margin / Force', 'FontName', fontName, ...
    'FontSize', axesFontSize, 'FontWeight', 'bold');
xlabel(tLayout, '\theta_{k} , \circ', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold');

%% Console summary
fprintf('\n========== REVERSE-PULLEY FLEXOR DATA ==========\n')
fprintf('result mat                       = %s\n', resultFile)
fprintf('transmission                     = %s\n', configNote)
fprintf('pulley gain G                    = %.6f\n', pulleyGain)
fprintf('tackle exit row index            = %d\n', ...
    pulleyConfig.pulleyExitIndex)
fprintf('Bowden envelope                  = boss %.1f mm, run clearance %.1f mm\n', ...
    1000*pulleyConfig.bowdenBossDia, 1000*pulleyConfig.bowdenRunClearance)
fprintf('rest / kmax                      = %.6f / %.6f m\n', rest, kmax)
fprintf('tendon length                    = %.6f m\n', tendon)
fprintf('Xi0 / Xi1 / Xi2 / Xi3            = %.6g / %.6g / %.6g / %.6g\n', ...
    ctx.Xi0, ctx.Xi1, ctx.Xi2, Xi3)
fprintf('wrap release found               = %d\n', routeInfo.releaseFound)
if routeInfo.releaseFound
    fprintf('first inactive wrap angle        = %+.6f deg\n', ...
        routeInfo.releaseAngleD)
end
fprintf('BPA 2 wrap release found         = %d\n', routeInfoB.releaseFound)
fprintf('infeasible frames (both routes)  = %d of %d\n', ...
    infeasibleCount, ctx.N)
fprintf('slack frames (both routes)       = %d of %d\n', slackCount, ctx.N)
fprintf('min Contraction/KMAX             = %+.6f\n', min(Contraction)/KMAX)
fprintf('max Contraction/KMAX             = %+.6f\n', max(Contraction)/KMAX)
fprintf('|F_ins| range, routes summed     = %.3f to %.3f N\n', ...
    min(FinsMag), max(FinsMag))
fprintf('||ReactionF|| range, routes sum  = %.3f to %.3f N\n', ...
    min(ReactionMag), max(ReactionMag))
fprintf('min route-bookkeeping margin     = %+.6f (%+.2f%%)\n', ...
    min(torqueMarginFraction(validHumanTorque)), ...
    100*min(torqueMarginFraction(validHumanTorque)))
fprintf('min INSERTION margin fraction    = %+.6f (%+.2f%%)\n', ...
    min(insertionMarginFraction(validHumanTorque)), ...
    100*min(insertionMarginFraction(validHumanTorque)))
fprintf('================================================\n')

%% Local: one styled panel (Ben's axis conventions from Opt_run.m)
function stylePanel(ax, fontName, axesFontSize, axesLineWidth, ...
        tickLength, xLimits, yLabel)
set(ax, 'FontName', fontName, 'FontSize', axesFontSize, ...
    'FontWeight', 'bold', 'LineWidth', axesLineWidth, ...
    'Box', 'off', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', tickLength, 'XLim', xLimits, ...
    'GridLineStyle', 'none')
xlabel(ax, '\theta_k, °', 'Interpreter', 'tex', ...
    'FontName', fontName, 'FontSize', axesFontSize, 'FontWeight', 'bold')
ylabel(ax, yLabel, 'FontName', fontName, ...
    'FontSize', axesFontSize, 'FontWeight', 'bold')
grid(ax, 'off')
end
