function export_sns_simscape_figures()
% Refresh the SNS-Simscape dissertation figures from the CURRENT library.
%
% Replaces the 2026-09-09-era exports (7-block pre-redesign library):
%   1. sns_library_overview  — composed grid of the 12 core-library blocks
%      using the REAL masked icons rendered at high dpi (names hidden on
%      canvas; block names and one-line functions are typeset text at
%      10.5-11 pt Arial at final print scale), grouped core vs additions.
%   2. simulink_neuron_detail — styled native internals of NonSpikingNeuron
%      (auto-arranged in memory so the strip prints with >=10 pt fonts).
%   3. components only: simulink_knee_reflex_arranged (styled) and
%      sns_winner_slice_simulink (the R8 winner slice).
%
% All styling is applied in memory; no source .slx is ever saved. Output
% text uses no dates and Oxford commas. Fonts: Arial, 10 pt floor.

repoRoot = fileparts(mfilename('fullpath'));
while ~strcmp(getLastFolder(repoRoot), 'Bipedal_Robot')
    parent = fileparts(repoRoot);
    assert(~strcmp(parent, repoRoot), 'Repository root not found.');
    repoRoot = parent;
end
srcDir = fullfile(repoRoot, 'Code', 'Matlab', 'SNS_Simscape');
figDir = fullfile(repoRoot, 'Documentation', 'Reports and Papers', ...
    'Dissertation', 'Figures', '60-sns');
compDir = fullfile(figDir, 'components');
tmpDir = fullfile(repoRoot, 'tmp', 'sns_fig_export_20261021');
if ~exist(tmpDir, 'dir'), mkdir(tmpDir); end
addpath(srcDir);
addpath(fullfile(srcDir, 'demos'));
fprintf('MATLAB %s\n', version);

%% ---------------------------------------------------------------- overview
blocks = {  % name, canvas w px, canvas h px, display name, summary
 'NonSpikingNeuron',  260, 260, 'NonSpiking\nNeuron',       'graded membrane\nvoltage, S(V) drive'
 'NonSpikingSynapse', 300, 150, 'NonSpiking\nSynapse',      'excitatory or\ninhibitory, Esyn sign'
 'SpikingLIFNeuron',  260, 260, 'SpikingLIF\nNeuron',       'integrate-and-fire\nwith threshold reset'
 'IaMuscleSpindle',   260, 260, 'IaMuscle\nSpindle',        'length and stretch\nvelocity encoder'
 'IbGolgiTendon',     260, 260, 'IbGolgi\nTendon',          'Golgi tendon organ\nforce encoder'
 'MuscleActivation',  260, 220, 'Muscle\nActivation',       'first-order\nactivation filter'
 'SynSum',            240, 170, 'SynSum',                   'sums up to eight\nsynapse lines'
 'BPA_10mm',          320, 116, 'BPA_10mm',                 '10 mm BPA force,\nFesto fit surface'
 'BPA_20mm',          320, 116, 'BPA_20mm',                 '20 mm BPA force,\nFesto fit surface'
 'BPA_40mm',          320, 116, 'BPA_40mm',                 '40 mm BPA force,\nFesto fit surface'
 'BioMuscle',         320, 123, 'BioMuscle',                'Thelen-style\nbiological muscle'
 'BPAForce',          280, 102, 'BPAForce',                 'generic BPA force\nplaceholder'};
for k = 1:size(blocks,1)
    blocks{k, 4} = strrep(blocks{k, 4}, '\n', newline);  % display name
    blocks{k, 5} = strrep(blocks{k, 5}, '\n', newline);  % summary
end

load_system(fullfile(srcDir, 'SNS_Library.slx'));
cleanupLib = onCleanup(@() closeClean('SNS_Library'));
wasLocked = get_param('SNS_Library', 'Lock');
set_param('SNS_Library', 'Lock', 'off');
set_param('SNS_Library', 'ScreenColor', 'white');

% snapshot original positions, then lay out the icon grid (names off);
% muscle-family icons render WITHOUT their in-Simulink disp labels here —
% the composed figure typesets the block names itself.
origPos = cell(size(blocks,1), 1);
pitchX = 440; pitchY = 360;
dispLabels = {'BPAForce', 'BPAForce'; 'BPA_10mm', 'BPA 10'; ...
              'BPA_20mm', 'BPA 20'; 'BPA_40mm', 'BPA 40'; ...
              'BioMuscle', 'BIO'};
for k = 1:size(blocks,1)
    full = ['SNS_Library/' blocks{k, 1}];
    origPos{k} = get_param(full, 'Position');
    x = 60 + pitchX * mod(k - 1, 3);
    y = 60 + pitchY * floor((k - 1) / 3);
    set_param(full, 'Position', ...
        [x y x + blocks{k, 2} y + blocks{k, 3}], 'ShowName', 'off');
end
for k = 1:size(dispLabels, 1)
    m = Simulink.Mask.get(['SNS_Library/' dispLabels{k, 1}]);
    m.Display = sns_icon_muscle(dispLabels{k, 2}, true);
end
open_system('SNS_Library');
drawnow; pause(1.0);   % let mask icons paint
canvasPng = fullfile(tmpDir, 'library_canvas.png');
print('-sSNS_Library', '-dpng', '-r400', canvasPng);

img = imread(canvasPng);
dpx = 400 / 72;   % canvas px -> device px
crops = cell(size(blocks,1), 1);
for k = 1:size(blocks,1)
    x = 60 + pitchX * mod(k - 1, 3);
    y = 60 + pitchY * floor((k - 1) / 3);
    w = blocks{k, 2}; h = blocks{k, 3};
    m = 28;
    c1 = max(1, floor((x - m) * dpx) + 1);
    c2 = min(size(img, 2), ceil((x + w + m) * dpx));
    r1 = max(1, floor((y - m) * dpx) + 1);
    r2 = min(size(img, 1), ceil((y + h + m) * dpx));
    crops{k} = img(r1:r2, c1:c2, :);
end

% restore the library display state exactly, then close without saving
for k = 1:size(blocks,1)
    set_param(['SNS_Library/' blocks{k, 1}], 'Position', origPos{k}, ...
              'ShowName', 'on');
end
for k = 1:size(dispLabels, 1)
    m = Simulink.Mask.get(['SNS_Library/' dispLabels{k, 1}]);
    m.Display = sns_icon_muscle(dispLabels{k, 2}, false);
end
set_param('SNS_Library', 'Lock', wasLocked);
clear cleanupLib;
close_system('SNS_Library', 0);

% ---- compose at dissertation print scale ----------------------------------
% tex: width=0.96\textwidth (6.5in -> 15.6 cm), height cap 0.42\textheight
% (~9.6 cm). A 15.3 x 9.9 cm figure prints at ~1.02 scale, so 10-11 pt
% fonts hold their face value at final size.
figW = 15.3; figH = 9.9;
f = figure('Units', 'centimeters', 'Position', [2 2 figW figH], ...
           'Color', 'w', 'Visible', 'off', 'PaperPositionMode', 'auto');
ax = axes('Parent', f, 'Units', 'normalized', 'Position', [0 0 1 1]);
axis(ax, 'off'); hold(ax, 'on');
set(ax, 'XLim', [0 figW], 'YLim', [0 figH], 'YDir', 'normal');

mL = 0.40; colW = (figW - 2 * mL) / 3;
headH = 0.50; rowH = (figH - 0.45 - 2 * headH) / 4;
groups = {[1 6], 'Core blocks'; [7 12], 'Full-library additions'};
for g = 1:size(groups, 1)
    rRow = floor((groups{g, 1}(1) - 1) / 3) + 1;   % 1-based grid row
    yh = figH - 0.30 - headH / 2 - (rRow - 1) * (rowH + headH);
    text(ax, mL, yh, groups{g, 2}, 'FontName', 'Arial', 'FontSize', 11, ...
         'FontWeight', 'bold', 'Color', [0.30 0.30 0.30], ...
         'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
end
for k = 1:size(blocks,1)
    col = mod(k - 1, 3); row = floor((k - 1) / 3);
    cx = mL + col * colW;
    cyTop = figH - 0.30 - headH - row * (rowH + headH);
    % icon: crop is pure canvas pixels; native cm = (px/dpx) * 2.54/72
    cw = size(crops{k}, 2); ch = size(crops{k}, 1);
    natW = (cw / dpx) * 2.54 / 72;   % cm
    natH = (ch / dpx) * 2.54 / 72;
    boxW = 1.60; boxH = rowH - 0.20;
    s = min(boxW / natW, boxH / natH);
    iw = natW * s; ih = natH * s;
    ix = cx + 0.10; iyTop = cyTop - 0.06;
    image(ax, 'CData', crops{k}, 'XData', [ix ix + iw], ...
          'YData', [iyTop iyTop - ih]);   % CData row 1 lands at cell top
    tx = cx + boxW + 0.25;
    ty = cyTop + 0.04;
    text(ax, tx, ty, blocks{k, 4}, 'FontName', 'Arial', ...
         'FontSize', 10.5, 'FontWeight', 'bold', 'Color', [0.10 0.10 0.10], ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'Interpreter', 'none');
    nl = numel(strsplit(blocks{k, 4}, newline));
    text(ax, tx, ty - 0.50 * nl - 0.14, blocks{k, 5}, 'FontName', 'Arial', ...
         'FontSize', 10, 'Color', [0.25 0.25 0.25], ...
         'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
         'Interpreter', 'none');
end
overviewStem = fullfile(figDir, 'sns_library_overview');
exportgraphics(f, [overviewStem '.pdf'], 'ContentType', 'vector');
print(f, '-dpng', '-r300', [overviewStem '.png']);
close(f);
writeOverviewAssets();
fprintf('overview written\n');

%% ------------------------------------------------------- neuron internals
load_system(fullfile(srcDir, 'SNS_Library.slx'));
cleanupLib2 = onCleanup(@() closeClean('SNS_Library'));
set_param('SNS_Library', 'Lock', 'off');
set_param('SNS_Library', 'ScreenColor', 'white');
sub = 'SNS_Library/NonSpikingNeuron';
open_system(sub, 'force');
drawnow; pause(0.8);
blks = find_system(sub, 'SearchDepth', 1, 'Type', 'Block');
origFont = nan(numel(blks), 1);
for b = 1:numel(blks)
    try
        origFont(b) = str2double(get_param(blks{b}, 'FontSize'));
        set_param(blks{b}, 'FontSize', '18');
    catch
    end
end
Simulink.BlockDiagram.arrangeSystem(sub);
drawnow; pause(0.5);
detStem = fullfile(figDir, 'simulink_neuron_detail');
print(['-s' sub], '-dpdf', '-bestfit', fullfile(tmpDir, 'neuron_raw.pdf'));
print(['-s' sub], '-dpng', '-r300', [detStem '.png']);
tightenPdf(fullfile(tmpDir, 'neuron_raw.pdf'), [detStem '.pdf']);
for b = 1:numel(blks)
    if ~isnan(origFont(b))
        set_param(blks{b}, 'FontSize', num2str(origFont(b)));
    end
end
clear cleanupLib2;
close_system('SNS_Library', 0);
fprintf('neuron detail written\n');

%% ------------------------------------------- components: reflex + winner
% styled in-memory exports via the library's own machinery; muscle-family
% instances are resized to the side-view aspect first (in memory only) so
% the wide icons do not squash.
exportWithMuscleAspect('KneeReflexDemo', true);
copyfile(fullfile(srcDir, 'figures', 'KneeReflexDemo_simulink.pdf'), ...
         fullfile(compDir, 'simulink_knee_reflex_arranged.pdf'));
copyfile(fullfile(srcDir, 'figures', 'KneeReflexDemo_simulink.png'), ...
         fullfile(compDir, 'simulink_knee_reflex_arranged.png'));

exportWithMuscleAspect('SNS_WinnerSlice_R8', false);
movefile(fullfile(compDir, 'SNS_WinnerSlice_R8_simulink.pdf'), ...
         fullfile(compDir, 'sns_winner_slice_simulink.pdf'));
movefile(fullfile(compDir, 'SNS_WinnerSlice_R8_simulink.png'), ...
         fullfile(compDir, 'sns_winner_slice_simulink.png'));
fprintf('components written\n');

%% ---------------------------------------------------------------- helpers
    function writeOverviewAssets()
        for k = 1:size(blocks,1)
            imwrite(crops{k}, fullfile(compDir, ['icon_' blocks{k, 1} '.png']));
        end
        alt = [ ...
'Alt text for sns_library_overview.pdf' newline newline ...
'Twelve blocks of the SNS-Simscape Simulink library arranged in a grid of ' ...
'three columns and four rows, each block shown as its masked icon with the ' ...
'block name and a one-line function summary to the right. The row group ' ...
'Core blocks contains the non-spiking neuron (open circle with a graded ' ...
'potential waveform), the non-spiking synapse (axon bar terminating in a ' ...
'white triangle for excitation), the spiking leaky integrate-and-fire ' ...
'neuron (circle with a spike train), the Ia muscle spindle encoder ' ...
'(spindle-shaped capsule), the Ib Golgi tendon encoder (capsule with ' ...
'braided strands), and the muscle activation filter (pentagon). The row ' ...
'group Full-library additions contains the synapse-sum junction, three ' ...
'diameter-specific BPA force elements drawn as side views of slightly ' ...
'inflated McKibben muscles with braided sleeves and clamp fittings for ' ...
'10, 20, and 40 mm, the Thelen-style biological muscle element drawn as ' ...
'a spindle with pennation-angle fibers running into tapering tendons, ' ...
'and the generic BPA force placeholder. The figure documents the library ' ...
'inventory and the icon language, not quantitative performance.' newline];
        fid = fopen(fullfile(compDir, 'sns_library_overview_alt.txt'), 'w');
        fwrite(fid, alt); fclose(fid);
    end
end

function exportWithMuscleAspect(mdl, arrange)
% Load mdl, resize muscle-family instances to the side-view icon aspect
% (display only, restored after), then run the styled export.
compDir = fullfile(fileparts(mfilename('fullpath')), '..', ...
    'Figures', '60-sns', 'components');
load_system(mdl);
blks = find_system(mdl, 'LookUnderMasks', 'all', 'Type', 'Block');
orig = {};
for b = 1:numel(blks)
    try
        mt = get_param(blks{b}, 'MaskType');
    catch
        mt = '';
    end
    asp = 0;
    if contains(mt, 'BPA'), asp = 2.75; end
    if contains(mt, 'Biological Muscle'), asp = 2.6; end
    if asp > 0
        p = get_param(blks{b}, 'Position');
        orig(end+1, :) = {blks{b}, p}; %#ok<AGROW>
        set_param(blks{b}, 'Position', ...
            [p(1) p(2) p(3) p(2) + round((p(3) - p(1)) / asp)]);
    end
end
sns_export_diagram(compDir, {mdl}, arrange);
for b = 1:size(orig, 1)
    set_param(orig{b, 1}, 'Position', orig{b, 2});
end
close_system(mdl, 0);
end

function tightenPdf(src, dst)
pc = fullfile(getenv('LOCALAPPDATA'), 'Programs', 'MiKTeX', 'miktex', ...
    'bin', 'x64', 'pdfcrop.exe');
if ~isfile(pc)
    copyfile(src, dst);
    fprintf('pdfcrop not found; %s copied untightened\n', dst);
    return;
end
system(['"' pc '" --margin 4 "' src '" "' dst '"']);
end

function closeClean(mdl)
if bdIsLoaded(mdl), close_system(mdl, 0); end
end

function leaf = getLastFolder(path)
[~, leaf] = fileparts(path);
end
