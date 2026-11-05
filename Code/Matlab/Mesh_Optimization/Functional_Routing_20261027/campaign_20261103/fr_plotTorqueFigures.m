function fr_plotTorqueFigures()
%FR_PLOTTORQUEFIGURES Journal-standard signed torque figures (Ben's
% figure specification: Tol palette, Arial, 10 pt floor, CVD-safe
% encoding = color AND linestyle AND marker, alt-text file per figure,
% exportgraphics PNG). Reads output/torque_curves/<group>.csv written by
% fr_deliverables; writes output/figures_signed/.
%
% Run headless: matlab -batch "addpath('<camp>'); fr_plotTorqueFigures"

curveDir = fullfile('output', 'metrics');
outDir = fullfile('output', 'figures_signed');
if ~exist(outDir, 'dir'); mkdir(outDir); end
d0 = dir(fullfile(curveDir, '*_rigid.mat'));
assert(~isempty(d0), 'fr_plotTorqueFigures:noCurves', ...
    'no torque_curves CSVs - run fr_deliverables first');

% Tol palette (Code/Matlab/Colors.m lineage), CVD-safe triples:
% human = black solid; rigid = Tol blue #4477AA, dashed, circle markers;
% compliant = Tol red #EE6677, dash-dot, square markers.
colHuman = [0 0 0];
colRigid = [0.267 0.467 0.667];
colComp = [0.933 0.400 0.467];
set(groot, 'DefaultAxesFontName', 'Arial');
set(groot, 'DefaultTextFontName', 'Arial');
set(groot, 'DefaultAxesFontSize', 10);
set(groot, 'DefaultLegendFontSize', 8);

for i = 1:numel(d0)
    S = load(fullfile(curveDir, d0(i).name));
    [~, base] = fileparts(d0(i).name);
    gname = strrep(base, '_rigid', '');
    tgt = S.tauT(:);
    rig = S.tP(:);
    deg = S.degV(:);
    ckC = fullfile(curveDir, sprintf('%s_compliant.mat', gname));
    assert(exist(ckC, 'file') == 2, 'missing compliant mat for %s', gname);
    S2 = load(ckC);
    comp = S2.tP(:);
    n = numel(deg);

    fig = figure('Visible', 'off', 'Units', 'centimeters', ...
        'Position', [2 2 13.2 8.1]);
    ax = axes(fig);
    hold(ax, 'on');
    plot(ax, deg, tgt, '-', 'Color', colHuman, 'LineWidth', 1.6, ...
        'DisplayName', 'Human target (stock model)');
    plot(ax, deg, rig, '--o', 'Color', colRigid, 'LineWidth', 1.4, ...
        'MarkerSize', 5, 'MarkerFaceColor', colRigid, ...
        'DisplayName', 'Robot rigid');
    plot(ax, deg, comp, '-.s', 'Color', colComp, 'LineWidth', 1.4, ...
        'MarkerSize', 5, 'MarkerFaceColor', colComp, ...
        'DisplayName', 'Robot compliant (Xi_0-Xi_3)');
    xlabel(ax, 'Sweep coordinate (deg)');
    ylabel(ax, 'Torque (N*m)');
    title(ax, gname, 'FontSize', 10, 'FontWeight', 'normal');
    grid(ax, 'on');
    ax.GridAlpha = 0.3;
    legend(ax, 'Location', 'best', 'Box', 'off');
    png = fullfile(outDir, [gname '.png']);
    exportgraphics(fig, png, 'Resolution', 150);

    % signed meet stats for the alt text
    prj = sign(tgt);
    prj(abs(tgt) <= 1e-12) = 0;
    meetR = mean(rig .* prj >= abs(tgt) - 1e-9);
    meetC = mean(comp .* prj >= abs(tgt) - 1e-9);
    alt = sprintf(['Line chart for group %s: signed torque about the ' ...
        'primary dof across %d pose cells. Human target (black solid) ' ...
        'vs optimized robot route, rigid branch (blue dashed, circles; ' ...
        'signed meet %d/%d cells) and compliant branch Xi0-Xi3 (red ' ...
        'dash-dot, squares; %d/%d cells). Wrong-sign torque counts as ' ...
        'a miss.'], gname, n, round(meetR * n), n, ...
        round(meetC * n), n);
    fid = fopen(fullfile(outDir, [gname '_alttext.txt']), 'w');
    fwrite(fid, alt);
    fclose(fid);
    close(fig);
    fprintf('figure %d/%d: %s (meet rigid %d/%d, compliant %d/%d)\n', ...
        i, numel(d0), gname, round(meetR * n), n, round(meetC * n), n);
end
fprintf('figures written to %s\n', outDir);
end
