%% sns_draw_circuit.m — journal-style redraw of the KneeReflexDemo circuit
%
% Produces the reader-facing circuit figure in the diagram language of
% Szczecinski et al. 2017 (Fig. 2) / Rybak & Shevtsova / Animatlab:
%   open circle = neuron, Ia/Ib circles = afferents, ellipses = muscles,
%   open triangle = excitatory contact, solid black dot = inhibitory contact,
%   dashed = sensory feedback. Shape-only coding = colorblind safe.
% Matches the Simulink topology in sns_build_demo.m 1:1.
% Outputs: figures/KneeReflex_circuit.{png,pdf,svg}

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
outDir = fullfile(cdto, 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 11.5], 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, [0 100 0 62]);
axis(ax, 'off');

rAff = 3.2;   % afferent radius
rMN  = 4.2;   % motoneuron radius

% ---- motoneurons ----
snsfig('neuron', 46, 42, rMN, 'MN', 'extensor');
snsfig('neuron', 46, 19, rMN, 'MN', 'flexor');

% ---- afferents (left column) ----
snsfig('afferent', 13, 47, rAff, 'Ia', 'ext');
snsfig('afferent', 13, 38, rAff, 'Ib', 'ext');
snsfig('afferent', 13, 24, rAff, 'Ia', 'flex');
snsfig('afferent', 13, 13, rAff, 'Ib', 'flex');

% ---- synaptic contacts (E/I markers at postsynaptic side) ----
% [source x y, MN x y, approach angle at MN (deg), type]
connList = { ...
    13, 47, 46, 42, 165, 'exc';   % Ia ext  + MN ext   (stretch reflex)
    13, 38, 46, 42, 185, 'inh';   % Ib ext  - MN ext   (autoinhibition)
    13, 24, 46, 42, 205, 'inh';   % Ia flex - MN ext   (reciprocal inhibition)
    13, 47, 46, 19, 150, 'inh';   % Ia ext  - MN flex  (reciprocal inhibition)
    13, 24, 46, 19, 185, 'exc';   % Ia flex + MN flex
    13, 13, 46, 19, 215, 'inh'};  % Ib flex - MN flex
rMk = 6.0;  % marker radius from MN center
for k = 1:size(connList, 1)
    drawConn(connList{k,1}, connList{k,2}, connList{k,3}, connList{k,4}, ...
        rAff, rMk, connList{k,5}, connList{k,6});
end

% ---- descending drive (bias) ----
snsfig('arrow', 46, 49.5, 46, 46.6);
snsfig('arrow', 46, 26.3, 46, 23.6);
snsfig('label', 52.5, 51.5, 'descending drive', 7.5, 'gi');

% ---- muscles + activation arrows ----
snsfig('arrow', 50.4, 43.0, 60.6, 44.7);
snsfig('arrow', 50.4, 18.0, 60.6, 16.3);
snsfig('label', 55.6, 45.4, 'A', 7.5, 'gi');
snsfig('label', 55.6, 14.6, 'A', 7.5, 'gi');
snsfig('muscle', 66, 45, 10, 4.6, 'BPA ext');
snsfig('muscle', 66, 16, 10, 4.6, 'BPA flex');

% ---- knee joint plant ----
snsfig('box', 86, 31, 11, 14, sprintf('knee\njoint'));
snsfig('arrow', 71.3, 44.2, 80.4, 34.6);
snsfig('arrow', 71.3, 16.8, 80.4, 27.4);
snsfig('label', 77.3, 41.4, 'T_{ext}', 7.5, 'gi');
snsfig('label', 77.3, 19.6, 'T_{flex}', 7.5, 'gi');
snsfig('label', 86, 22.2, '\theta  (1-DOF)', 7.5, 'gi');

% ---- sensory feedback (dashed) ----
snsfig('edge', 86, 24, 86, 6.5, '--');
snsfig('edge', 86, 6.5, 13, 6.5, '--');
snsfig('edge', 13, 6.5, 13, 9.6, '--');
snsfig('arrow', 13, 6.5, 13, 9.6, '--');
snsfig('edge', 13, 6.5, 8, 6.5, '--');
snsfig('edge', 8, 6.5, 8, 47, '--');
snsfig('arrow', 8, 47, 9.6, 47, '--');
snsfig('arrow', 8, 24, 9.6, 24, '--');
snsfig('edge', 13, 6.5, 4, 6.5, '--');
snsfig('edge', 4, 6.5, 4, 38, '--');
snsfig('arrow', 4, 38, 9.6, 38, '--');
snsfig('label', 55, 7.9, 'sensory feedback  (\theta, \theta'', F)', 7.5, 'gi');

% ---- key (top right) ----
rectangle('Position', [64 50.5 34.5 10], 'EdgeColor', [0.6 0.6 0.6], ...
    'FaceColor', 'w', 'LineWidth', 0.6);
snsfig('exc', 67.3, 57.9, 0, 2.4);
snsfig('label', 79, 57.9, 'excitatory', 7.5, 'n');
snsfig('inh', 67.3, 55.2, 2.4);
snsfig('label', 79, 55.2, 'inhibitory', 7.5, 'n');
snsfig('edge', 65.3, 52.5, 69.3, 52.5, '--');
snsfig('label', 79, 52.5, 'sensory feedback', 7.5, 'n');

exportgraphics(fig, fullfile(outDir, 'KneeReflex_circuit.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'KneeReflex_circuit.pdf'));
exportgraphics(fig, fullfile(outDir, 'KneeReflex_circuit.svg'));
close(fig);
fprintf('Wrote figures/KneeReflex_circuit.{png,pdf,svg}\n');

%% ---- local helper: line + E/I marker at postsynaptic side ----
function drawConn(sx, sy, nx, ny, rSrc, rMk, angDeg, type)
    % marker sits at radius rMk from neuron center, along angDeg
    mx = nx + rMk*cosd(angDeg);
    my = ny + rMk*sind(angDeg);
    % trim line start at source circle edge
    d = [mx - sx, my - sy]; L = norm(d); u = d / L;
    sx2 = sx + u(1)*rSrc; sy2 = sy + u(2)*rSrc;
    snsfig('edge', sx2, sy2, mx, my);
    % marker points INTO the neuron
    angInto = atan2d(ny - my, nx - mx);
    switch type
        case 'exc'
            snsfig('exc', mx, my, angInto, 2.6);
        case 'inh'
            snsfig('inh', mx, my, 2.6);
        otherwise
            error('unknown contact type %s', type);
    end
end
