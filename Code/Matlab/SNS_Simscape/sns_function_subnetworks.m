%% sns_function_subnetworks.m — function-subnetwork primitives panel (journal figure)
%
% Six canonical SNS subnetworks that perform arithmetic on graded (non-spiking)
% signals, after Szczecinski et al. 2017, "A Functional Subnetwork Approach to
% Designing Synthetic Nervous Systems That Control Legged Robot Locomotion"
% (Front. Neurorobotics 11:37, DOI 10.3389/fnbot.2017.00037; Fig. 2 networks
% A-D and Sec. 4 for differentiation/integration), extended to the generalized
% linear integrate-and-fire neuron by Szczecinski & co. 2020
% (DOI 10.3389/fnbot.2020.577804). Drawn in our diagram language (same as
% sns_draw_circuit.m): open triangle = excitatory, solid dot = inhibitory,
% dashed = modulatory pathway. Shape-only coding = colorblind safe.
% Outputs: figures/SNS_function_subnetworks.{png,pdf,svg}

cdto = fileparts(mfilename('fullpath'));
cd(cdto);
outDir = fullfile(cdto, 'figures');
if ~exist(outDir, 'dir'), mkdir(outDir); end

fig = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 18 11], 'Color', 'w');
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

panels = {'(a) addition', '(b) subtraction', '(c) division', ...
          '(d) multiplication', '(e) differentiation', '(f) integration'};
for k = 1:6
    ax = nexttile(tl);
    hold(ax, 'on'); axis(ax, 'equal'); axis(ax, [0 30 0 20]); axis(ax, 'off');
    title(ax, panels{k}, 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
    switch k
        case 1, panelAdd(ax);
        case 2, panelSub(ax);
        case 3, panelDiv(ax);
        case 4, panelMul(ax);
        case 5, panelDiff(ax);
        case 6, panelInt(ax);
    end
end

exportgraphics(fig, fullfile(outDir, 'SNS_function_subnetworks.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(outDir, 'SNS_function_subnetworks.pdf'));
exportgraphics(fig, fullfile(outDir, 'SNS_function_subnetworks.svg'));
close(fig);
fprintf('Wrote figures/SNS_function_subnetworks.{png,pdf,svg}\n');

%% ---------------- panels ----------------
function panelAdd(ax)
    axes(ax);
    stub(3, 14, 'U_1');  syn(5.6, 14, 18, 10, 'exc');
    stub(3, 6,  'U_2');  syn(5.6, 6,  18, 10, 'exc');
    snsfig('neuron', 18, 10, 2.6, '');
    snsfig('arrow', 20.6, 10, 24.2, 10);
    snsfig('label', 26.6, 11.8, 'U_1 + U_2', 7.5, 'b');
    snsfig('label', 15, 1.6, 'two excitatory signal-transmission synapses', 6.5, 'gi');
end

function panelSub(ax)
    axes(ax);
    stub(3, 14, 'U_1');  syn(5.6, 14, 18, 10, 'exc');
    stub(3, 6,  'U_2');  syn(5.6, 6,  18, 10, 'inh');
    snsfig('neuron', 18, 10, 2.6, '');
    snsfig('arrow', 20.6, 10, 24.2, 10);
    snsfig('label', 26.6, 11.8, 'U_1 - U_2', 7.5, 'b');
    snsfig('label', 15, 1.6, 'inhibitory synapse weight = -k_{syn,1}\DeltaE_{s,1}/\DeltaE_{s,2}', 6.5, 'gi');
end

function panelDiv(ax)
    axes(ax);
    stub(3, 14, 'U_1');  syn(5.6, 14, 18, 10, 'exc');
    stub(3, 6,  'U_2');
    snsfig('edge', 5.6, 6, 14.7, 8.9, '--');
    snsfig('inh', 15.0, 9.0, 2.2);
    snsfig('neuron', 18, 10, 2.6, '');
    snsfig('arrow', 20.6, 10, 24.2, 10);
    snsfig('label', 25.4, 11.9, 'U_1/(1+c_{syn}U_2)', 6.5, 'b');
    snsfig('label', 15, 1.6, 'shunting: \DeltaE_{s,2} = 0; dashed = modulatory', 6.5, 'gi');
end

function panelMul(ax)
    axes(ax);
    % U1 excites the output neuron; U2 inhibits a tonically-active (Iapp = R)
    % mid neuron, which inhibits the output neuron (disinhibitory cascade).
    stub(3, 14.5, 'U_1');  syn(5.6, 14.5, 20, 12, 'exc');
    stub(3, 4, 'U_2');     syn(5.6, 4, 12.5, 5, 'inh');
    snsfig('neuron', 12.5, 5, 2.4, '');
    snsfig('label', 17.6, 2.7, 'I_{app} = R', 7, 'gi');
    snsfig('edge', 14.3, 6.6, 17.6, 9.8);
    snsfig('inh', 17.8, 10.0, 2.2);
    snsfig('neuron', 20, 12, 2.4, '');
    snsfig('arrow', 22.4, 12, 24.6, 12);
    snsfig('label', 27.2, 13.8, 'U_1U_2 / R', 7.5, 'b');
    snsfig('label', 15, 0.9, 'disinhibitory cascade: modulatory synapses in series', 6.5, 'gi');
end

function panelDiff(ax)
    axes(ax);
    % Same input into two neurons with different membrane capacities;
    % subtract fast copy (Cm,1) from slow copy (Cm,2) -> derivative (Reichardt).
    stub(3, 10, 'U');
    syn(5.6, 10, 16, 13.5, 'exc');
    syn(5.6, 10, 16, 6.5, 'exc');
    snsfig('neuron', 16, 13.5, 2.2, '');
    snsfig('neuron', 16, 6.5, 2.2, '');
    snsfig('label', 16, 16.3, 'C_{m,1}', 7, 'gi');
    snsfig('label', 16, 3.4, 'C_{m,2} > C_{m,1}', 7, 'gi');
    syn(18.2, 13.8, 23.5, 10, 'exc');
    syn(18.2, 6.2, 23.5, 10, 'inh');
    snsfig('neuron', 23.5, 10, 2.2, '');
    snsfig('arrow', 25.7, 10, 28, 10);
    snsfig('label', 27.4, 12.3, 'k_d dU/dt', 7, 'b');
    snsfig('label', 13, 1.6, 'subtract: k_d = C_{m,2} - C_{m,1}', 6.5, 'gi');
end

function panelInt(ax)
    axes(ax);
    % Self-disinhibition line attractor: signal u into U1, tonic R into U2,
    % mutual inhibition -> no leak; U1 ramps at k_i * u while driven.
    stub(3, 14, 'u');    syn(5.6, 14, 16, 13.5, 'exc');
    stub(3, 4, 'R');     syn(5.6, 4, 16, 5.5, 'exc');
    snsfig('neuron', 16, 13.5, 2.4, '');
    snsfig('neuron', 16, 5.5, 2.4, '');
    snsfig('edge', 14.2, 11.2, 14.2, 7.7);
    snsfig('inh', 14.2, 8.4, 2.2);
    snsfig('edge', 17.8, 7.8, 17.8, 11.3);
    snsfig('inh', 17.8, 10.6, 2.2);
    snsfig('arrow', 18.4, 13.5, 24, 13.5);
    snsfig('label', 24.8, 15.3, 'dU_1/dt = k_i u', 6.5, 'b');
    snsfig('label', 13, 1.6, 'mutual disinhibition = line attractor (no leak)', 6.5, 'gi');
end

%% ---------------- helpers ----------------
function stub(x, y, lbl)
    % short input arrow + bold signal label
    snsfig('arrow', x, y, x + 2.6, y);
    snsfig('label', x + 1.2, y + 1.5, lbl, 8, 'b');
end

function syn(sx, sy, nx, ny, type)
    % synaptic edge from (sx,sy) ending in an E/I marker near neuron (nx,ny)
    d = [nx - sx, ny - sy]; u = d / norm(d);
    rMk = 3.2;
    mx = nx - u(1)*rMk; my = ny - u(2)*rMk;
    snsfig('edge', sx, sy, mx - u(1)*0.3, my - u(2)*0.3);
    angInto = atan2d(ny - my, nx - mx);
    switch type
        case 'exc'
            snsfig('exc', mx, my, angInto, 2.2);
        case 'inh'
            snsfig('inh', mx, my, 2.2);
        otherwise
            error('unknown contact type %s', type);
    end
end
