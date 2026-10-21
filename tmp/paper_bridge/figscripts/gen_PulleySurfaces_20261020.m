% gen_PulleySurfaces_20261020.m -- EDIT_SPEC section-5 figures 1-2:
%   PulleyTorqueSurfaces_20261020.eps  |tau|(theta, p) contourf per site
%                                     (G = 1/1, recorded sizes) + the
%                                     |tau| = |tau_req| zero-margin curve.
%   PulleyMarginSurfaces_20261020.eps  margin (tau_r - tau_req)/tau_req vs
%                                     theta for gains {1/1, 2/1, 3/1} at
%                                     620 kPa; RoM shading where margin >= 0.
% Campaign style per EDIT_SPEC 20261020 section 0 (Arial 12 bold, LineWidth
% 2, PlotPrompt palette, legend 8 pt, panel letters on the title line).
% QC PNGs (150 dpi) to tmp/paper_bridge/figscripts.
root = 'D:/GitHub/Bipedal_Robot';
resDir = fullfile(root, 'Code', 'Matlab', 'Mesh_Optimization', 'Results');
outdir = fullfile(root, 'Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/Figures');
qcdir = fullfile(root, 'tmp/paper_bridge/figscripts');

c1 = [1 0.843 0];          % gold
c3 = [0.980 0.529 0.459];  % light orange  - gain 1/1 (tackle disabled)
c5 = [0.804 0.204 0.710];  % magenta       - gain 2/1
c7 = [0 0 1];              % blue          - gain 3/1
c8 = [0 0 0];
cgray = [0.4 0.4 0.4];

sty = {'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial', ...
    'LineWidth', 2, 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickLength', [0.025 0.05], 'GridLineStyle', 'none'};
gains = {'1/1', '2/1', '3/1'};
gcol = {c3, c5, c7};
glst = {'-', '--', '-.'};
gmrk = {'o', 's', '^'};

sites = {'soleus_r', 'tib_ant_r', 'bifemsh_r', 'med_gas_r', 'rect_fem_r'};
jtik = struct('soleus_r', [-50 -25 0 20], 'tib_ant_r', [-50 -25 0 20], ...
    'bifemsh_r', [-120 -90 -60 -30 0 10], 'med_gas_r', [-120 -90 -60 -30 0 10], ...
    'rect_fem_r', [-15 16 48 79 110]);
D = struct();
for s = 1:numel(sites)
    D.(sites{s}) = load(fullfile(resDir, sprintf('PulleySurfaces_%s_20261020.mat', sites{s})));
end

%% ================= FIG 1: torque surfaces tau(theta, p) =====================
fig1 = figure('Name', 'PulleyTorqueSurfaces', 'Color', 'w');
fig1.Position = [60 40 1150 760];
t1 = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose');
ttl1 = {'Soleus (ankle)', 'Tibialis anterior (ankle)', 'Biceps femoris sh. (knee)', ...
    'Medial gastrocnemius (knee)', 'Rectus femoris (hip)'};
axP1 = gobjects(1, 5);
for s = 1:5
    nexttile(s); hold on
    d = D.(sites{s});
    Z = abs(d.tauSurfP);              % N1 x nP (G = 1/1)
    X = d.theta1Slice_deg; Y = d.pressSurf;
    contourf(X, Y, Z', 14, 'LineColor', 'none');
    colormap(parula); cb = colorbar; cb.Label.String = '|tau|, N{\cdot}m';
    cb.Label.FontSize = 12; cb.Label.FontWeight = 'bold'; cb.Label.FontName = 'Arial';
    cb.FontSize = 11; cb.FontName = 'Arial'; cb.FontWeight = 'bold';
    % zero-margin curve: |tau(theta,p)| = |tau_req(theta)|
    req = abs(d.tauReqSlice);
    F = Z - req;        % N1 x nP broadcast over pressure
    contour(X, Y, F', [0 0], 'k--', 'LineWidth', 2);
    title(['\bf ' ttl1{s}]); set(gca, sty{:});
    axP1(s) = gca;
    xlim([min(X) max(X)]); axP1(s).XTick = jtik.(sites{s});
    axP1(s).FontSize = 11;
end
% info tile
nexttile(6); axis off; hold on
txt = { ...
    '\bfTorque surfaces, gain 1/1', ...
    'contour: |\tau| (\theta, p), recorded sizes', ...
    'knee / ankle / hip flexion negative', ...
    'black dashed: |\tau| = |\tau_{req}| (human target)', ...
    'above the dashed curve: target met', ...
    '', ...
    'moving_via routing (bowden never used)', ...
    'sources: PulleySurfaces_*_20261020.mat'};
text(0.02, 0.97, txt, 'Units', 'normalized', 'FontSize', 11, ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'VerticalAlign', 'top');
ylabel(t1, '\bf Pressure, kPa'); xlabel(t1, '\bf joint angle, \circ');
set(t1, 'TileSpacing', 'loose');
drawnow;
letters = {'(A)', '(B)', '(C)', '(D)', '(E)'};
xp = [0.003 0.335 0.668 0.003 0.335];
for j = 1:5
    axj = axP1(j); p4 = axj.Position;
    tcY = 0.935;
    try
        Tj = axj.Title; Tj.Units = 'normalized';
        te = Tj.Extent;
        cand = p4(2) + (te(2) + te(4) / 2) * p4(4);
        if isfinite(cand) && cand > 0.3, tcY = cand; end
    catch
    end
    a1 = annotation(fig1, 'textbox', [xp(j) tcY - 0.025 0.05 0.05], ...
        'String', ['\bf ' letters{j}], 'FontSize', 12, 'FontName', 'Arial', ...
        'FontWeight', 'bold', 'EdgeColor', 'none');
    try
        drawnow; ae = a1.Extent;
        a1.Position(2) = a1.Position(2) + (tcY - (ae(2) + ae(4) / 2));
    catch
    end
end
exportgraphics(fig1, fullfile(outdir, 'PulleyTorqueSurfaces_20261020.eps'), 'ContentType', 'vector');
exportgraphics(fig1, fullfile(qcdir, 'QC_PulleyTorqueSurfaces_20261020.png'), 'Resolution', 150);
fprintf('WROTE PulleyTorqueSurfaces_20261020.eps\n'); close(fig1);

%% ================= FIG 2: margin surfaces over (theta x gain) ==============
fig2 = figure('Name', 'PulleyMarginSurfaces', 'Color', 'w');
fig2.Position = [60 40 1150 760];
t2 = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose');
axP2 = gobjects(1, 5);
legDone = false;
for s = 1:5
    nexttile(s); hold on
    d = D.(sites{s});
    X = d.theta1Slice_deg;
    pp = 3;   % 620 kPa
    yl = [-105 320];
    % shade the union region where any gain meets the target (drawn first)
    okAny = any(reshape(isfinite(d.marginSlice(:, :, pp)) & ...
        d.marginSlice(:, :, pp) >= 0, N1r(d), 3), 2);
    area(X, okAny * (yl(2) - yl(1)) + yl(1), 'BaseValue', yl(1), ...
        'FaceColor', 0.65 + 0.35 * c1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    for g = 1:3
        m = d.marginSlice(:, g, pp);
        mm = m; mm(~isfinite(mm)) = -105;
        mm = min(100 * mm, 300);           % display clip at +300%
        plot(X, mm, glst{g}, 'Color', gcol{g}, 'LineWidth', 2, ...
            'DisplayName', ['gain ' gains{g}]);
        mk = 1:3:numel(X);
        plot(X(mk), mm(mk), gmrk{g}, 'Color', gcol{g}, ...
            'MarkerFaceColor', gcol{g}, 'MarkerSize', 4.5, ...
            'HandleVisibility', 'off');
    end
    plot(xlim, [0 0], 'k-', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    title(['\bf ' ttl1{s}]); set(gca, sty{:});
    ylim(yl); yticks([-100 -50 0 50 100 150 200 250 300]);
    axP2(s) = gca;
    xlim([min(X) max(X)]); axP2(s).XTick = jtik.(sites{s});
    axP2(s).FontSize = 11;
    if ~legDone
        lg = legend('Location', 'northwest'); lg.FontSize = 8; legDone = true;
    end
end
nexttile(6); axis off; hold on
txt2 = { ...
    '\bfMargin over human target, 620 kPa', ...
    'margin = (|\tau_{robot}| - |\tau_{req}|) / |\tau_{req}|', ...
    'gains = tackle line parts at crossing 1', ...
    'sizes fixed at the recorded 1/1 optimum', ...
    'gold band: some gain meets the target', ...
    '-105%: infeasible or zero-torque cell', ...
    'curves clipped at +300%', ...
    '', ...
    '\tau_{req}: OpenSim Gait2392 (soleus, tib ant.,', ...
    'med. gas, rect. fem.; moment arm x fmax);', ...
    'bifemsh: recorded mono-run human torque'};
text(0.02, 0.97, txt2, 'Units', 'normalized', 'FontSize', 11, ...
    'FontName', 'Arial', 'FontWeight', 'bold', 'VerticalAlign', 'top');
ylabel(t2, '\bf margin, %'); xlabel(t2, '\bf joint angle, \circ');
set(t2, 'TileSpacing', 'loose');
drawnow;
for j = 1:5
    axj = axP2(j); p4 = axj.Position;
    tcY = 0.935;
    try
        Tj = axj.Title; Tj.Units = 'normalized';
        te = Tj.Extent;
        cand = p4(2) + (te(2) + te(4) / 2) * p4(4);
        if isfinite(cand) && cand > 0.3, tcY = cand; end
    catch
    end
    a1 = annotation(fig2, 'textbox', [xp(j) tcY - 0.025 0.05 0.05], ...
        'String', ['\bf ' letters{j}], 'FontSize', 12, 'FontName', 'Arial', ...
        'FontWeight', 'bold', 'EdgeColor', 'none');
    try
        drawnow; ae = a1.Extent;
        a1.Position(2) = a1.Position(2) + (tcY - (ae(2) + ae(4) / 2));
    catch
    end
end
exportgraphics(fig2, fullfile(outdir, 'PulleyMarginSurfaces_20261020.eps'), 'ContentType', 'vector');
exportgraphics(fig2, fullfile(qcdir, 'QC_PulleyMarginSurfaces_20261020.png'), 'Resolution', 150);
fprintf('WROTE PulleyMarginSurfaces_20261020.eps\n'); close(fig2);

function n = N1r(d)
n = numel(d.theta1Slice_deg);
end
