function sns_demo_figures_v4(which)
% V4 dissertation/defense demo figures (Ben 2026-10-26 round-2 spec):
%  - synapse markers at the END of each connection (postsynaptic side)
%  - color coding: muscles/BPAs salmon-red; afferent feedback BLUE dashed;
%    neural efferent BLACK solid; consistent legend in every figure
%  - Beer: mechanical model shown TWICE (empty cup, start level; full cup,
%    sagged), muscles visible along the arm, ON vs OFF sag clearly
%    annotated (final values printed + labeled on the plot)
%  - KneeCPG: cyclic kinematics (2 s window + full-run range note) + a
%    leg-motion strip (four poses across a cycle) + architecture + the
%    four-corner half-center core
%  - KneeReflex: leg mechanism (explicit axes) + kinematics + cross-
%    software comparison (Simulink vs numpy twin of the same equations)

if nargin < 1, which = 'all'; end
here = fileparts(mfilename('fullpath'));
cd(here);
res = fullfile(here, '..', 'results');
pics = fullfile(res, 'pictures', 'v4');
if ~exist(pics, 'dir'), mkdir(pics); end

CLR = struct('muscle', [0.90 0.35 0.30], 'affer', [0.00 0.45 0.74], ...
    'effer', [0.10 0.10 0.10], 'mech', [0.55 0.55 0.55]);

if any(strcmp(which, {'all', 'knee'}))
    s = load(fullfile(res, 'sns_demo_results.mat'));
    tw = load(fullfile(res, 'knee_twin_ref.mat'));
    knee_figure(s.out.log_th, tw, CLR, fullfile(pics, 'fig_kneereflex_v4.png'));
end
if any(strcmp(which, {'all', 'cpg'}))
    s = load(fullfile(res, 'sns_cpg_results.mat'));
    cpg_figure(s.out.log_th, CLR, fullfile(pics, 'fig_kneecpg_v4.png'));
end
if any(strcmp(which, {'all', 'beer'}))
    s = load(fullfile(res, 'sns_beer_results.mat'));
    beer_figure(s.runs, CLR, fullfile(pics, 'fig_beer_v4.png'));
end
end

%% ---------------- glyphs ----------------
function capsule(ax, p1, p2, w1, w2, col)
d = p2 - p1; L = norm(d); u_ = d / L; n_ = [-u_(2), u_(1)];
c1 = p1 + n_ * w1 / 2; c2 = p2 + n_ * w2 / 2;
c3 = p2 - n_ * w2 / 2; c4 = p1 - n_ * w1 / 2;
patch(ax, 'XData', [c1(1) c2(1) c3(1) c4(1)], ...
      'YData', [c1(2) c2(2) c3(2) c4(2)], ...
      'FaceColor', col, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1.0);
th_ = linspace(0, 2 * pi, 24)';
patch(ax, 'XData', p1(1) + w1 / 2 * cos(th_), 'YData', p1(2) + w1 / 2 * sin(th_), ...
      'FaceColor', col, 'EdgeColor', 'none');
patch(ax, 'XData', p2(1) + w2 / 2 * cos(th_), 'YData', p2(2) + w2 / 2 * sin(th_), ...
      'FaceColor', col, 'EdgeColor', 'none');
end

function joint(ax, p, r)
th_ = linspace(0, 2 * pi, 24)';
patch(ax, 'XData', p(1) + r * cos(th_), 'YData', p(2) + r * sin(th_), ...
      'FaceColor', [0.80 0.80 0.84], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
end

function bpa(ax, p1, p2, act, lw)
if nargin < 5, lw = 0.028; end
d = p2 - p1; L = norm(d); u_ = d / L; n_ = [-u_(2), u_(1)];
col = [1.0, 0.55 + 0.25 * (1 - act), 0.45];
t = linspace(0.15, 0.85, 18)';
bx = lw / 2 * sin(pi * (t - 0.15) / 0.7);
top = p1 + u_ .* (L * t) + n_ .* bx;
bot = p1 + u_ .* (L * t) - n_ .* bx;
patch(ax, 'XData', [top(:, 1); flipud(bot(:, 1))], ...
      'YData', [top(:, 2); flipud(bot(:, 2))], ...
      'FaceColor', col, 'EdgeColor', [0.35 0.2 0.15], 'LineWidth', 0.8);
for k = 1:9
    s1 = 2 * (k / 10) - 1; s2 = -(2 * ((k + 0.5) / 10) - 1);
    q1 = p1 + u_ * (L * (0.15 + 0.7 * k / 10)) + n_ * (lw / 2 * s1);
    q2 = p1 + u_ * (L * (0.15 + 0.7 * (k + 0.5) / 10)) + n_ * (lw / 2 * s2);
    plot(ax, [q1(1) q2(1)], [q1(2) q2(2)], 'Color', [0.45 0.3 0.2], 'LineWidth', 0.6);
end
capsule(ax, p1, p1 + u_ * 0.018, 0.018, 0.018, [0.65 0.65 0.68]);
capsule(ax, p2 - u_ * 0.018, p2, 0.018, 0.018, [0.65 0.65 0.68]);
end

function draw_leg(ax, thd, ghost)
% Leg with muscles routed OUTSIDE the limb silhouette; ghost = motion strip
if nargin < 3, ghost = false; end
if ghost, face = [1 1 1]; edge = [0.6 0.6 0.6]; lwF = 0.7; else
    face = [0.92 0.92 0.94]; edge = [0.25 0.25 0.25]; lwF = 1.0; end
Lf = 0.40; Ls = 0.42;
hip = [0, 0]; knee = [0, -Lf];
a = deg2rad(-90 - thd);
su = [cos(a), sin(a)]; ankle = knee + Ls * su;
nrm = [-su(2), su(1)];
capsule(ax, hip, knee, 0.115, 0.085, face); hold(ax, 'on');
capsule(ax, knee, ankle, 0.075, 0.050, face * 0.96);
capsule(ax, ankle - 0.035 * nrm + 0.045 * su, ankle + 0.10 * nrm + 0.045 * su, ...
    0.052, 0.058, face * 0.96);
joint(ax, hip, 0.030); joint(ax, knee, 0.026); joint(ax, ankle, 0.022);
if ghost, return; end
qO = [0.075, -0.10]; qI = ankle - 0.36 * su + 0.030 * nrm;
pat = knee + 0.028 * nrm;
bpa(ax, qO, pat, 0.6, 0.030);
plot(ax, [pat(1) qI(1)], [pat(2) qI(2)], 'Color', [0.85 0.82 0.70], 'LineWidth', 3);
plot(ax, qO(1), qO(2), 'k.', 'MarkerSize', 12);
plot(ax, qI(1), qI(2), 'k.', 'MarkerSize', 12);
hO = [-0.075, -0.10]; hI = ankle - 0.34 * su - 0.030 * nrm;
bpa(ax, hO, hI, 0.35, 0.026);
plot(ax, hO(1), hO(2), 'k.', 'MarkerSize', 12);
plot(ax, hI(1), hI(2), 'k.', 'MarkerSize', 12);
text(ax, 0.10, -0.13, 'extensor BPA', 'FontSize', 9, 'Color', [0.55 0.15 0.10]);
text(ax, -0.19, -0.13, 'flexor BPA', 'FontSize', 9, 'Color', [0.55 0.15 0.10]);
text(ax, knee(1) + 0.05, knee(2) - 0.02, sprintf('knee \\theta = %.0f°', thd), 'FontSize', 9);
axis(ax, 'equal'); xlim(ax, [-0.40 0.30]); ylim(ax, [-0.95 0.10]); axis(ax, 'off');
end

function neuron_glyph(ax, p, r, lbl)
th_ = linspace(0, 2 * pi, 20)';
patch(ax, 'XData', p(1) + r * cos(th_), 'YData', p(2) + r * sin(th_), ...
      'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.3);
if nargin > 3 && ~isempty(lbl)
    text(ax, p(1), p(2), lbl, 'HorizontalAlignment', 'center', 'FontSize', 8.5);
end
end

function connect_glyph(ax, p1, p2, kind, CLR)
% kind: 'effer' (black), 'affer' (blue dashed) ; synapse marker AT THE END
% (on the postsynaptic neuron), exc = open triangle, inh = filled circle.
if nargin < 5, CLR = struct('affer', [0 0.45 0.74], 'effer', [0.1 0.1 0.1]); end %#ok<INUSD>
if strcmp(kind, 'affer')
    col = [0.00 0.45 0.74]; lst = '--';
else
    col = [0.10 0.10 0.10]; lst = '-';
end
plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], 'Color', col, 'LineStyle', lst, 'LineWidth', 1.4);
end

function synapse_marker(ax, p1, p2, exc)
% marker at the END of the p1->p2 synapse, just outside the target.
d = p2 - p1; u_ = d / norm(d);
m = p2 - u_ * 0.075;                    % just before the target soma
if exc
    w = u_ * 0.030; n_ = [-u_(2), u_(1)] * 0.020;
    patch(ax, 'XData', [m(1) + w(1), m(1) - w(1) + n_(1), m(1) - w(1) - n_(1)], ...
          'YData', [m(2) + w(2), m(2) - w(2) + n_(2), m(2) - w(2) - n_(2)], ...
          'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.0);
else
    th_ = linspace(0, 2 * pi, 12)';
    patch(ax, 'XData', m(1) + 0.018 * cos(th_), 'YData', m(2) + 0.018 * sin(th_), ...
          'FaceColor', [0.1 0.1 0.1], 'EdgeColor', 'k');
end
end

function colorkey(ax, CLR) %#ok<INUSD>
plot(ax, nan, nan, 'Color', [0.90 0.35 0.30], 'LineWidth', 6);
hold(ax, 'on');
plot(ax, nan, nan, 'Color', [0.00 0.45 0.74], 'LineStyle', '--', 'LineWidth', 1.8);
plot(ax, nan, nan, 'Color', [0.10 0.10 0.10], 'LineWidth', 1.8);
lg = legend(ax, {'muscle / BPA', 'afferent feedback (Ia/Ib)', ...
    'efferent neural pathway'}, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'FontSize', 8);
lg.Box = 'off';
end

%% ---------------- Figure: knee reflex ----------------
function knee_figure(th, twin, CLR, fname) %#ok<INUSD>
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 1500 520], ...
    'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% A kinematics + cross-software comparison
ax = nexttile(tl);
t = th.Time(:); v = rad2deg(th.Data(:));
plot(ax, t, v, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8); hold(ax, 'on');
tv = twin.t(:); vv = rad2deg(twin.theta(:));
plot(ax, tv, vv, 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.4);
vvi = interp1(tv, vv, t, 'linear', 'extrap');
d = max(abs(v - vvi));
legend(ax, {'Simulink (SNS\_Library)', 'numpy twin (same equations)'}, ...
    'Location', 'southeast', 'FontSize', 8.5);
grid(ax, 'on');
xlabel(ax, 'time (s)'); ylabel(ax, 'knee flexion \theta (deg)');
title(ax, sprintf('A  Kinematics: settle %.1f° (twin %.1f°, max dev %.1f°)', ...
    mean(v(end - 200:end)), vv(end), d), 'FontWeight', 'bold');
% B mechanism with leg
ax = nexttile(tl);
draw_leg(ax, 30, false);
title(ax, 'B  Leg model + BPA routing', 'FontWeight', 'bold');
% C pathway, markers at ENDS, color-coded
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
ia_e = [0.28, 0.86]; ib_e = [0.28, 0.66]; ia_f = [0.28, 0.32];
mn_e = [0.75, 0.62]; mn_f = [0.75, 0.20];
neuron_glyph(ax, ia_e, 0.045, 'Ia_{ext}');
neuron_glyph(ax, ib_e, 0.045, 'Ib_{ext}');
neuron_glyph(ax, ia_f, 0.045, 'Ia_{flx}');
neuron_glyph(ax, mn_e, 0.055, 'MN_{ext}');
neuron_glyph(ax, mn_f, 0.055, 'MN_{flx}');
connect_glyph(ax, ia_e, mn_e, 'affer'); synapse_marker(ax, ia_e, mn_e, true);
connect_glyph(ax, ib_e, mn_e, 'affer'); synapse_marker(ax, ib_e, mn_e, false);
connect_glyph(ax, ia_f, mn_e, 'affer'); synapse_marker(ax, ia_f, mn_e, false);
connect_glyph(ax, ia_f, mn_f, 'affer'); synapse_marker(ax, ia_f, mn_f, true);
text(ax, 0.02, 0.95, 'stretch (Ia, exc)', 'FontSize', 8, 'Color', [0 0.45 0.74]);
text(ax, 0.02, 0.74, 'autogenic (Ib, inh)', 'FontSize', 8, 'Color', [0 0.45 0.74]);
text(ax, 0.02, 0.40, 'reciprocal (Ia, inh)', 'FontSize', 8, 'Color', [0 0.45 0.74]);
title(ax, 'C  Reflex pathway', 'FontWeight', 'bold');
exportgraphics(fig, fname, 'Resolution', 300);
close(fig);
fprintf('wrote %s (twin settle %.2f deg, max dev %.2f deg)\n', fname, vv(end), d);
end

%% ---------------- Figure: knee CPG ----------------
function cpg_figure(th, CLR, fname) %#ok<INUSD>
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 1500 560], ...
    'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% A cyclic kinematics: 2 s window (cycles legible) + full-run note
ax = nexttile(tl);
t = th.Time(:); v = rad2deg(th.Data(:));
w = t >= 3 & t <= 5;
plot(ax, t(w), v(w), 'Color', [0.64 0.08 0.18], 'LineWidth', 1.4); grid(ax, 'on');
ipk = find(islocalmax(v, 'MinProminence', 1));
ipk = ipk(t(ipk) > 1);
if numel(ipk) >= 2
    per = mean(diff(t(ipk)));
else
    % robust fallback: zero crossings of the mean-removed trace
    z = v - mean(v);
    ix = find(z(1:end - 1) .* z(2:end) < 0 & t(1:end - 1) > 0.5);
    per = 2 * mean(diff(t(ix)));
end
xlabel(ax, 'time (s)'); ylabel(ax, 'knee flexion \theta (deg)');
if isnan(per)
    title(ax, sprintf('A  Cyclic response: range %.0f-%.0f° (hysteresis switching)', ...
        min(v), max(v)), 'FontWeight', 'bold');
else
    title(ax, sprintf('A  Cyclic response: period %.3f s (%.1f Hz), range %.0f-%.0f°', ...
        per, 1 / per, min(v), max(v)), 'FontWeight', 'bold');
end
% B leg motion strip: four poses across one cycle
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [-0.55 0.45]); ylim(ax, [-0.98 0.12]);
poses = [12, 29, 46, 29];
for k = 1:4
    draw_leg_at(ax, poses(k) + (k - 1) * 0.10, -0.18 * (k - 1) + 0.06, ...
        k == 2 || k == 4);
end
text(ax, -0.52, 0.10, 'one cycle: flexion swings the shank posteriorly', ...
    'FontSize', 8.5);
title(ax, 'B  Leg motion (one cycle)', 'FontWeight', 'bold');
% C architecture + four-corner inset
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
blks = {0.06, 'RG'; 0.30, 'PF'; 0.54, 'MN'; 0.78, 'BPA'};
for k = 1:4
    rectangle(ax, 'Position', [blks{k, 1}, 0.70, 0.16, 0.13], ...
        'Curvature', 0.3, 'EdgeColor', 'k', 'FaceColor', [0.94 0.94 0.96]);
    text(ax, blks{k, 1} + 0.08, 0.765, blks{k, 2}, 'HorizontalAlignment', ...
        'center', 'FontSize', 9, 'FontWeight', 'bold');
    if k < 4
        plot(ax, [blks{k, 1} + 0.16, blks{k + 1, 1}], [0.765 0.765], 'k-', ...
            'LineWidth', 1.3);
    end
end
text(ax, 0.30, 0.93, 'architecture: RG \rightarrow PF \rightarrow MN \rightarrow BPA \rightarrow knee', ...
    'FontSize', 8.5);
% four-corner core below
hE = [0.26, 0.44]; hF = [0.74, 0.44];
iE = [0.26, 0.14]; iF = [0.74, 0.14];
neuron_glyph(ax, hE, 0.05, 'HC_E');
neuron_glyph(ax, hF, 0.05, 'HC_F');
neuron_glyph(ax, iE, 0.045, 'IN_E');
neuron_glyph(ax, iF, 0.045, 'IN_F');
connect_glyph(ax, hE, iE, 'effer'); synapse_marker(ax, hE, iE, true);
connect_glyph(ax, hF, iF, 'effer'); synapse_marker(ax, hF, iF, true);
connect_glyph(ax, iE, hF, 'effer'); synapse_marker(ax, iE, hF, false);
connect_glyph(ax, iF, hE, 'effer'); synapse_marker(ax, iF, hE, false);
plot(ax, [hE(1) + 0.045, hF(1) - 0.045], [hE(2) - 0.03, hF(2) - 0.03], 'k:', ...
    'LineWidth', 1.0);
text(ax, 0.50, 0.335, 'weak cross-excitation', 'HorizontalAlignment', 'center', ...
    'FontSize', 7.5);
title(ax, 'C  Architecture + half-center core', 'FontWeight', 'bold');
exportgraphics(fig, fname, 'Resolution', 300);
close(fig);
fprintf('wrote %s (period %.4f s)\n', fname, per);
end

function draw_leg_at(ax, thd, dx, dy)
Lf = 0.40; Ls = 0.42;
hip = [dx, dy]; knee = hip + [0, -Lf];
a = deg2rad(-90 - thd);
su = [cos(a), sin(a)]; ankle = knee + Ls * su;
nrm = [-su(2), su(1)];
capsule(ax, hip, knee, 0.105, 0.080, [0.99 0.99 1.0]);
capsule(ax, knee, ankle, 0.070, 0.046, [0.99 0.99 1.0]);
capsule(ax, ankle - 0.035 * nrm + 0.045 * su, ankle + 0.10 * nrm + 0.045 * su, ...
    0.050, 0.056, [0.99 0.99 1.0]);
plot(ax, [hip(1) knee(1) ankle(1)], [hip(2) knee(2) ankle(2)], ...
    'Color', [0.35 0.35 0.35], 'LineWidth', 1.6);
end

%% ---------------- Figure: beer ----------------
function beer_figure(runs, CLR, fname) %#ok<INUSD>
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 1680 560], ...
    'Color', 'w');
tl = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% Panel 1: mechanism, EMPTY cup, forearm level (start)
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [-0.40 0.40]); ylim(ax, [-0.85 0.15]);
draw_arm(ax, 0, 0.0, 0.35);
text(ax, -0.36, 0.08, 'start: empty cup, forearm level', 'FontSize', 8.5);
title(ax, 'B1  Empty cup (start)', 'FontWeight', 'bold');
% Panel 2: mechanism, FULL cup, sagged (reflex ON final; OFF shown in A)
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [-0.40 0.40]); ylim(ax, [-0.85 0.15]);
sagON = rad2deg(runs(1).th.Data(end));
sagOFF = rad2deg(runs(2).th.Data(end));
draw_arm(ax, sagON, 0.5, 0.35);
text(ax, -0.36, 0.08, sprintf('end: full cup, sag %.1f° (ON)', sagON), 'FontSize', 8.5);
title(ax, 'B2  Full cup (end)', 'FontWeight', 'bold');
% Panel 3: kinematics ON vs OFF, ELEVATION convention (negative = below
% level) so the physically higher glass plots HIGHER (Ben 2026-10-26:
% positive-up sag read as "OFF raises the glass"). ON = -2.0 deg sits
% above OFF = -7.8 deg; more droop = lower line = lower glass.
ax = nexttile(tl);
t1 = runs(1).th.Time(:); v1 = -rad2deg(runs(1).th.Data(:));
t2 = runs(2).th.Time(:); v2 = -rad2deg(runs(2).th.Data(:));
hold(ax, 'on'); grid(ax, 'on');
plot(ax, t1, v1, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.8);
plot(ax, t2, v2, 'Color', [0.55 0.55 0.55], 'LineStyle', '--', 'LineWidth', 1.6);
text(ax, t1(end) + 0.15, v1(end), sprintf('  ON %.2f° (glass held higher)', v1(end)), ...
    'FontSize', 9, 'Color', [0.00 0.45 0.74]);
text(ax, t2(end) + 0.15, v2(end), sprintf('  OFF %.2f° (droops)', v2(end)), ...
    'FontSize', 9, 'Color', [0.35 0.35 0.35]);
yline(ax, 0, 'k:', 'level');
ylim(ax, [floor(min(v2) * 1.25) 1]);
xlabel(ax, 'time (s)'); ylabel(ax, 'forearm elevation (deg; 0 = level)');
title(ax, sprintf('A  Sag: ON droops %.2f° vs OFF %.2f° (%.1fx stiffer)', ...
    abs(v1(end)), abs(v2(end)), abs(v2(end) / v1(end))), 'FontWeight', 'bold');
% Panel 4: pathway, markers at ENDS, color-coded feedback vs efferent
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
ia_bi = [0.26, 0.86]; ib_bi = [0.26, 0.66]; ia_tri = [0.26, 0.32];
mn_bi = [0.76, 0.66]; mn_tri = [0.76, 0.22];
neuron_glyph(ax, ia_bi, 0.045, 'Ia_{bi}');
neuron_glyph(ax, ib_bi, 0.045, 'Ib_{bi}');
neuron_glyph(ax, ia_tri, 0.045, 'Ia_{tri}');
neuron_glyph(ax, mn_bi, 0.055, 'MN_{bi}');
neuron_glyph(ax, mn_tri, 0.055, 'MN_{tri}');
connect_glyph(ax, ia_bi, mn_bi, 'affer'); synapse_marker(ax, ia_bi, mn_bi, true);
connect_glyph(ax, ib_bi, mn_bi, 'affer'); synapse_marker(ax, ib_bi, mn_bi, false);
connect_glyph(ax, ia_tri, mn_bi, 'affer'); synapse_marker(ax, ia_tri, mn_bi, false);
connect_glyph(ax, ia_tri, mn_tri, 'affer'); synapse_marker(ax, ia_tri, mn_tri, true);
connect_glyph(ax, ia_bi, mn_tri, 'affer'); synapse_marker(ax, ia_bi, mn_tri, false);
text(ax, 0.02, 0.95, 'feedback = blue dashed; markers at the' , 'FontSize', 7.5, ...
    'Color', [0 0.45 0.74]);
text(ax, 0.02, 0.90, 'postsynaptic ending', 'FontSize', 7.5, 'Color', [0 0.45 0.74]);
title(ax, 'C  Reflex pathway', 'FontWeight', 'bold');
exportgraphics(fig, fname, 'Resolution', 300);
close(fig);
fprintf('wrote %s (ON %.3f deg, OFF %.3f deg)\n', fname, v1(end), v2(end));
end

function draw_arm(ax, sag, fillFrac, scale)
% Natural posture: humerus VERTICAL, forearm horizontal; sag rotates the
% forearm below level. Muscles drawn ALONG the upper arm (anterior biceps,
% posterior triceps), routed outside the limb silhouette.
s = scale;
shoulder = [0, 0];
elbow = shoulder + s * [0, -1];                       % humerus down
fa = deg2rad(sag);
hand = elbow + s * [cos(-fa), sin(-fa)];
capsule(ax, shoulder, elbow, 0.10 * s, 0.08 * s, [0.92 0.92 0.94]);
capsule(ax, elbow, hand, 0.075 * s, 0.055 * s, [0.88 0.88 0.91]);
capsule(ax, shoulder + [-0.10 * s, 0.03], shoulder + [0.10 * s, 0.03], ...
    0.11 * s, 0.11 * s, [0.90 0.90 0.93]);            % shoulder block
joint(ax, shoulder, 0.030);
joint(ax, elbow, 0.026);
% biceps: anterior upper-arm compartment -> forearm insertion
uU = (elbow - shoulder) / norm(elbow - shoulder);
nU = [-uU(2), uU(1)];                                  % anterior side
biO = shoulder + 0.12 * uU + 0.075 * nU;
biI = elbow + 0.09 * (hand - elbow) / norm(hand - elbow) + 0.05 * nU;
bpa(ax, biO, biI, 0.45, 0.028);
plot(ax, biO(1), biO(2), 'k.', 'MarkerSize', 12);
plot(ax, biI(1), biI(2), 'k.', 'MarkerSize', 12);
% triceps: posterior upper-arm compartment -> olecranon
trO = shoulder + 0.10 * uU - 0.075 * nU;
trI = elbow - 0.03 * nU;
bpa(ax, trO, trI, 0.40, 0.024);
plot(ax, trO(1), trO(2), 'k.', 'MarkerSize', 12);
plot(ax, trI(1), trI(2), 'k.', 'MarkerSize', 12);
text(ax, biO(1) + 0.02, biO(2) + 0.05, 'biceps', 'FontSize', 8, ...
    'Color', [0.55 0.15 0.10]);
text(ax, trO(1) - 0.09, trO(2) + 0.02, 'triceps', 'FontSize', 8, ...
    'Color', [0.55 0.15 0.10]);
% pint glass with level beer surface
gn = [cos(fa + pi / 2), sin(fa + pi / 2)];
uf = (hand - elbow) / norm(hand - elbow);
base = hand + 0.015 * gn;
b1 = base - 0.075 * uf; b2 = base + 0.075 * uf;
r1 = b1 + 0.17 * gn; r2 = b2 + 0.17 * gn;
patch(ax, 'XData', [b1(1) b2(1) r2(1) r1(1)], 'YData', [b1(2) b2(2) r2(2) r1(2)], ...
    'FaceColor', [0.94 0.96 1.0], 'EdgeColor', [0.55 0.62 0.70], ...
    'LineWidth', 1.6, 'FaceAlpha', 0.55);
if fillFrac > 0
    yF = base(2) + (0.08 + 0.56 * fillFrac) * 0.17 * 1.0 + 0.0;  % approx along axis
    yF = base(2) + (0.08 + 0.56 * fillFrac) * 0.17;
    wL = seg_pt(b1, r1, yF); wR = seg_pt(b2, r2, yF);
    if ~any(isnan(wL)) && ~any(isnan(wR))
        patch(ax, 'XData', [wL(1) wR(1) b2(1) b1(1)], ...
              'YData', [wL(2) wR(2) b2(2) b1(2)], ...
              'FaceColor', [0.96 0.74 0.16], 'EdgeColor', 'none');
        plot(ax, [wL(1) wR(1)], [wL(2) wR(2)], 'w-', 'LineWidth', 3);
    end
end
end

function p = seg_pt(a, b_, y)
t = (y - a(2)) / (b_(2) - a(2));
if t < 0 || t > 1, p = [NaN NaN]; else, p = a + t * (b_ - a); end
end
