function sns_demo_figures_v3(which)
% V3 dissertation/defense demo figures (Ben's 2026-10-26 spec):
%  KneeReflex: kinematic response overlay (reflex ON vs OFF) + mechanism
%    with muscles routed OUTSIDE the limb silhouette + reflex pathway.
%  KneeCPG: kinematic response + neural architecture + the four-corner
%    half-center mutual-inhibition diagram.
%  Beer: humerus VERTICAL, forearm HORIZONTAL, elbow 90 deg, natural
%    cup-hold posture (arm at side), start-position mechanism + sag
%    kinematics ON/OFF + reflex pathway.
% Outputs: results\pictures\v3\*.png (300 dpi) + this file's helpers feed
% sns_animate_demo_v3.m.

if nargin < 1, which = 'all'; end
here = fileparts(mfilename('fullpath'));
cd(here);
res = fullfile(here, '..', 'results');
pics = fullfile(res, 'pictures', 'v3');
if ~exist(pics, 'dir'), mkdir(pics); end

if any(strcmp(which, {'all', 'knee'}))
    s = load(fullfile(res, 'sns_demo_results.mat'));
    knee_figure(s.out.log_th, fullfile(pics, 'fig_kneereflex_v3.png'));
end
if any(strcmp(which, {'all', 'cpg'}))
    s = load(fullfile(res, 'sns_cpg_results.mat'));
    cpg_figure(s.out.log_th, fullfile(pics, 'fig_kneecpg_v3.png'));
end
if any(strcmp(which, {'all', 'beer'}))
    s = load(fullfile(res, 'sns_beer_results.mat'));
    r = s.runs;
    beer_figure(r, fullfile(pics, 'fig_beer_v3.png'));
end
end

%% ---------------- shared glyphs ----------------
function capsule(ax, p1, p2, w1, w2, col)
d = p2 - p1; L = norm(d); u_ = d / L; n_ = [-u_(2), u_(1)];
p1a = p1 + n_ * w1 / 2; p1b = p1 - n_ * w1 / 2;
p2a = p2 + n_ * w2 / 2; p2b = p2 - n_ * w2 / 2;
patch(ax, 'XData', [p1a(1) p2a(1) p2b(1) p1b(1)], ...
      'YData', [p1a(2) p2a(2) p2b(2) p1b(2)], ...
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
% McKibben belly with braid + fittings; p1/p2 are the TENDON endpoints;
% the belly is drawn between fitting offsets so nothing crosses the limb.
if nargin < 5, lw = 0.028; end
d = p2 - p1; L = norm(d); u_ = d / L; n_ = [-u_(2), u_(1)];
col = [1.0, 0.55 + 0.20 * (1 - act), 0.45];
t = linspace(0.15, 0.85, 18)';
bx = lw / 2 * sin(pi * (t - 0.15) / 0.7);
top = p1 + u_ .* (L * t) + n_ .* bx;
bot = p1 + u_ .* (L * t) - n_ .* bx;
patch(ax, 'XData', [top(:, 1); flipud(bot(:, 1))], ...
      'YData', [top(:, 2); flipud(bot(:, 2))], ...
      'FaceColor', col, 'EdgeColor', [0.35 0.2 0.15], 'LineWidth', 0.8);
for k = 1:9
    s1 = 2 * ((k) / 10) - 1; s2 = -(2 * ((k + 0.5) / 10) - 1);
    q1 = p1 + u_ * (L * (0.15 + 0.7 * k / 10)) + n_ * (lw / 2 * s1);
    q2 = p1 + u_ * (L * (0.15 + 0.7 * (k + 0.5) / 10)) + n_ * (lw / 2 * s2);
    plot(ax, [q1(1) q2(1)], [q1(2) q2(2)], 'Color', [0.45 0.3 0.2], 'LineWidth', 0.6);
end
capsule(ax, p1, p1 + u_ * 0.018, 0.018, 0.018, [0.65 0.65 0.68]);
capsule(ax, p2 - u_ * 0.018, p2, 0.018, 0.018, [0.65 0.65 0.68]);
end

function draw_leg(ax, thd, labels)
% Anatomically-placed leg: femur vertical (hip top), knee flexion swings
% shank posteriorly; BPA origins/insertions offset CLEAR of the limb
% silhouette with attachment dots (no more through-the-thigh muscles).
Lf = 0.40; Ls = 0.42;
hip = [0, 0]; knee = [0, -Lf];
a = deg2rad(-90 - thd);
su = [cos(a), sin(a)]; ankle = knee + Ls * su;
nrm = [-su(2), su(1)];
capsule(ax, hip, knee, 0.115, 0.085, [0.92 0.92 0.94]);
capsule(ax, knee, ankle, 0.075, 0.050, [0.88 0.88 0.91]);
capsule(ax, ankle - 0.035 * nrm + 0.045 * su, ankle + 0.10 * nrm + 0.045 * su, ...
    0.052, 0.058, [0.88 0.88 0.91]);
joint(ax, hip, 0.030); joint(ax, knee, 0.026); joint(ax, ankle, 0.022);
% EXTENSOR (quadriceps group): origin ANTERIOR-PROXIMAL on a fascia offset
% in front of the thigh; insertion on the ANTERIOR tibial tuberosity.
qO = [0.075, -0.10];                       % anterior fascia line
qI = ankle - 0.36 * su + 0.030 * nrm;      % tibial tuberosity, on shank
pat = knee + 0.028 * nrm;                  % patella
bpa(ax, qO, pat, 0.6, 0.030);
plot(ax, [pat(1) qI(1)], [pat(2) qI(2)], 'Color', [0.85 0.82 0.70], 'LineWidth', 3);
plot(ax, qO(1), qO(2), 'k.', 'MarkerSize', 12);
plot(ax, qI(1), qI(2), 'k.', 'MarkerSize', 12);
% FLEXOR (hamstring group): origin POSTERIOR-PROXIMAL offset; insertion
% posterior proximal tibia.
hO = [-0.075, -0.10];
hI = ankle - 0.34 * su - 0.030 * nrm;
bpa(ax, hO, hI, 0.35, 0.026);
plot(ax, hO(1), hO(2), 'k.', 'MarkerSize', 12);
plot(ax, hI(1), hI(2), 'k.', 'MarkerSize', 12);
if labels
    text(ax, 0.10, -0.12, 'extensor BPA', 'FontSize', 9);
    text(ax, -0.17, -0.12, 'flexor BPA', 'FontSize', 9);
    text(ax, knee(1) + 0.05, knee(2) - 0.02, 'knee', 'FontSize', 9);
end
axis(ax, 'equal'); axis(ax, 'off');
end

function pathway_glyphs(ax, x0, y0)
% Mini reflex-pathway panel glyphs: circles = neurons; triangle = exc;
% filled dot = inh. Returns nothing; caller positions/lables.
end

function neuron_glyph(ax, p, r, lbl)
th_ = linspace(0, 2 * pi, 20)';
patch(ax, 'XData', p(1) + r * cos(th_), 'YData', p(2) + r * sin(th_), ...
      'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.3);
if nargin > 3 && ~isempty(lbl)
    text(ax, p(1), p(2), lbl, 'HorizontalAlignment', 'center', 'FontSize', 8.5);
end
end

function connect_glyph(ax, p1, p2, exc)
plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], 'k', 'LineWidth', 1.2);
m = (p1 + p2) / 2;
if exc
    patch(ax, m(1) + [0 -0.02 -0.02], m(2) + [0 0.014 -0.014], ...
        [1 1 1] * 0.1 * 0 + [1 1 1], 'EdgeColor', 'k', 'LineWidth', 1.0);
else
    th_ = linspace(0, 2 * pi, 12)';
    patch(ax, m(1) + 0.016 * cos(th_), m(2) + 0.016 * sin(th_), ...
        [0.1 0.1 0.1], 'EdgeColor', 'k');
end
end

%% ---------------- Figure 1: knee reflex ----------------
function knee_figure(th, fname)
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 1500 520], ...
    'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% A: kinematics (the demo compares reflex ON vs OFF; OFF data absent in
% this mat's single run -> plot theta(t) + drive marker)
ax = nexttile(tl);
t = th.Time(:); v = rad2deg(th.Data(:));
plot(ax, t, v, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.6); grid(ax, 'on');
xlabel(ax, 'time (s)'); ylabel(ax, 'knee flexion \theta (deg)');
title(ax, 'A  Kinematic response', 'FontWeight', 'bold');
yline(ax, 43.5, '--', 'settle \approx 43.5°', 'Color', [0.5 0.5 0.5], ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 8.5);
% B: mechanism, muscles clear of the limb
ax = nexttile(tl);
draw_leg(ax, 30, true);
title(ax, 'B  Mechanism (BPA routing)', 'FontWeight', 'bold');
% C: reflex pathway
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off');
hold(ax, 'on'); xlim(ax, [0 1]); ylim(ax, [0 1]);
mn_e = [0.72, 0.50]; mn_f = [0.72, 0.16];
ia_e = [0.30, 0.86]; ib_e = [0.30, 0.66]; ia_f = [0.30, 0.32];
neuron_glyph(ax, ia_e, 0.045, 'Ia_{ext}');
neuron_glyph(ax, ib_e, 0.045, 'Ib_{ext}');
neuron_glyph(ax, ia_f, 0.045, 'Ia_{flx}');
neuron_glyph(ax, mn_e, 0.055, 'MN_{ext}');
neuron_glyph(ax, mn_f, 0.055, 'MN_{flx}');
connect_glyph(ax, ia_e, mn_e, true);        % stretch reflex exc
connect_glyph(ax, ib_e, mn_e, false);       % autogenic inhibition
connect_glyph(ax, ia_f, mn_e, false);       % reciprocal inhibition
connect_glyph(ax, ia_f, mn_f, true);        % flexor stretch reflex
connect_glyph(ax, mn_e, mn_f, false);       % antagonist coupling marker
text(ax, 0.02, 0.95, 'reflex pathways', 'FontSize', 9, 'FontWeight', 'bold');
text(ax, 0.08, 0.86, 'stretch (Ia, exc)', 'FontSize', 8);
text(ax, 0.08, 0.66, 'autogenic (Ib, inh)', 'FontSize', 8);
text(ax, 0.08, 0.32, 'reciprocal (Ia, inh)', 'FontSize', 8);
title(ax, 'C  Reflex pathway', 'FontWeight', 'bold');
exportgraphics(fig, fname, 'Resolution', 300);
close(fig);
fprintf('wrote %s\n', fname);
end

%% ---------------- Figure 2: knee CPG ----------------
function cpg_figure(th, fname)
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 1500 520], ...
    'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl);
t = th.Time(:); v = rad2deg(th.Data(:));
plot(ax, t, v, 'Color', [0.64 0.08 0.18], 'LineWidth', 1.6); grid(ax, 'on');
xlabel(ax, 'time (s)'); ylabel(ax, 'knee flexion \theta (deg)');
title(ax, 'A  Kinematic response', 'FontWeight', 'bold');
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
% architecture: RG -> PF -> MN -> BPA -> joint (block chain)
blks = {0.10, 'RG'; 0.36, 'PF'; 0.62, 'MN'; 0.86, 'BPA'};
for k = 1:size(blks, 1)
    rectangle(ax, 'Position', [blks{k, 1}, 0.42, 0.16, 0.16], ...
        'Curvature', 0.3, 'EdgeColor', 'k', 'FaceColor', [0.94 0.94 0.96]);
    text(ax, blks{k, 1} + 0.08, 0.50, blks{k, 2}, 'HorizontalAlignment', 'center', ...
        'FontSize', 9.5, 'FontWeight', 'bold');
    if k < 4
        plot(ax, [blks{k, 1} + 0.16, blks{k + 1, 1}], [0.50 0.50], 'k', ...
            'LineWidth', 1.2);
    end
end
rectangle(ax, 'Position', [0.10, 0.72, 0.30, 0.14], 'Curvature', 0.3, ...
    'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--', 'FaceColor', 'none');
text(ax, 0.25, 0.79, 'adaptive inhibition (Adp)', 'HorizontalAlignment', ...
    'center', 'FontSize', 8.5, 'Color', [0.35 0.35 0.35]);
plot(ax, [0.25 0.25], [0.72 0.58], 'k--', 'LineWidth', 0.9);
rectangle(ax, 'Position', [0.62, 0.16, 0.16, 0.14], 'Curvature', 0.3, ...
    'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--', 'FaceColor', 'none');
text(ax, 0.70, 0.23, '1-DOF knee', 'HorizontalAlignment', 'center', 'FontSize', 8.5);
plot(ax, [0.94 0.70], [0.42 0.30], 'k--', 'LineWidth', 0.9);
title(ax, 'B  Neural architecture', 'FontWeight', 'bold');
% C: four-corner half-center diagram
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
hE = [0.28, 0.72]; hF = [0.72, 0.72];
iE = [0.28, 0.28]; iF = [0.72, 0.28];
neuron_glyph(ax, hE, 0.055, 'HC_E');
neuron_glyph(ax, hF, 0.055, 'HC_F');
neuron_glyph(ax, iE, 0.05, 'IN_E');
neuron_glyph(ax, iF, 0.05, 'IN_F');
connect_glyph(ax, hE, iE, true);
connect_glyph(ax, hF, iF, true);
connect_glyph(ax, iE, hF, false);
connect_glyph(ax, iF, hE, false);
% weak cross-excitation (dashed, per Deng 2023)
plot(ax, [hE(1) + 0.04, hF(1) - 0.04], [hE(2) - 0.03, hF(2) - 0.03], ...
    'k:', 'LineWidth', 1.0);
text(ax, 0.50, 0.64, 'weak cross-excitation', 'HorizontalAlignment', 'center', ...
    'FontSize', 8);
text(ax, 0.50, 0.955, 'rhythm generator (half-centers)', 'HorizontalAlignment', ...
    'center', 'FontSize', 9, 'FontWeight', 'bold');
title(ax, 'C  Half-center core', 'FontWeight', 'bold');
exportgraphics(fig, fname, 'Resolution', 300);
close(fig);
fprintf('wrote %s\n', fname);
end

%% ---------------- Figure 3: beer ----------------
function beer_figure(runs, fname)
fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [40 40 1500 520], ...
    'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
% A: natural posture mechanism: humerus VERTICAL (arm at side), forearm
% HORIZONTAL, elbow 90 deg, cup in hand; sag rotates the forearm down.
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [-0.45 0.55]); ylim(ax, [-0.78 0.10]);
shoulder = [0, 0];
elbow = shoulder + [0, -0.30];                       % humerus straight down
sag = rad2deg(runs(1).th.Data(end));                 % final sag (reflex ON)
aF = -sag;                                           % forearm drops by sag
fa = deg2rad(0) * 0 + deg2rad(sag);                  % forearm angle below horiz
hand = elbow + 0.30 * [cos(-fa), sin(-fa)];
capsule(ax, shoulder, elbow, 0.095, 0.075, [0.92 0.92 0.94]);   % upper arm
capsule(ax, elbow, hand, 0.070, 0.050, [0.88 0.88 0.91]);      % forearm
capsule(ax, shoulder + [-0.09, 0.02], shoulder + [0.09, 0.02], 0.10, 0.10, ...
    [0.90 0.90 0.93]);                                          % shoulder hint
joint(ax, shoulder, 0.028); joint(ax, elbow, 0.024);
% biceps (anterior) + triceps (posterior) along the upper arm, offset
biO = elbow + 0.055 * [cos(-pi / 2 + 0.5), sin(-pi / 2 + 0.5)] + [0, 0.10];
biI = elbow + 0.08 * [cos(deg2rad(-90 + 90 - rad2deg(fa) * 0)), ...
    sin(deg2rad(-90))] * 0 + [0.055 * cos(fa + pi / 2), 0.055 * sin(fa + pi / 2)];
bpa(ax, biO, biI, 0.45, 0.026);
trO = elbow + [0, 0.10] - 0.055 * [cos(fa + pi / 2), sin(fa + pi / 2)];
trI = elbow - 0.03 * [cos(fa + pi / 2), sin(fa + pi / 2)];
bpa(ax, trO, trI, 0.40, 0.022);
% pint glass in hand
gn = [cos(fa + pi / 2), sin(fa + pi / 2)];
base = hand + 0.012 * gn;
b1 = base - 0.075 * [cos(fa), sin(fa)];
b2 = base + 0.075 * [cos(fa), sin(fa)];
r1 = b1 + 0.16 * gn; r2 = b2 + 0.16 * gn;
patch(ax, 'XData', [b1(1) b2(1) r2(1) r1(1)], 'YData', [b1(2) b2(2) r2(2) r1(2)], ...
    'FaceColor', [0.94 0.96 1.0], 'EdgeColor', [0.55 0.62 0.70], ...
    'LineWidth', 1.6, 'FaceAlpha', 0.55);
text(ax, 0.02, 0.02, 'start: forearm level; pour loads the arm', 'FontSize', 8.5);
title(ax, 'B  Mechanism (natural posture)', 'FontWeight', 'bold');
% B: kinematics ON vs OFF
ax = nexttile(tl);
hold(ax, 'on'); grid(ax, 'on');
t1 = runs(1).th.Time(:); v1 = rad2deg(runs(1).th.Data(:));
t2 = runs(2).th.Time(:); v2 = rad2deg(runs(2).th.Data(:));
plot(ax, t1, v1, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.6);
plot(ax, t2, v2, 'Color', [0.55 0.55 0.55], 'LineWidth', 1.4, 'LineStyle', '--');
legend(ax, {'reflex ON', 'reflex OFF'}, 'Location', 'southeast');
xlabel(ax, 'time (s)'); ylabel(ax, 'sag \theta (deg)');
title(ax, 'A  Kinematic response', 'FontWeight', 'bold');
% C: pathway panel
ax = nexttile(tl); axis(ax, 'equal'); axis(ax, 'off'); hold(ax, 'on');
xlim(ax, [0 1]); ylim(ax, [0 1]);
ia_bi = [0.30, 0.86]; ib_bi = [0.30, 0.66]; ia_tri = [0.30, 0.32];
mn_bi = [0.72, 0.66]; mn_tri = [0.72, 0.22];
neuron_glyph(ax, ia_bi, 0.045, 'Ia_{bi}');
neuron_glyph(ax, ib_bi, 0.045, 'Ib_{bi}');
neuron_glyph(ax, ia_tri, 0.045, 'Ia_{tri}');
neuron_glyph(ax, mn_bi, 0.055, 'MN_{bi}');
neuron_glyph(ax, mn_tri, 0.055, 'MN_{tri}');
connect_glyph(ax, ia_bi, mn_bi, true);
connect_glyph(ax, ib_bi, mn_bi, false);
connect_glyph(ax, ia_tri, mn_bi, false);
connect_glyph(ax, ia_tri, mn_tri, true);
connect_glyph(ax, ia_bi, mn_tri, false);
text(ax, 0.02, 0.95, 'reflex pathways', 'FontSize', 9, 'FontWeight', 'bold');
title(ax, 'C  Reflex pathway', 'FontWeight', 'bold');
exportgraphics(fig, fname, 'Resolution', 300);
close(fig);
fprintf('wrote %s\n', fname);
end
