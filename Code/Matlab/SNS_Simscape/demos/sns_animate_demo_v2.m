function sns_animate_demo_v2(which)
% V2 mechanism animations (Ben, 2026-10-24: "make it look like an actual
% leg and an actual cup - the image in the dissertation is terrible").
% Replays the SAME logged joint motion as sns_animate_demo but draws
% anatomically shaped limbs: tapered limb capsules with joint markers,
% braided McKibben BPA bellies (crosshatch) + tendons, a foot; the beer
% arm holds a pint glass whose beer SURFACE stays level while the arm
% sags. GIFs to results\animations\v2\, stills to results\pictures\v2\.

if nargin < 1, which = 'all'; end
here = fileparts(mfilename('fullpath'));
cd(here);
res = fullfile(here, '..', 'results');
anim = fullfile(res, 'animations', 'v2');
pics = fullfile(res, 'pictures', 'v2');
if ~exist(anim, 'dir'), mkdir(anim); end
if ~exist(pics, 'dir'), mkdir(pics); end

if any(strcmp(which, {'all', 'knee'}))
    f = fullfile(res, 'sns_demo_results.mat');
    if exist(f, 'file')
        s = load(f); th = s.out.log_th;
        animate_knee(th.Time, th.Data * 180 / pi, ...
            'SNS knee reflex', fullfile(anim, 'knee_reflex_v2.gif'), ...
            fullfile(pics, 'knee_reflex_v2_still.png'));
    end
    f = fullfile(res, 'sns_cpg_results.mat');
    if exist(f, 'file')
        s = load(f); th = s.out.log_th;
        animate_knee(th.Time, th.Data * 180 / pi, ...
            'SNS CPG knee', fullfile(anim, 'knee_cpg_v2.gif'), ...
            fullfile(pics, 'knee_cpg_v2_still.png'));
    end
end
if any(strcmp(which, {'all', 'beer'}))
    f = fullfile(res, 'sns_beer_results.mat');
    if exist(f, 'file')
        s = load(f); r = s.runs;
        animate_beer(r(1).th.Time, r(1).th.Data * 180 / pi, ...
            r(1).mCup.Data, 'reflex ON', ...
            fullfile(anim, 'beer_reflex_on_v2.gif'), ...
            fullfile(pics, 'beer_reflex_on_v2_still.png'));
        animate_beer(r(2).th.Time, r(2).th.Data * 180 / pi, ...
            r(2).mCup.Data, 'reflex OFF', ...
            fullfile(anim, 'beer_reflex_off_v2.gif'), ...
            fullfile(pics, 'beer_reflex_off_v2_still.png'));
    end
end
end

%% ------------- shared drawing helpers -------------
function capsule(ax, p1, p2, w1, w2, col)
% Tapered limb segment: quad strip perpendicular to the axis.
d = p2 - p1; L = norm(d); u_ = d / L; n_ = [-u_(2), u_(1)];
p1a = p1 + n_ * w1 / 2; p1b = p1 - n_ * w1 / 2;
p2a = p2 + n_ * w2 / 2; p2b = p2 - n_ * w2 / 2;
patch(ax, 'XData', [p1a(1) p2a(1) p2b(1) p1b(1)], ...
          'YData', [p1a(2) p2a(2) p2b(2) p1b(2)], ...
          'FaceColor', col, 'EdgeColor', [0.25 0.25 0.25], 'LineWidth', 1.0);
% rounded ends
for q = 1:2
    pts = {p1, p2}; ws = {w1, w2};
    th_ = linspace(0, 2 * pi, 24)';
    patch(ax, 'XData', pts{q}(1) + ws{q} / 2 * cos(th_), ...
              'YData', pts{q}(2) + ws{q} / 2 * sin(th_), ...
              'FaceColor', col, 'EdgeColor', 'none');
end
end

function joint(ax, p, r, col)
th_ = linspace(0, 2 * pi, 24)';
patch(ax, 'XData', p(1) + r * cos(th_), 'YData', p(2) + r * sin(th_), ...
      'FaceColor', col, 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
end

function bpa(ax, p1, p2, act, lw)
% McKibben BPA belly: braided barrel + end fittings; act in [0,1] colors.
if nargin < 5, lw = 0.030; end
d = p2 - p1; L = norm(d); u_ = d / L; n_ = [-u_(2), u_(1)];
col = [1.00, 0.55 + 0.20 * (1 - act), 0.45];      % salmon -> red with act
th_ = linspace(0, pi, 20)';
bx = lw / 2 * sin(th_);                            % barrel profile
top = p1 + u_ .* (L * th_ / pi) + n_ .* bx;
bot = p1 + u_ .* (L * th_ / pi) - n_ .* bx;
patch(ax, 'XData', [top(:, 1); flipud(bot(:, 1))], ...
          'YData', [top(:, 2); flipud(bot(:, 2))], ...
          'FaceColor', col, 'EdgeColor', [0.35 0.2 0.15], 'LineWidth', 0.8);
% braid: zigzag between top and bottom edges
m = 12;
for k = 1:m
    t1 = k / m; t2 = (k + 0.5) / m;
    s1 = 2 * (k / m) - 1; s2 = -(2 * ((k + 0.5) / m) - 1);   % alternate sides
    q1 = p1 + u_ * (L * t1) + n_ * (lw / 2 * s1);
    q2 = p1 + u_ * (L * t2) + n_ * (lw / 2 * s2);
    plot(ax, [q1(1) q2(1)], [q1(2) q2(2)], 'Color', [0.45 0.3 0.2], 'LineWidth', 0.6);
end
% clamp fittings
capsule(ax, p1 - u_ * 0.012, p1 + u_ * 0.012, 0.020, 0.020, [0.65 0.65 0.68]);
capsule(ax, p2 - u_ * 0.012, p2 + u_ * 0.012, 0.020, 0.020, [0.65 0.65 0.68]);
end

function tendon(ax, p1, p2)
plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], 'Color', [0.85 0.82 0.70], ...
     'LineWidth', 3.5);
plot(ax, [p1(1) p2(1)], [p1(2) p2(2)], 'Color', [0.55 0.5 0.35], ...
     'LineWidth', 1.0, 'LineStyle', ':');
end

function frames = grab(fig, n)
f = getframe(fig);
if n == 1
    frames = im2frame(f.cdata, f.colormap); %#ok<INUSD>
end
end

%% ------------- knee mechanism (v2) -------------
function animate_knee(t, thd_, ttl, gifName, stillName)
Lf = 0.40; Ls = 0.42;
hip = [0, 0]; knee = [0, -Lf];
skin = [0.92 0.92 0.94]; skin2 = [0.88 0.88 0.91];
fig = figure('Visible', 'off', 'Units', 'pixels', ...
             'Position', [100 100 560 660], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');
n = numel(t); skip = max(1, floor(n / 110));
idx = [1:skip:n, n];
F = {}; MAP = []; frameMaps = {};
for k = idx
    cla(ax); hold(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-0.42 0.30]); ylim(ax, [-0.95 0.10]); axis(ax, 'off');
    % kinematics: theta flexion swings shank posteriorly
    a = deg2rad(-90 - thd_(k));
    shankU = [cos(a), sin(a)];
    ankle = knee + Ls * shankU;
    nrm = [-shankU(2), shankU(1)];                  % anterior side
    footToe = ankle + 0.10 * nrm + 0.045 * shankU;  % toes forward-down
    footHeel = ankle - 0.035 * nrm + 0.045 * shankU;
    % limb segments
    capsule(ax, hip, knee, 0.115, 0.085, skin);            % thigh
    capsule(ax, knee, ankle, 0.075, 0.050, skin2);         % shank
    capsule(ax, footHeel, footToe, 0.052, 0.058, skin2);   % foot
    % knee indicates flexion in the title; draw patella + joints
    patella = knee + 0.030 * nrm;
    joint(ax, hip, 0.030, [0.75 0.75 0.78]);
    joint(ax, knee, 0.026, [0.80 0.80 0.84]);
    joint(ax, ankle, 0.022, [0.75 0.75 0.78]);
    % actuators: routed OUTSIDE the limb silhouette (v3 anatomy fix):
    % extensor origin on an anterior fascia offset, flexor posterior.
    quadO = [0.075, -0.10];
    flex = thd_(k) / 90;                                     % 0..~0.55
    bpa(ax, quadO, patella, min(1, max(0, 0.35 + 0.6 * (1 - flex))), 0.034);
    tibTub = ankle - 0.36 * shankU + 0.030 * nrm;            % on shank
    tendon(ax, patella, tibTub);
    % hamstring BPA (post. thigh -> post. proximal tibia)
    hamO = [-0.075, -0.10];
    hamI = ankle - 0.34 * shankU - 0.024 * nrm;
    bpa(ax, hamO, hamI, min(1, max(0, 0.25 + 0.7 * flex)), 0.030);
    title(ax, sprintf('%s — knee \\theta = %.1f° (t = %.2f s)', ttl, thd_(k), t(k)), ...
        'FontWeight', 'bold', 'FontSize', 10.5);
    f = getframe(fig);
    [A, MAP] = rgb2ind(f.cdata, 256); %#ok<ASGLU>
    F{end + 1} = A; frameMaps{end + 1} = MAP; %#ok<AGROW>
    if k == idx(round(numel(idx) * 0.45))
        exportgraphics(fig, stillName, 'Resolution', 220);
    end
end
imwrite(F{1}, frameMaps{1}, gifName, 'WriteMode', 'overwrite', 'DelayTime', 0.03);
for k = 2:numel(F)
    imwrite(F{k}, frameMaps{k}, gifName, 'WriteMode', 'append', 'DelayTime', 0.03);
end
fprintf('wrote %s (+ still)\n', gifName);
close(fig);
end

%% ------------- beer arm mechanism (v2) -------------
function animate_beer(t, thd_, mcup, mode, gifName, stillName)
% Forearm holds a pint glass; arm sags by theta (deg, + = sag down).
Lua = 0.30; Lfa = 0.30;
shoulder = [0, 0];
elbow = shoulder + Lua * [cos(deg2rad(-20 - thd_(1) * 0)), sin(deg2rad(-20))];
fig = figure('Visible', 'off', 'Units', 'pixels', ...
             'Position', [100 100 620 460], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');
n = numel(t); skip = max(1, floor(n / 110));
idx = [1:skip:n, n];
F = {}; frameMaps = {};
for k = idx
    cla(ax); hold(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-0.35 0.75]); ylim(ax, [-0.55 0.28]); axis(ax, 'off');
    % arm in NATURAL posture (v3 fix, Ben 2026-10-26): humerus VERTICAL
    % (arm at side, elbow ~10 deg abduction), forearm horizontal, cup in
    % hand; sag rotates the forearm down from level.
    aU = deg2rad(-80);
    aF = deg2rad(0 - thd_(k));                     % sag = forearm drops
    elbow = shoulder + Lua * [cos(aU), sin(aU)];
    hand = elbow + Lfa * [cos(aF), sin(aF)];
    skin = [0.92 0.92 0.94];
    capsule(ax, shoulder, elbow, 0.095, 0.075, skin);
    capsule(ax, elbow, hand, 0.070, 0.050, skin);
    joint(ax, shoulder, 0.028, [0.75 0.75 0.78]);
    joint(ax, elbow, 0.024, [0.80 0.80 0.84]);
    % muscles: biceps BPA (ant. upper arm -> forearm), triceps post.
    biO = shoulder + 0.05 * [cos(aU + pi / 2), sin(aU + pi / 2)] + ...
          0.10 * [cos(aU), sin(aU)];
    biI = elbow + 0.085 * [cos(aF), sin(aF)];
    sag = max(0, thd_(k)) / 8;
    bpa(ax, biO, biI, min(1, 0.35 + 0.6 * sag), 0.030);
    trO = shoulder - 0.05 * [cos(aU + pi / 2), sin(aU + pi / 2)] + ...
          0.08 * [cos(aU), sin(aU)];
    trI = elbow - 0.03 * [cos(aF + pi / 2), sin(aF + pi / 2)];
    bpa(ax, trO, trI, max(0, 0.40 - 0.5 * sag), 0.024);
    % pint glass in the hand, glass tilted with the forearm
    gl = 0.075;                                   % glass half-width
    gh = 0.16;                                    % glass height
    gdir = [cos(aF), sin(aF)]; gn = [-gdir(2), gdir(1)];
    base = hand + 0.012 * gn;
    rim = base + gh * gn; rim2 = base + gh * gn - 0.010 * gdir;  % slight taper
    b1 = base - gl * gdir; b2 = base + gl * gdir;
    r1 = rim - 0.92 * gl * gdir; r2 = rim + 0.92 * gl * gdir;
    patch(ax, 'XData', [b1(1) b2(1) r2(1) r1(1)], ...
          'YData', [b1(2) b2(2) r2(2) r1(2)], ...
          'FaceColor', [0.94 0.96 1.00], 'EdgeColor', [0.55 0.62 0.70], ...
          'LineWidth', 1.6, 'FaceAlpha', 0.55);
    % beer fill: level stays HORIZONTAL in the tilted glass
    frac = min(1, mcup(k) / 0.5);
    hF = (0.12 + 0.62 * frac) * gh;               % fill height along glass axis
    yFill = base(2) + hF;
    % intersect walls with y = yFill
    wL = wall_pt(b1, r1, yFill); wR = wall_pt(b2, r2, yFill);
    if ~any(isnan(wL)) && ~any(isnan(wR))
        % beer body polygon: surface line down to the (tilted) base
        patch(ax, 'XData', [wL(1) wR(1) b2(1) b1(1)], ...
              'YData', [wL(2) wR(2) b2(2) b1(2)], ...
              'FaceColor', [0.96 0.74 0.16], 'EdgeColor', 'none');
        % foam
        plot(ax, [wL(1) wR(1)], [wL(2) wR(2)], 'Color', [1 1 1], 'LineWidth', 4);
    end
    title(ax, sprintf('Beer-cup reflex (%s) — sag \\theta = %.2f°, m = %.2f kg (t = %.1f s)', ...
        mode, thd_(k), mcup(k), t(k)), 'FontWeight', 'bold', 'FontSize', 10.5);
    f = getframe(fig);
    [A, MAP] = rgb2ind(f.cdata, 256); %#ok<ASGLU>
    F{end + 1} = A; frameMaps{end + 1} = MAP; %#ok<AGROW>
    if k == idx(round(numel(idx) * 0.75))
        exportgraphics(fig, stillName, 'Resolution', 220);
    end
end
imwrite(F{1}, frameMaps{1}, gifName, 'WriteMode', 'overwrite', 'DelayTime', 0.03);
for k = 2:numel(F)
    imwrite(F{k}, frameMaps{k}, gifName, 'WriteMode', 'append', 'DelayTime', 0.03);
end
fprintf('wrote %s (+ still)\n', gifName);
close(fig);
end

function p = wall_pt(a, b_, y)
% point on segment a->b at height y (NaN if outside)
t = (y - a(2)) / (b_(2) - a(2));
if t < 0 || t > 1, p = [NaN NaN]; else, p = a + t * (b_ - a); end
end
