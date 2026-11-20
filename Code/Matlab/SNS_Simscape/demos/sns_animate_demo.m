%% sns_animate_demo.m — animate the mechanical parts of the SNS demos
%
% The demos are reduced-order (1-DOF), so there is no Simscape Multibody
% geometry yet. This script replays the LOGGED joint motion in a 2D
% mechanism drawing (journal-style heavy lines) and saves a GIF per demo.
%
% KNEE GEOMETRY (fixed 2026-09-22 — the old drawing swung the shank UP
% alongside the femur, so the leg read as folded 180 deg at theta~0):
%   sagittal view, anterior = +x (right). Femur vertical, hip at top, knee
%   at bottom. theta = knee FLEXION (0 = straight leg); flexion swings the
%   shank POSTERIORLY, shank direction angle = -(90 + theta) deg.
%   Quadriceps: belly from the femur (front, upper) to the patella (fixed
%   just anterior of the knee) + patellar tendon to the tibial tuberosity
%   (body-fixed on the shank, rotates with theta) — the tendon visibly
%   lengthens/realigns as the knee flexes. Hamstring: belly from the femur
%   (back, upper) to a body-fixed insertion on the posterior proximal
%   tibia — shortens with flexion.
%
% ELBOW (beer) GEOMETRY (fixed 2026-09-22): the biceps ORIGIN is on the
% upper arm at the shoulder (it was drawn floating 10 cm lateral to the
% shoulder), inserting ~20% along the forearm; the TRICEPS runs posterior
% along the upper arm to the olecranon at the elbow. Cup fill stays LEVEL
% (gravity) while the arm sags — the whole point of the demo.
%
% Usage:  sns_animate_demo            % animates all
%         sns_animate_demo('knee')    % KneeReflexDemo + BPACPGLegDemo runs
%         sns_animate_demo('beer')

function sns_animate_demo(which)
if nargin < 1, which = 'all'; end
here = fileparts(mfilename('fullpath'));
cd(here);
res = fullfile(here, '..', 'results');
anim = fullfile(res, 'animations');
if ~exist(anim, 'dir'), mkdir(anim); end

if any(strcmp(which, {'all', 'knee'}))
    f = fullfile(res, 'sns_demo_results.mat');
    if exist(f, 'file')
        s = load(f);
        th = s.out.log_th;
        animate_knee(th.Time, th.Data*180/pi, ...
            'SNS knee reflex — Ia/Ib reflexes + antagonist BPAs', ...
            fullfile(anim, 'sns_knee_animation.gif'));
    else
        fprintf('knee: sns_demo_results.mat not found - run sns_run_demo first\n');
    end
    f = fullfile(res, 'sns_cpg_results.mat');
    if exist(f, 'file')
        s = load(f);
        th = s.out.log_th;
        animate_knee(th.Time, th.Data*180/pi, ...
            'SNS CPG knee — half-center RG driving antagonist BPAs', ...
            fullfile(anim, 'sns_cpg_knee_animation.gif'));
    else
        fprintf('cpg: sns_cpg_results.mat not found - run sns_run_cpg_demo first\n');
    end
end
if any(strcmp(which, {'all', 'beer'}))
    f = fullfile(res, 'sns_beer_results.mat');
    if exist(f, 'file')
        s = load(f);
        r = s.runs;   % struct array: (1) reflex ON, (2) OFF
        animate_beer(r(1).th.Time, r(1).th.Data*180/pi, r(1).mCup.Data, ...
            fullfile(anim, 'sns_beer_animation.gif'));
    else
        fprintf('beer: sns_beer_results.mat not found - run sns_run_beer_demo first\n');
    end
end
end

%% ---------------- knee mechanism ----------------
function animate_knee(t, th_deg, ttl, gifName)
% Geometry (m): femur 0.40 vertical, hip at top; shank 0.42 from knee.
% theta = knee flexion (deg, 0 = full extension). Flexion swings the shank
% POSTERIORLY (foot back): shank direction = -(90 + theta) deg from +x.
Lf = 0.40; Ls = 0.42;
kneeX = 0; kneeY = -Lf;
shankA = @(deg) deg2rad(-90 - deg);        % direction from knee toward ankle
u = @(a) [cos(a), sin(a)];                 % shank axis unit vector
v = @(a) [-sin(a), cos(a)];                % anterior side of the shank

fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 640 640], 'Color', 'w');
ax = axes(fig); hold(ax, 'on'); axis(ax, 'equal');
xlim(ax, [-0.55 0.30]); ylim(ax, [-0.90 0.12]);
axis(ax, 'off');

% muscle anchors
patella = [0.033, kneeY + 0.005];              % fixed, just anterior of the knee
quadOrigin  = [0.045, kneeY + 0.27];           % femur, anterior, upper
hamOrigin   = [-0.045, kneeY + 0.27];          % femur, posterior, upper

n = numel(t); skip = max(1, floor(n/150));   % ~150 frames
first = true;
for k = 1:skip:n
    cla(ax); hold(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-0.55 0.30]); ylim(ax, [-0.90 0.12]); axis(ax, 'off');
    title(ax, sprintf('%s — \\theta = %+.1f deg  (t = %.2f s)', ttl, th_deg(k), t(k)), ...
        'FontWeight', 'bold', 'FontSize', 10);

    a = shankA(th_deg(k));
    ua = u(a); va = v(a);                      % shank axis / anterior side
    ank = [kneeX, kneeY] + Ls*ua;
    % femur + hip block
    rectangle(ax, 'Position', [-0.06 0.02 0.12 0.06], 'FaceColor', [0.8 0.8 0.8], ...
        'EdgeColor', 'k', 'LineWidth', 1.8);
    plot(ax, [0 0], [0.02 kneeY], 'k-', 'LineWidth', 3.5);
    % shank + foot (foot points to the anterior side of the shank)
    plot(ax, [kneeX ank(1)], [kneeY ank(2)], 'k-', 'LineWidth', 3.5);
    toe = ank + 0.10*va;
    plot(ax, [ank(1) toe(1)], [ank(2) toe(2)], 'k-', 'LineWidth', 3);
    plot(ax, [ank(1) ank(1) - 0.03*ua(1)], [ank(2) ank(2) - 0.03*ua(2)], 'k-', 'LineWidth', 3);
    % knee pivot
    plot(ax, kneeX, kneeY, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 1.8);

    % quadriceps: belly (femur -> patella) + patellar tendon (patella -> tuberosity)
    tuber = [kneeX, kneeY] + 0.075*ua + 0.030*va;   % tibial tuberosity, rotates
    draw_muscle(ax, quadOrigin, patella, 'quad', [0.94 0.76 0.74]);
    plot(ax, [patella(1) tuber(1)], [patella(2) tuber(2)], 'k-', 'LineWidth', 2);
    % hamstring: belly (femur -> posterior proximal tibia, rotates)
    hamIns = [kneeX, kneeY] + 0.075*ua - 0.030*va;
    draw_muscle(ax, hamOrigin, hamIns, 'hamstrings', [0.94 0.76 0.74]);

    drawnow;
    frame = getframe(fig);
    if first
        [A, map] = rgb2ind(frame2im(frame), 256);
        imwrite(A, map, gifName, 'gif', 'LoopCount', 0, 'DelayTime', 0.05);
        first = false;
    else
        imwrite(rgb2ind(frame2im(frame), 256), gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end
close(fig);
fprintf('wrote %s (%d frames)\n', gifName, numel(skip:skip:n));
end

%% ---------------- elbow + cup mechanism ----------------
function animate_beer(t, th_deg, mCup, gifName)
% Geometry (m): shoulder at origin, upper arm 0.30 vertical down, forearm
% 0.28 horizontal at level; theta = sag angle (deg, + = cup tips down).
% BICEPS: origin on the upper arm just distal of the shoulder (anterior),
% insertion ~20% along the forearm (radial tuberosity). TRICEPS: origin on
% the posterior upper arm, insertion at the olecranon (elbow).
Lua = 0.30; Lfa = 0.28; Hcup = 0.10; Wcup = 0.075;

fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 640 480], 'Color', 'w');
ax = axes(fig); hold(ax, 'on'); axis(ax, 'equal');
xlim(ax, [-0.20 0.45]); ylim(ax, [-0.45 0.10]);
axis(ax, 'off');

biOrigin = [0.018, -0.05];     % anterior upper arm, just below the shoulder
triOrigin = [-0.018, -0.05];   % posterior upper arm

n = numel(t); skip = max(1, floor(n/150));
first = true;
cmapBeer = [1 1 1; 0 0 0; 0.85 0.45 0.15; 0.94 0.76 0.74; 0.8 0.8 0.8; 0.5 0.75 0.95];
for k = 1:skip:n
    cla(ax); hold(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-0.20 0.45]); ylim(ax, [-0.45 0.10]); axis(ax, 'off');
    title(ax, sprintf('beer = %.0f g   \\theta = %+.1f deg  (t = %.1f s)', ...
        1000*mCup(k), th_deg(k), t(k)), 'FontWeight', 'bold');

    th = deg2rad(th_deg(k));
    sh = [0, 0]; el = [0, -Lua];
    fa = [cos(-th), sin(-th)];               % forearm direction (+x at level)
    faUp = [-fa(2), fa(1)];                  % forearm "up" side
    hn = el + Lfa*fa;                        % hand
    % triceps FIRST (behind the bones): posterior upper arm -> olecranon
    triIns = el + 0.015*fa - 0.022*faUp;
    draw_muscle(ax, triOrigin, triIns, 'triceps', [0.90 0.66 0.64], 0.028);
    % bones + joints
    plot(ax, [sh(1) el(1)], [sh(2) el(2)], 'k-', 'LineWidth', 4);
    plot(ax, [el(1) hn(1)], [el(2) hn(2)], 'k-', 'LineWidth', 4);
    plot(ax, sh(1), sh(2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.8);
    plot(ax, el(1), el(2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.8);
    % hand stub gripping the cup
    plot(ax, [hn(1) hn(1)+0.045*fa(1)], [hn(2) hn(2)+0.045*fa(2)], 'k-', 'LineWidth', 4);
    % biceps AFTER the bones (anterior): shoulder origin, crosses the elbow,
    % inserts on the proximal forearm (lengthens with sag)
    biIns = el + 0.17*Lfa*fa + 0.024*faUp;
    draw_muscle(ax, biOrigin, biIns, 'biceps', [0.94 0.76 0.74], 0.036);

    % cup below the hand, axis along the forearm; beer surface stays HORIZONTAL
    cupDir = [sin(-th), -cos(-th)];      % down along cup axis
    side = fa;
    c1 = hn + side*(Wcup/2);
    c2 = hn - side*(Wcup/2);
    b1 = c1 + cupDir*Hcup; b2 = c2 + cupDir*Hcup;
    patch(ax, [c1(1) c2(1) b2(1) b1(1)], [c1(2) c2(2) b2(2) b1(2)], [1 1 1], 'EdgeColor', 'k', 'LineWidth', 1.8);
    fillFrac = min(mCup(k)/0.62, 1);
    if fillFrac > 0.02
        hBeer = fillFrac*Hcup;
        s1 = c1 + cupDir*(Hcup - hBeer);
        s2 = c2 + cupDir*(Hcup - hBeer);
        q1 = s1 + cupDir*hBeer; q2 = s2 + cupDir*hBeer;
        patch(ax, [s1(1) s2(1) q2(1) q1(1)], [s1(2) s2(2) q2(2) q1(2)], ...
            [0.85 0.45 0.15], 'EdgeColor', 'none');
    end
    plot(ax, [-0.18 0.43], [hn(2) hn(2)], 'k--', 'LineWidth', 0.8);

    drawnow;
    frame = getframe(fig);
    if first
        A = rgb2ind(frame2im(frame), cmapBeer);
        imwrite(A, cmapBeer, gifName, 'gif', 'LoopCount', 0, 'DelayTime', 0.06);
        first = false;
    else
        imwrite(rgb2ind(frame2im(frame), cmapBeer), gifName, 'gif', 'WriteMode', 'append', 'DelayTime', 0.06);
    end
end
close(fig);
fprintf('wrote %s (%d frames)\n', gifName, numel(skip:skip:n));
end

%% ---------------- fusiform muscle between two anchor points ----------------
function draw_muscle(ax, p1, p2, lbl, fill, wmax)
if nargin < 5, fill = [0.94 0.76 0.74]; end
if nargin < 6, wmax = 0.032; end
d = p2 - p1; L = norm(d);
if L < 1e-6, return; end
ang = atan2d(d(2), d(1));
w = min(wmax, L*0.35);
c = (p1 + p2)/2;
th = linspace(0, 2*pi, 49);
locx = (L/2)*cos(th); locy = w*sin(th);
wx = c(1) + locx*cosd(ang) - locy*sind(ang);
wy = c(2) + locx*sind(ang) + locy*cosd(ang);
patch(ax, wx, wy, fill, 'EdgeColor', 'k', 'LineWidth', 1.8);
text(ax, c(1), c(2), lbl, 'HorizontalAlignment', 'center', 'FontSize', 7);
end
