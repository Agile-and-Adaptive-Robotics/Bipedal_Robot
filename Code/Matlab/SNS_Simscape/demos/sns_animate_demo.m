%% sns_animate_demo.m — animate the mechanical parts of the SNS demos
%
% The demos are reduced-order (1-DOF), so there is no Simscape Multibody
% geometry yet. This script replays the LOGGED joint motion in a 2D
% mechanism drawing (journal-style heavy lines) and saves a GIF per demo:
%
%   knee (KneeReflexDemo / BPACPGLegDemo): femur fixed vertical, shank
%       swings by theta, BPA muscles drawn as fusiform actuators on each
%       side, beer-free.
%   elbow (BeerCupReflexDemo): upper arm vertical, forearm rotates by the
%       sag angle theta, cup in hand with a beer fill that stays LEVEL
%       (gravity) while the cup tilts — the whole point of the demo.
%
% Usage:  sns_animate_demo            % animates both, writes GIFs
%         sns_animate_demo('beer')    % just the beer demo
%         sns_animate_demo('knee')
%
% For real 3D mechanical parts: open a Simscape Multibody model (the
% imported OpenSim/SolidWorks models, see sns_osim_import.m) — Mechanics
% Explorer opens automatically on simulation. Blocked on this machine's
% license until SimMechanics is available.

function sns_animate_demo(which)
if nargin < 1, which = 'both'; end
here = fileparts(mfilename('fullpath'));
cd(here);

if any(strcmp(which, {'both', 'knee'}))
    f = fullfile(here, '..', 'results', 'sns_cpg_results.mat');
    if exist(f, 'file')
        s = load(f);
        th = s.out.log_th;
        animate_knee(th.Time, th.Data*180/pi, fullfile(here, '..', 'results', 'animations', 'sns_knee_animation.gif'));
    else
        fprintf('knee: sns_cpg_results.mat not found - run sns_run_cpg_demo first\n');
    end
end
if any(strcmp(which, {'both', 'beer'}))
    f = fullfile(here, '..', 'results', 'sns_beer_results.mat');
    if exist(f, 'file')
        s = load(f);
        r = s.runs;   % struct array: (1) reflex ON, (2) OFF
        animate_beer(r(1).th.Time, r(1).th.Data*180/pi, r(1).mCup.Data, ...
            fullfile(here, '..', 'results', 'animations', 'sns_beer_animation.gif'));
    else
        fprintf('beer: sns_beer_results.mat not found - run sns_run_beer_demo first\n');
    end
end
end

%% ---------------- knee mechanism ----------------
function animate_knee(t, th_deg, gifName)
% geometry (m): femur 0.40 vertical; shank 0.42 from knee; theta measured
% from full extension. BPA muscles drawn between femur and shank anchors.
Lf = 0.40; Ls = 0.42; rk = 0.05;      % knee pivot offset from femur axis
kneeX = 0; kneeY = -Lf;
shankA = @(deg) deg2rad(90 + deg);    % shank direction from knee

fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 640 640], 'Color', 'w');
ax = axes(fig); hold(ax, 'on'); axis(ax, 'equal');
xlim(ax, [-0.45 0.45]); ylim(ax, [-0.75 0.15]);
axis(ax, 'off');
title(ax, 'SNS CPG knee demo — BPA antagonist pair on 1-DOF knee', 'FontWeight', 'bold');

% static: femur + hip
rectangle(ax, 'Position', [-0.06 0.05 0.12 0.05], 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'k', 'LineWidth', 1.8);
plot(ax, [0 0], [0.05 kneeY], 'k-', 'LineWidth', 3);

% BPA actuator anchors (femur side, knee side)
extF = [0.05, -0.06]; extK = [0.035, kneeY + 0.10];    % extensor (front)
flxF = [-0.05, -0.06]; flxK = [-0.035, kneeY + 0.10];  % flexor (back)

n = numel(t); skip = max(1, floor(n/150));   % ~150 frames
first = true;
for k = 1:skip:n
    cla(ax); hold(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-0.45 0.45]); ylim(ax, [-0.75 0.15]); axis(ax, 'off');
    title(ax, sprintf('SNS CPG knee — \\theta = %+.1f deg  (t = %.2f s)', th_deg(k), t(k)), 'FontWeight', 'bold');

    % shank + foot
    a = shankA(th_deg(k));
    ank = [kneeX + Ls*cos(a), kneeY + Ls*sin(a)];
    plot(ax, [kneeX ank(1)], [kneeY ank(2)], 'k-', 'LineWidth', 3);
    plot(ax, [ank(1)-0.06 ank(1)+0.02], [ank(2)-0.02 ank(2)-0.02], 'k-', 'LineWidth', 3);
    % knee pivot
    plot(ax, kneeX, kneeY, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.8);

    % muscles: fusiform ellipses along the anchor lines (red = actuator)
    draw_muscle(ax, extF, extK, 'ext');
    draw_muscle(ax, flxF, flxK, 'flex');
    % muscle anchor lines to the shank (move with theta)
    sh1 = [kneeX + 0.10*cos(a - 0.35), kneeY + 0.10*sin(a - 0.35)];
    sh2 = [kneeX - 0.10*cos(a - 0.35), kneeY - 0.10*sin(a - 0.35) + 0.0];
    plot(ax, [extK(1) sh1(1)], [extK(2) sh1(2)], 'k-', 'LineWidth', 1.8);
    plot(ax, [flxK(1) sh2(1)], [flxK(2) sh2(2)], 'k-', 'LineWidth', 1.8);

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
% geometry (m): shoulder at origin, upper arm 0.30 vertical down, forearm
% 0.28 horizontal at level; theta = sag angle (deg, + = cup tips down).
% Beer surface stays horizontal; fill height from mCup.
Lua = 0.30; Lfa = 0.28; Lcup = 0.085; Wcup = 0.075; Hcup = 0.10;

fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 640 480], 'Color', 'w');
ax = axes(fig); hold(ax, 'on'); axis(ax, 'equal');
xlim(ax, [-0.12 0.45]); ylim(ax, [-0.42 0.10]);
axis(ax, 'off');
title(ax, 'Beer-cup reflex demo — Ia/Ib reflexes hold the cup level while beer pours', 'FontWeight', 'bold');

n = numel(t); skip = max(1, floor(n/150));
first = true;
cmapBeer = [1 1 1; 0 0 0; 0.85 0.45 0.15; 0.94 0.76 0.74; 0.8 0.8 0.8; 0.5 0.75 0.95];
for k = 1:skip:n
    cla(ax); hold(ax, 'on'); axis(ax, 'equal');
    xlim(ax, [-0.12 0.45]); ylim(ax, [-0.42 0.10]); axis(ax, 'off');
    title(ax, sprintf('beer = %.0f g   \\theta = %+.1f deg  (t = %.1f s)', 1000*mCup(k), th_deg(k), t(k)), 'FontWeight', 'bold');

    th = deg2rad(th_deg(k));
    sh = [0, 0]; el = [0, -Lua];
    hn = [el(1) + Lfa*cos(-th), el(2) + Lfa*sin(-th)];   % hand end, -th = sag
    % bones
    plot(ax, [sh(1) el(1)], [sh(2) el(2)], 'k-', 'LineWidth', 4);
    plot(ax, [el(1) hn(1)], [el(2) hn(2)], 'k-', 'LineWidth', 4);
    plot(ax, sh(1), sh(2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.8);
    plot(ax, el(1), el(2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9, 'LineWidth', 1.8);
    % biceps: shoulder-anchored fusiform to mid-forearm
    mid = [(el(1)+hn(1))/2, (el(2)+hn(2))/2];
    draw_muscle(ax, [0.10, -0.02], mid, 'biceps');

    % cup hangs below hand (cup axis along forearm direction)
    cupDir = [sin(-th), -cos(-th)];      % down along cup axis
    side = [cos(-th), sin(-th)];
    c1 = hn + side*(Wcup/2);
    c2 = hn - side*(Wcup/2);
    b1 = c1 + cupDir*Hcup; b2 = c2 + cupDir*Hcup;
    patch(ax, [c1(1) c2(1) b2(1) b1(1)], [c1(2) c2(2) b2(2) b1(2)], [1 1 1], 'EdgeColor', 'k', 'LineWidth', 1.8);
    % beer fill: fraction of cup height, surface stays HORIZONTAL
    fillFrac = min(mCup(k)/0.62, 1);
    if fillFrac > 0.02
        hBeer = fillFrac*Hcup;
        s1 = c1 + cupDir*(Hcup - hBeer);
        s2 = c2 + cupDir*(Hcup - hBeer);
        q1 = s1 + cupDir*hBeer; q2 = s2 + cupDir*hBeer;
        patch(ax, [s1(1) s2(1) q2(1) q1(1)], [s1(2) s2(2) q2(2) q1(2)], ...
            [0.85 0.45 0.15], 'EdgeColor', 'none');
    end
    % horizontal reference dashed line at hand height
    plot(ax, [-0.10 0.44], [hn(2) hn(2)], 'k--', 'LineWidth', 0.8);

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
function draw_muscle(ax, p1, p2, lbl)
d = p2 - p1; L = norm(d); u = d/L; nrm = [-u(2), u(1)];
w = min(0.035, L*0.35);
c = (p1 + p2)/2;
ang = atan2d(d(2), d(1));
th = linspace(0, 2*pi, 49);
ex = c(1) + (L/2)*cosd(ang)*cos(th) - 0*w*sin(th)*(-sind(ang));
% ellipse in local frame, rotated: param points
locx = (L/2)*cos(th); locy = w*sin(th);
wx = c(1) + locx*cosd(ang) - locy*sind(ang);
wy = c(2) + locx*sind(ang) + locy*cosd(ang);
patch(ax, wx, wy, [0.94 0.76 0.74], 'EdgeColor', 'k', 'LineWidth', 1.8);
text(ax, c(1), c(2), lbl, 'HorizontalAlignment', 'center', 'FontSize', 7);
end
