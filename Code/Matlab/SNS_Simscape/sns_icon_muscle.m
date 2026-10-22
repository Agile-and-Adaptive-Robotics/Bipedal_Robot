function s = sns_icon_muscle(kind, nolabel)
% SNS_ICON_MUSCLE  Mask-icon Display code for the muscle-family blocks.
%
% Emits the restricted numbers-only drawing commands Simulink masks allow
% (see the gotcha note in sns_build_library.m): patch(...,[r g b]) fills,
% color('black') + plot(x, y) strokes, disp() for the centered label.
% Geometry comes from sns_icon_geom.m (side views, Ben 2026-10-21):
%   - BPA family: slightly inflated McKibben muscle, braided sleeve,
%     clamp fittings at both ends
%   - BIO: spindle belly with pennation-angle fibers into tendons
% nolabel (optional, true) omits the disp() label — used when composing
% figure crops that typeset the block name separately.

if nargin < 2, nolabel = false; end
g = sns_icon_geom(kind);
cmd = {};

% body: black backing patch (heavier outline), then the fill patch
cmd{end+1} = sprintf('patch([%s], [%s], [0 0 0]);', ...
    vecstr(g.bodyX), vecstr(g.bodyY + edgeSign(g.bodyY, g.edge)));
cmd{end+1} = sprintf('patch([%s], [%s], [%.3g %.3g %.3g]);', ...
    vecstr(g.bodyX), vecstr(g.bodyY), g.fill);

if isfield(g, 'braid')
    cmd{end+1} = 'color(''black'');';
    for i = 1:size(g.braid, 1)
        cmd{end+1} = sprintf('plot([%.4g %.4g], [%.4g %.4g]);', g.braid(i, :));
    end
    caps = {g.capR, g.capL};
    for c = 1:2
        cp = caps{c};
        cmd{end+1} = sprintf('patch([%s], [%s], [0 0 0]);', ...
            vecstr(cp(:,1)), vecstr(cp(:,2) + 0.018 * sign(cp(:,2))));
        cmd{end+1} = sprintf('patch([%s], [%s], [%.3g %.3g %.3g]);', ...
            vecstr(cp(:,1)), vecstr(cp(:,2)), g.cap);
    end
elseif isfield(g, 'fibers')
    cmd{end+1} = 'color(''black'');';
    for i = 1:size(g.fibers, 1)
        cmd{end+1} = sprintf('plot([%.4g %.4g], [%.4g %.4g]);', g.fibers(i, :));
    end
    tns = {g.tendonR, g.tendonL};
    for c = 1:2
        tp = tns{c};
        cmd{end+1} = sprintf('patch([%s], [%s], [0 0 0]);', ...
            vecstr(tp(:,1)), vecstr(tp(:,2) + 0.018 * sign(tp(:,2))));
        cmd{end+1} = sprintf('patch([%s], [%s], [%.3g %.3g %.3g]);', ...
            vecstr(tp(:,1)), vecstr(tp(:,2)), g.tendon);
    end
end
if ~nolabel
    cmd{end+1} = sprintf('disp(''%s'');', g.label);
end
s = strjoin(cmd, newline);
end

function t = vecstr(v)
t = strjoin(arrayfun(@(x) sprintf('%.4g', x), v, ...
    'UniformOutput', false), ' ');
end

function y = edgeSign(y, e)
% expand the closed body polygon outward by e (sign-aware)
y = y + e * sign(y);
y(y == 0) = e;   % degenerate flat points still get the outline
end
