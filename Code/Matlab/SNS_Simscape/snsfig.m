function h = snsfig(cmd, varargin)
% SNSFIG  Drawing primitives for journal-style SNS circuit diagrams.
%
% Diagram language (Szczecinski et al. 2017 Fig. 2 / Rybak & Shevtsova /
% Animatlab conventions, per Ben 2026-09-09):
%   neuron      = open circle, black edge, white fill
%   afferent    = open circle with Ia/Ib label
%   muscle      = fusiform ellipse (light green, Okabe-Ito tint)
%   EXCITATORY  = white triangle with black edges, tip pointing into target
%   INHIBITORY  = solid black circle
% Shape is the primary code (colorblind-safe by construction); accents use
% the Okabe-Ito CVD-safe palette. Draw in data units with axis equal.
%
% Commands:
%   snsfig('neuron', x, y, r, label)          open circle + centered label
%   snsfig('neuron', x, y, r, label, sub)     + sublabel below circle
%   snsfig('afferent', x, y, r, 'Ia'|'Ib', sub)
%   snsfig('muscle', x, y, w, ht, label)      fusiform ellipse
%   snsfig('box', x, y, w, ht, label)         rounded rectangle (plant etc.)
%   snsfig('edge', x1,y1, x2,y2, style)       plain segment, style '-','--'
%   snsfig('edgeR', xs,ys, xd,yd, r)          segment trimmed at circle radius r
%   snsfig('exc', x, y, angDeg, s)            open triangle pointing angDeg
%   snsfig('inh', x, y, s)                    filled black dot
%   snsfig('arrow', x1,y1, x2,y2, style)      segment with arrowhead at end
%   snsfig('label', x, y, txt, size, style)   text, style 'n','b','i','g'
% All return the graphics handle of the main object.

switch lower(cmd)
    case 'neuron'
        x = varargin{1}; y = varargin{2}; r = varargin{3};
        lbl = varargin{4}; sub = '';
        if numel(varargin) >= 5, sub = varargin{5}; end
        th = linspace(0, 2*pi, 73);
        h = patch(x + r*cos(th), y + r*sin(th), [1 1 1], 'EdgeColor', 'k', 'LineWidth', 1.1);
        if ~isempty(lbl)
            text(x, y, lbl, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 8);
        end
        if ~isempty(sub)
            text(x, y - r - 1.6, sub, 'HorizontalAlignment', 'center', 'FontSize', 7.5, 'FontAngle', 'italic');
        end
    case 'afferent'
        x = varargin{1}; y = varargin{2}; r = varargin{3};
        lbl = varargin{4}; sub = '';
        if numel(varargin) >= 5, sub = varargin{5}; end
        th = linspace(0, 2*pi, 73);
        h = patch(x + r*cos(th), y + r*sin(th), [0.95 0.95 0.95], 'EdgeColor', 'k', 'LineWidth', 1.1);
        text(x, y, lbl, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 8);
        if ~isempty(sub)
            text(x, y - r - 1.6, sub, 'HorizontalAlignment', 'center', 'FontSize', 7.5, 'FontAngle', 'italic');
        end
    case 'muscle'
        x = varargin{1}; y = varargin{2}; w = varargin{3}; ht = varargin{4};
        lbl = varargin{5};
        th = linspace(0, 2*pi, 73);
        h = patch(x + (w/2)*cos(th), y + (ht/2)*sin(th), [0.82 0.92 0.87], ...
            'EdgeColor', 'k', 'LineWidth', 1.1);
        text(x, y, lbl, 'HorizontalAlignment', 'center', 'FontSize', 7.5);
    case 'box'
        x = varargin{1}; y = varargin{2}; w = varargin{3}; ht = varargin{4};
        lbl = varargin{5};
        h = rectangle('Position', [x - w/2, y - ht/2, w, ht], 'Curvature', 0.18, ...
            'FaceColor', [0.93 0.93 0.93], 'EdgeColor', 'k', 'LineWidth', 1.1);
        if contains(lbl, newline)
            text(x, y, lbl, 'HorizontalAlignment', 'center', 'FontSize', 8);
        else
            text(x, y, lbl, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 8);
        end
    case 'edge'
        x1 = varargin{1}; y1 = varargin{2}; x2 = varargin{3}; y2 = varargin{4};
        st = '-';
        if numel(varargin) >= 5, st = varargin{5}; end
        h = line([x1 x2], [y1 y2], 'Color', 'k', 'LineStyle', st, 'LineWidth', 1.1);
    case 'edgeR'
        % segment from source point, stopping at radius r around destination
        xs = varargin{1}; ys = varargin{2}; xd = varargin{3}; yd = varargin{4};
        r = varargin{5}; st = '-'; lw = 1.1;
        if numel(varargin) >= 6, st = varargin{6}; end
        if numel(varargin) >= 7, lw = varargin{7}; end
        d = [xd - xs, yd - ys]; L = norm(d); u = d / L;
        xe = xd - u(1)*r; ye = yd - u(2)*r;
        h = line([xs xe], [ys ye], 'Color', 'k', 'LineStyle', st, 'LineWidth', lw);
    case 'exc'
        % white triangle with black edges. INVERTED per Ben 2026-09-09:
        % the tip points OPPOSITE to angDeg — i.e. back along the connection,
        % with the flat base toward the target (postsynaptic) side.
        x = varargin{1}; y = varargin{2}; ang = varargin{3}; s = varargin{4};
        u = [cosd(ang), sind(ang)]; n = [-sind(ang), cosd(ang)];
        tip = [x, y] - u*(s/2);
        b1 = [x, y] + u*(s/2) + n*(s*0.55);
        b2 = [x, y] + u*(s/2) - n*(s*0.55);
        h = patch([tip(1) b1(1) b2(1)], [tip(2) b1(2) b2(2)], [1 1 1], ...
            'EdgeColor', 'k', 'LineWidth', 1.1);
    case 'inh'
        x = varargin{1}; y = varargin{2}; s = varargin{3};
        th = linspace(0, 2*pi, 49);
        h = patch(x + (s/2.2)*cos(th), y + (s/2.2)*sin(th), [0 0 0], 'EdgeColor', 'k');
    case 'arrow'
        x1 = varargin{1}; y1 = varargin{2}; x2 = varargin{3}; y2 = varargin{4};
        st = '-'; col = 'k';
        if numel(varargin) >= 5, st = varargin{5}; end
        if numel(varargin) >= 6, col = varargin{6}; end
        h = annotation_ne_arrow(x1, y1, x2, y2, st, col);
    case 'label'
        x = varargin{1}; y = varargin{2}; txt = varargin{3};
        fs = 8; st = 'n';
        if numel(varargin) >= 4, fs = varargin{4}; end
        if numel(varargin) >= 5, st = varargin{5}; end
        opts = {'HorizontalAlignment', 'center', 'FontSize', fs, 'Interpreter', 'tex'};
        switch st
            case 'b', opts = [opts, {'FontWeight', 'bold'}];
            case 'i', opts = [opts, {'FontAngle', 'italic'}];
            case 'g', opts = [opts, {'Color', [0.35 0.35 0.35]}];
            case 'gi', opts = [opts, {'Color', [0.35 0.35 0.35], 'FontAngle', 'italic'}];
        end
        h = text(x, y, txt, opts{:});
    otherwise
        error('snsfig: unknown command %s', cmd);
end
end

function h = annotation_ne_arrow(x1, y1, x2, y2, st, col)
% Small filled arrowhead in data units (annotation objects avoided on purpose).
d = [x2 - x1, y2 - y1]; L = norm(d); u = d / L; n = [-u(2), u(1)];
ah = min(1.8, L/4); aw = ah * 0.55;
bx = x2 - u(1)*ah; by = y2 - u(2)*ah;
h = line([x1 bx], [y1 by], 'Color', col, 'LineStyle', st, 'LineWidth', 1.1);
patch([x2 bx + n(1)*aw bx - n(1)*aw], [y2 by + n(2)*aw by - n(2)*aw], col, 'EdgeColor', col);
end
