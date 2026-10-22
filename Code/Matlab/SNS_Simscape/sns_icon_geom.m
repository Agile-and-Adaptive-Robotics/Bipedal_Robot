function g = sns_icon_geom(kind)
% SNS_ICON_GEOM  Side-view geometry for the muscle-family library icons.
%
% Ben, 2026-10-21: the muscle symbols must read as SIDE VIEWS, not
% cross-sections. The BPA is a slightly inflated McKibben muscle: a
% barrel-shaped braided sleeve with clamp fittings at both ends. The
% biological muscle is a spindle belly with a pennation angle: oblique
% fibers running between aponeuroses into tapering tendons.
%
% Output is plain vertex data in mask-autoscale coordinates; both the
% Simulink mask-icon emitter (sns_icon_muscle.m) and the GUI canvas draw
% from it so the symbol stays identical everywhere.
%
% ASPECT MATTERS: mask icons stretch to the block rectangle, so the
% muscle blocks must be WIDE (aspect ~2.7-2.8, the drawing's natural
% extent) or the side view squashes back into a blob.

switch upper(kind)
    case 'BIO'
        g.fill = [0.87 0.62 0.60];              % darker red, existing BIO
        g.tendon = [0.88 0.88 0.88];
        g.blockAspect = 2.6;                    % w/h for unstretched look
        % spindle belly, rounded taper
        xb = linspace(-0.58, 0.58, 49);
        rb = 0.335 * (1 - (xb / 0.58).^2).^0.72;
        % pennation fibers: 28 deg from vertical, inset to stay inside
        f = [];
        for i = 1:7
            x1 = -0.46 + 0.155 * (i - 1);
            dx = 2 * rOf(x1) * tan(28 * pi / 180);
            x2 = min(x1 + dx, 0.52);
            f = [f; x1 -0.90 * rOf(x1) x2 0.90 * rOf(x2)]; %#ok<AGROW>
        end
        g.fibers = f;
        % tapering tendons (half height 0.13 -> 0.05)
        g.tendonR = [0.58 0.95 0.95 0.58; -0.13 -0.05 0.05 0.13]';
        g.tendonL = [-0.58 -0.95 -0.95 -0.58; -0.13 -0.05 0.05 0.13]';
        g.bodyX = [xb, fliplr(xb)];
        g.bodyY = [rb, -fliplr(rb)];
        g.edge = 0.022;                          % outline thickness
        g.label = 'BIO';
    otherwise                                    % McKibben BPA family
        g.fill = [0.94 0.76 0.74];              % light salmon, existing BPA
        g.cap = [0.55 0.55 0.55];
        g.blockAspect = 2.75;
        % slightly inflated barrel: radius swells toward mid-length
        xb = linspace(-0.70, 0.70, 49);
        rb = 0.185 + 0.115 * (1 - (xb / 0.70).^2);
        % braided sleeve: two crossing helix families; endpoints inset
        % 10% so the strokes stay clear of the outline ring
        br = [];
        for i = 1:7
            x1 = -0.64 + 0.175 * (i - 1);
            x2 = min(x1 + 0.30, 0.66);
            br = [br; x1 -0.90 * rOf(x1) x2 0.90 * rOf(x2)];  %#ok<AGROW>
            br = [br; x1  0.90 * rOf(x1) x2 -0.90 * rOf(x2)]; %#ok<AGROW>
        end
        g.braid = br;
        % clamp fittings at both ends
        g.capR = [0.70 0.90 0.90 0.70; -0.20 -0.20 0.20 0.20]';
        g.capL = [-0.70 -0.90 -0.90 -0.70; -0.20 -0.20 0.20 0.20]';
        g.bodyX = [xb, fliplr(xb)];
        g.bodyY = [rb, -fliplr(rb)];
        g.edge = 0.025;
        g.label = kind;                           % 'BPA', 'BPA 10', ...
end

    function r = rOf(xq)                          % radius at arbitrary x
        switch upper(kind)
            case 'BIO'
                a = 0.58; r0 = 0.335; p = 0.72;
                r = r0 * max(1 - (xq / a).^2, 0).^p;
            otherwise
                a = 0.70; r0 = 0.185; b = 0.115;
                r = r0 + b * max(1 - (xq / a).^2, 0);
        end
    end
end
