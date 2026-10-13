function [xFull, ctxCam] = camToExtDesign20mm(xCam, ctx, xFixed)
% camToExtDesign20mm
%
% Map the 8-variable clip-on-cam design onto the 8-variable extensor
% family design + context.
%
% xCam   = [p6x p6y p7x p7y p8x p8y tendon p3dx]  (p6:p8 t1 frame, m;
%          p3dx = femur-frame +X offset of the p3 clip, m, Ben
%          2026-10-13: 25-50 mm)
% xFixed = [p1(1:3), pEnd(1:3), rest]             (loaded Opt_run_Ext winner)
%
% The cam points enter as explicit seed overrides (geo.seedP6XY /
% seedP7XY / seedP8XY / seedP3XY) consumed by adjustedSeed20mm inside
% buildDistalRingLocation20mm; the stock circle/wall construction and
% the fcec-rule p3 are bypassed for those rows. p4:p5 still follow the
% clearance-driven femur-chain repair. Cam designs also set
% geo.tibiaWrapCircumcircle so the tibia wrap-loss radius is the
% circumcircle of the actual p6:p8 contacts (Ben, 2026-10-13:
% recalculated bend radius), not the legacy fixed CAD center.

ctxCam = ctx;
ctxCam.geo.seedP6XY = reshape(xCam(1:2),1,2);
ctxCam.geo.seedP7XY = reshape(xCam(3:4),1,2);
ctxCam.geo.seedP8XY = reshape(xCam(5:6),1,2);

p3Base = ctx.routeSeed(3,1:2);
ctxCam.geo.seedP3XY = [p3Base(1) + xCam(8), p3Base(2)];

ctxCam.geo.tibiaWrapCircumcircle = true;

xFull = [reshape(xFixed,1,[]), xCam(7)];

end
