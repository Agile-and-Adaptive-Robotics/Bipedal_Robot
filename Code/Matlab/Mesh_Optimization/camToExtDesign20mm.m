function [xFull, ctxCam] = camToExtDesign20mm(xCam, ctx, xFixed)
% camToExtDesign20mm
%
% Map the 7-variable clip-on-cam design onto the 8-variable extensor
% family design + context.
%
% xCam   = [p6x p6y p7x p7y p8x p8y tendon]  (p6:p8 t1 frame, meters)
% xFixed = [p1(1:3), pEnd(1:3), rest]         (loaded Opt_run_Ext winner)
%
% The cam points enter as explicit t1 seed overrides (geo.seedP6XY /
% seedP7XY / seedP8XY) consumed by adjustedSeed20mm inside
% buildDistalRingLocation20mm; the stock circle/wall construction is
% bypassed for those rows. p3:p5 still follow p6 through the
% clearance-driven femur-chain repair.

ctxCam = ctx;
ctxCam.geo.seedP6XY = reshape(xCam(1:2),1,2);
ctxCam.geo.seedP7XY = reshape(xCam(3:4),1,2);
ctxCam.geo.seedP8XY = reshape(xCam(5:6),1,2);

xFull = [reshape(xFixed,1,[]), xCam(7)];

end
