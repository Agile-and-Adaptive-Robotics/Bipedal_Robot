function [c, ceq] = nonlconExtCam20mm(xCam, ctx, xFixed, vtugFloorKPa)
% nonlconExtCam20mm
%
% Hard constraints for the clip-on tibia cam design. c <= 0 feasible.
%
% xCam   = [p6x p6y p7x p7y p8x p8y tendon]
% xFixed = [p1(1:3), pEnd(1:3), rest]
% vtugFloorKPa = Festo VTUG valve-manifold floor (default 150 kPa per
%                the datasheet, internal pilot air supply).
%
% Family set on the mapped design (nonlconExt20mm): tendon geometry
% limits, p8 -x feasibility, cRestLength (the NO-STRETCH rule: at 0
% pressure the BPA sits at or below rest length everywhere in the RoM),
% and contraction <= KMAX.
%
% Cam additions (Ben, 2026-10-12 / 2026-10-13):
%   - p6/p7/p8 stay ordered proximal -> distal in y;
%   - no cam point may sit inside the stock tibia clearance circles;
%   - p6 and p7 sit within 15 mm of each other in X, and p6 sits 20-45
%     mm above p7 (rounded-corner clip manufacturability);
%   - p6 must not interfere with the femoral condyles at 10 deg
%     extension (the t1-side point in the femur frame at the
%     most-extended solver pose, against the inflated condyle
%     clearance envelope);
%   - VTUG FLOOR: the design must not exceed the human target AT the
%     floor pressure (below it the manifold cannot throttle).
%
% RETIRED (Ben, 2026-10-13): the theta(p6) angular restriction, the
% p3 clip (p3:p5 are back on their fcec-rule seeds), the p8
% contact-forcing binary and chord-interference window (p8 is a CAD
% design point that carries wrap loss when contacted, nothing more),
% and the hard arm-step cap (a cap cannot fix release math -- release
% continuity is verified in the winner's arm curve instead; the soft
% arm-step penalty stays in the objective as seed doctrine).
%
% p8's Y lives in the driver's +/-5 mm box around its constructed seed.

if nargin < 4 || isempty(vtugFloorKPa)
    vtugFloorKPa = 150;
end

[xFull, ctxCam] = camToExtDesign20mm(xCam, ctx, xFixed);

[cBase, ~] = nonlconExt20mm(xFull, ctxCam);

p6 = xCam(1:2);
p7 = xCam(3:4);
p8 = xCam(5:6);

cOrder = [p7(2) - p6(2);
          p8(2) - p7(2)];

geo = ctx.geo;

cInsideUpper = [ ...
    geo.tibiaUpperClearRadius - norm(p6 - geo.tibiaUpperCenter);
    geo.tibiaUpperClearRadius - norm(p7 - geo.tibiaUpperCenter)];

cInsideLower = [ ...
    geo.tibiaLowerClearRadius - norm(p8 - geo.tibiaLowerCenter);
    geo.tibiaLowerClearRadius - norm(p7 - geo.tibiaLowerCenter)];

% p6/p7 X window: |p7x - p6x| <= 15 mm, two smooth rows.
xWindow = 0.015;
cP6P7X = [p7(1) - p6(1) - xWindow;
          p6(1) - p7(1) - xWindow];

% p6/p7 height band: p6 sits 20-45 mm above p7.
yBand = [0.020, 0.045];
cP6P7Y = [yBand(1) - (p6(2) - p7(2));
          (p6(2) - p7(2)) - yBand(2)];

% p8 (Ben, 2026-10-14): NOT a constraint. p8 is a CAD design point
% whose slight protrusion past the pEnd -> p7 chord (he sketches
% roughly 1-2 mm of added path length) makes the BPA bend at it and
% contributes wrap/arc-length loss. Whether the route model contacts
% it follows from the geometry and the release logic -- nothing here
% forces it.

%% p6 vs femoral condyles at 10 deg extension
[~, idxExtPose] = max(ctx.phiD);
T6ext = ctx.T_Pam(:,:,idxExtPose)*ctx.T_ICR_t1(:,:,idxExtPose);
p6FemurExt = RowVecTrans(T6ext, [p6, 0]);
cP6Condyle = pointClearOfPolygon(p6FemurExt, geo.femurOffsetBoundary);

%% VTUG floor overshoot (scaled by ctx.torqueScale; <= 0 feasible)
cVtug = 1;
if vtugFloorKPa > 0
    try
        restV = xFull(7);
        kmaxV = restV*(1 - ctxCam.KMAX);
        [LocationV, bendV] = buildDistalRingLocation20mm( ...
            xFull(1:3), xFull(4:6), xFull(8), ctxCam);
        bpaV = MonoPamDataExplicit_balanceX3(ctxCam.Name, LocationV, ...
            ctxCam.CrossPoint, ctxCam.Dia, ctxCam.T_Pam, restV, kmaxV, ...
            xFull(8), ctxCam.fitting, vtugFloorKPa, ctxCam.Xi0, ...
            ctxCam.Xi1, ctxCam.Xi2, ctxCam.Xi3, ctxCam.wraps, ...
            ctxCam.phiD, ctxCam.BPAcount, bendV);
        humanAbsV = interp1(ctx.humanAngleD, ctx.humanTorqueAbs, ...
            ctx.phiD, 'pchip', 'extrap');
        hasTgtV = isfinite(humanAbsV) & humanAbsV > 0;
        tqV = abs(bpaV.Torque_p(:,3));
        if any(hasTgtV) && all(isfinite(tqV(hasTgtV)))
            cVtug = max((tqV(hasTgtV) - humanAbsV(hasTgtV))./ctx.torqueScale);
        end
    catch
        cVtug = 1;   % route/class failure = infeasible, family style
    end
end

ceq = [];

% Tail order is load-bearing for Opt_run_ExtCam's report print:
% cP6Condyle = cBest(end-1), cVtug = cBest(end).
c = [cBase(:); cOrder; cInsideUpper; cInsideLower; ...
    cP6P7X; cP6P7Y; cP6Condyle; cVtug];

end


function c = pointClearOfPolygon(pt, poly)
% <= 0 when the point lies outside the closed polygon (binary, like the
% family's cInside constraints; on-boundary tangency counts as clear).

[in, on] = inpolygon(pt(1), pt(2), poly(:,1), poly(:,2));

if in && ~on
    c = 1;
else
    c = -1e-6;
end

end
