function [c, ceq] = nonlconExtCam20mm(xCam, ctx, xFixed, thetaP6BaseD, vtugFloorKPa)
% nonlconExtCam20mm
%
% Hard constraints for the clip-on tibia cam design. c <= 0 feasible.
%
% xCam   = [p6x p6y p7x p7y p8x p8y tendon]
% xFixed = [p1(1:3), pEnd(1:3), rest]
% thetaP6BaseD = t1-frame polar angle of the CONSTRUCTED baseline p6.
% vtugFloorKPa = Festo VTUG valve-manifold floor (default 200 kPa).
%
% Family set on the mapped design (nonlconExt20mm): tendon geometry
% limits, p8 -x feasibility, cRestLength (the NO-STRETCH rule: at 0
% pressure the BPA sits at or below rest length everywhere in the RoM,
% so it cannot produce torque unpressurized), and contraction <= KMAX.
%
% Cam additions (Ben, 2026-10-12):
%   - the t1-frame polar angle of the line p6-t1 must not increase;
%   - p6/p7/p8 stay ordered proximal -> distal in y;
%   - no cam point may sit inside the stock tibia clearance circles;
%   - VTUG FLOOR: the manifold does not like to open/close below
%     ~200 kPa, so the design must not exceed the human target AT the
%     floor pressure -- if it does, the required pressure at that angle
%     falls below the floor and the valve cannot throttle down to the
%     target (this is what limits the cam size, not the search boxes).

if nargin < 5 || isempty(vtugFloorKPa)
    vtugFloorKPa = 200;
end

[xFull, ctxCam] = camToExtDesign20mm(xCam, ctx, xFixed);

[cBase, ~] = nonlconExt20mm(xFull, ctxCam);

p6 = xCam(1:2);
p7 = xCam(3:4);
p8 = xCam(5:6);

cTheta = atan2d(p6(2),p6(1)) - thetaP6BaseD;

cOrder = [p7(2) - p6(2);
          p8(2) - p7(2)];

geo = ctx.geo;

cInsideUpper = [ ...
    geo.tibiaUpperClearRadius - norm(p6 - geo.tibiaUpperCenter);
    geo.tibiaUpperClearRadius - norm(p7 - geo.tibiaUpperCenter)];

cInsideLower = [ ...
    geo.tibiaLowerClearRadius - norm(p8 - geo.tibiaLowerCenter);
    geo.tibiaLowerClearRadius - norm(p7 - geo.tibiaLowerCenter)];

%% VTUG floor overshoot (scaled by ctx.torqueScale; <= 0 feasible)
cVtug = 1;   % infeasible unless the floor check runs clean
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

c = [cBase(:); cTheta; cOrder; cInsideUpper; cInsideLower; cVtug];

end
