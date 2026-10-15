function J = objective_KneeExtCam20mm(xCam, ctx, xFixed, refPts)
% objective_KneeExtCam20mm
%
% Clip-on tibia cam objective (Ben, 2026-10-12: "find a cam shape that
% maximizes moment arm and torque when I need it").
%
% xCam   = [p6x p6y p7x p7y p8x p8y tendon]
% xFixed = [p1(1:3), pEnd(1:3), rest]
% refPts = 3x2 first-cut reference contacts (the 25.4 mm radial push of
%          p6/p7; p8 at its constructed position).
%
% Composition: the full family objective_KneeExt20mm on the mapped
% design (torque requirement vs the human target, strain/contraction
% penalties, hard-term shadows), plus
%   (1) a capped reward on the worst-angle torque margin vs the
%       REQUIRED target so the solver keeps buying moment arm instead
%       of stopping at the +5% requirement;
%   (2) a moment-arm step penalty across adjacent solver poses (Ben's
%       release-discontinuity doctrine; 2.2 mm accepted reference);
%   (3) a mild regularizer toward the first-cut push neighborhood.
% The family objective's small anti-overshoot term (weight 1e-3) stays
% inside the base value; it is six orders below the arm reward and
% cannot oppose the cam study's purpose.
%
% Cost note: this calls objective_KneeExt20mm AND predictKneeExt20mm,
% i.e. two full predictions per evaluation (the family driver already
% spends one more in nonlconExt20mm).

[xFull, ctxCam] = camToExtDesign20mm(xCam, ctx, xFixed);

Jbase = objective_KneeExt20mm(xFull, ctxCam);

pred = predictKneeExt20mm(xFull, ctxCam);

if ~pred.ok || any(~isfinite([pred.TorqueZ(:); pred.momentArm(:)]))
    J = 1e12;
    return
end

humanAbs = interp1( ...
    ctx.humanAngleD, ...
    ctx.humanTorqueAbs, ...
    ctx.phiD, ...
    'pchip', ...
    'extrap');

robotAbs = abs(pred.TorqueZ);

requiredTorque = (1 + ctx.requiredTorqueMargin)*humanAbs(:);

marginVsRequired = (robotAbs(:) - requiredTorque(:))./ctx.torqueScale;

% Reward the worst-angle margin, capped at +100% of the torque scale.
Jarm = -1e2*min([1; marginVsRequired]);

% Moment-arm discontinuity penalty (pose spacing 1.31 deg; the accepted
% release step from the 9/22 sweep is ~2.2 mm).
armStepRef = 0.0022;
maxArmStep = max(abs(diff(pred.momentArm)));
JarmStep = 1e3*max(0, maxArmStep - armStepRef)^2/(armStepRef^2);

% Mild tie-break toward the first-cut push neighborhood.
camPts = reshape(xCam(1:6),3,2);
Jreg = 1e-2*mean(sum(((camPts - refPts)/0.0254).^2, 2));

J = Jbase + Jarm + JarmStep + Jreg;

if ~isfinite(J)
    J = 1e12;
end

end
