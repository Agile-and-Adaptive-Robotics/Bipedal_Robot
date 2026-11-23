function [force, fiberRate] = thelen_reflex_fiber(activation, pathLength, fiberLength, p)
% Elastic tendon and fixed-height pennation. Units are N, m, and s.
% OpenSim 4.4 Thelen2003Muscle::calcfse, calcfpe, and calcdlceN.
% p = F0, lopt, lslack, pennation angle (rad), vmax, tendon strain,
% passive strain, passive and active curve shape, Af, Flen, FV extrapolation, amin.
%#codegen
height = p(2)*sin(p(4));
minFiber = max(0.01*p(2),height/sin(acos(0.1)));
fiber = max(minFiber,fiberLength);
projected = sqrt(max(fiber^2-height^2,0));
cosPennation = projected/fiber;
strain = (pathLength-projected)/p(3)-1;
toe = 99*exp(3)/(166*exp(3)-67)*p(6);
if strain <= 0
    tendon = 0;
elseif strain < toe
    tendon = 0.33*expm1(3*strain/toe)/expm1(3);
else
    tendon = 0.33+0.67*(strain-toe)/(p(6)-toe);
end
force = p(1)*tendon;
normalizedLength = fiber/p(2);
activeLength = exp(-(normalizedLength-1)^2/p(9));
passive = expm1(p(8)*max(normalizedLength-1,0)/p(7))/expm1(p(8));
a = min(1,max(p(13),activation));
activeIsometric = a*activeLength;
activeForce = tendon/cosPennation-passive;
% Outside the original concentric/eccentric branches, OpenSim extends
% the tangent from its boundary to keep velocity finite.
threshold = p(12)*p(11)*activeIsometric;
edgeForce = min(threshold,max(0,activeForce));
if edgeForce <= activeIsometric
    denom = activeIsometric+edgeForce/p(10);
    deriv = 1/p(10);
else
    scale = (2+2/p(10))/(p(11)-1);
    denom = scale*(activeIsometric*p(11)-edgeForce);
    deriv = -scale;
end
speedScale = 0.25+0.75*a;
edgeSpeed = speedScale*(edgeForce-activeIsometric)/denom;
slope = (speedScale-edgeSpeed*deriv)/denom;
fiberRate = (edgeSpeed+slope*(activeForce-edgeForce))*p(5)*p(2);
if fiberLength <= minFiber && fiberRate < 0
    fiberRate = 0;
end
end
