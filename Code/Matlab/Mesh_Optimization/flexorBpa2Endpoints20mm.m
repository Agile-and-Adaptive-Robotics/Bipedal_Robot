function [p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2)
%FLEXORBPA2ENDPOINTS20MM Second-BPA attachment points (Ben, 2026-09-21).
% The second BPA keeps the distal attachment of the old z-mirror scheme
% (pEnd{2}z = -pEnd{1}z), but its origin is no longer mirrored about the
% xy plane.  Instead the p1-to-pEnd z offset of BPA 2 equals BPA 1's:
%   p1{2}z = p1{1}z - (pEnd{1}z - pEnd{2}z)
% so both origins sit off the same side of the knee and the two routes no
% longer cross.  p1 is in the femur frame and p2 is in the t1 frame, but
% every frame in this model (femur, ICR, t1) shares the z axis: the knee
% rotates about z and all translations carry z = 0, so z coordinates are
% identical in every frame and pEnd{1}z = p2(3), pEnd{2}z = -p2(3).
% The second route itself is rebuilt by buildKneeFlexorRoute20mm so its
% wrap point and collision state are solved for this asymmetric path.

p1B = p1(:).';
p2B = p2(:).';

p1B(3) = p1(3) - 2*p2(3);   % = p1z - (pEnd1z - pEnd2z)
p2B(3) = -p2(3);            % pEnd{2} unchanged from the mirrored scheme

end
