function R = fr_xyzExtrinsic(a)
%FR_XYZEXTRINSIC frame-fixed x-y-z rotation sequence (OpenSim
% PhysicalOffsetFrame orientation convention): fixed axes x, then y, then z.
if numel(a) < 3; a = [0 0 0]; end
R = fr_rotAxis([0;0;1], a(3)) * fr_rotAxis([0;1;0], a(2)) * ...
    fr_rotAxis([1;0;0], a(1));
end
