function R = fr_rotAxis(axis, theta)
%FR_ROTAXIS Rodrigues rotation about axis through the origin.
a = axis(:) / norm(axis);
K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end
