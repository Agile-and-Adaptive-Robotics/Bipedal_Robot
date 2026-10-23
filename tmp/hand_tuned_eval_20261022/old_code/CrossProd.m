function ma = CrossProd(A,B,n,Rot)
% Takes cross product of A and B and returns scalar value for dimension n
% Rot is a rotation matrix for joints like Ankle, Subtalar, and MTP

v = cross(A,B);
switch nargin
 case 3
    ma = v(n);
 case 4
    u = Rot*v;
    ma = u(n);
 otherwise
    ma = 'Error calculating moment arm';
    error(ma)
 end
end