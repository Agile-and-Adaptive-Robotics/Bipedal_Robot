All of the functions herein are from the Modern Robotics Toolbox. https://github.com/NxRLab/ModernRobotics/tree/master/packages/MATLAB/mr
Adjoint: is used to change the reference frame of a Wrench.

RowVecTrans:
% Takes 4x4 transformation matrix T and 3x1 row vector p.
% Adds row to the column matrix p to get r
% Multiplies transform T and vector r together to get vector s
% Returns 3x1 column taking first three rows of vector s

RpToTrans takes rotation matrix R and position p. Returns the corresponding homogeneous transformation matrix T.

TransInv gives the inverse transformation matrix.

TransToRp takes transformation matrix and returns two outputs: the rotation matrix and the corresponding position vector.

VeCToso3 is used by the Adjoint function.:
% Takes a 3-vector (angular velocity).
% Returns the skew symmetric matrix in so(3).
