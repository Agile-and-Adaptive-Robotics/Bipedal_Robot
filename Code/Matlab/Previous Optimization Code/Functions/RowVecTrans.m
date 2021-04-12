function v = RowVecTrans(T,p)
% Takes 4x4 transformation matrix R and 3x1 column vector p.
% Adds row to the column matrix p to get r
% Multiplies transform T and vector r together to get vector s
% Returns 3x1 column taking first three rows of vector s

r = [p, 1];
s = T*r';
v = [s(1), s(2), s(3)];
end