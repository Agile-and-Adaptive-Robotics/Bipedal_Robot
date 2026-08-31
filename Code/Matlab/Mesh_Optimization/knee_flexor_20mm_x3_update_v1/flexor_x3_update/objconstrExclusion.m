function out = objconstrExclusion(x, obj, geo, ctx, idxP2)
%OBJCONSTREXCLUSION Objective/constraint wrapper for surrogateopt.
%
% surrogateopt wants:
%   out.Fval
%   out.Ineq <= 0

out.Fval = obj(x);

[c, ~] = nonlconExclusion(x, geo, ctx, idxP2);

out.Ineq = c;

end
