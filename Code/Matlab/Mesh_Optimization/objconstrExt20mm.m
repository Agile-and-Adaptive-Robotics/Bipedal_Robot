function out = objconstrExt20mm(x, obj, nonlcon)
% objconstrExt20mm
%
% surrogateopt objective/constraint wrapper.

out.Fval = obj(x);

[c, ~] = nonlcon(x);

out.Ineq = c(:);

end