function [c, ceq] = nonlconExt20mm(x, ctx)
% nonlconExt20mm
%
% Nonlinear geometry constraints for the 20 mm extensor route.
%
% c <= 0 is feasible.
%
% x = [p1(1:3), p8(1:3), rest, tendon]

p8 = x(4:6);
tendon = x(8);

geo = ctx.geo;

ceq = [];


%% Geometry-dependent tendon maximum

[~, tendonMaxGeom, ~, ~, cInside] = ...
    tendonLimit20mm(p8, geo);


%% Constraints
%
% 1. Tendon cannot extend farther than the point where its straight
%    line reaches the BPA clearance envelope.
%
% 2. There must be enough geometrically available length for the
%    minimum 25 mm tendon.
%
% 3. If a circular tangent is being used, p8 must permit a real
%    external tangent.

cTendonMax = tendon - tendonMaxGeom;

cMinimumPossible = geo.tendonMin - tendonMaxGeom;

c = [ ...
    cTendonMax;
    cMinimumPossible;
    cInside];

end