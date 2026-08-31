function [c, ceq] = nonlconExt20mm(x, ctx)
% nonlconExt20mm
%
% Nonlinear geometry constraints for the 20 mm extensor route.
%
% c <= 0 is feasible.
%
% x = [p1(1:3), pEnd(1:3), rest, tendon]

pEnd = x(4:6);
tendon = x(8);

geo = ctx.geo;

ceq = [];


%% Geometry-dependent tendon maximum

[~, tendonMaxGeom, ~, ~, cInside] = ...
    tendonLimit20mm(pEnd, geo);

pred = predictKneeExt20mm(x, ctx);

if pred.ok
    p8T1 = pred.routeInfo.seed20(8,1:2);
    cP8MinX = p8MinusXLimit20mm(p8T1, geo);

    % Hard length feasibility: the longest undeformed route must not exceed
    % the zero-strain musculotendon length.
    cRestLength = pred.cRestLength;

    % Hard maximum-contraction feasibility using the class's original
    % geometric Contraction definition, exactly as requested.
    relativeContraction = pred.bpa.Contraction(:)/pred.KMAX;
    cContractionMax = max(relativeContraction) - 1;
else
    cP8MinX = 1;
    cRestLength = 1;
    cContractionMax = 1;
end


%% Constraints
%
% 1. Tendon cannot extend farther than the point where its straight
%    line reaches the BPA clearance envelope.
%
% 2. There must be enough geometrically available length for the
%    minimum 25 mm tendon.
%
% 3. If a circular tangent is being used, pEnd must permit a real
%    external tangent.
%
% 4. The longest undeformed route must fit within the zero-strain series
%    length.
%
% 5. bpa.Contraction/KMAX must not exceed one at any knee position.

cTendonMax = tendon - tendonMaxGeom;

cMinimumPossible = geo.tendonMin - tendonMaxGeom;

c = [ ...
    cTendonMax;
    cMinimumPossible;
    cInside;
    cP8MinX;
    cRestLength;
    cContractionMax];

end


function c = p8MinusXLimit20mm(p8T1, geo)
% Positive c means p8 is too far in the -X direction.  The limit is a
% horizontal offset from the sloped CAD reference line in the t1 frame.

P1 = geo.p8MinXLineT1(1,:);
P2 = geo.p8MinXLineT1(2,:);

if abs(P2(2) - P1(2)) < 1e-12
    xLine = max(P1(1), P2(1));
else
    t = (p8T1(2) - P1(2)) / (P2(2) - P1(2));
    xLine = P1(1) + t*(P2(1) - P1(1));
end

xMin = xLine - geo.p8MinXMinusOffset;
c = xMin - p8T1(1);

end
