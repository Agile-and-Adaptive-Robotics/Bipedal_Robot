function ok = fr_wrappingCheck(C, kin, cases)
%FR_WRAPPINGCHECK Gate check (c): wrapping transitions.
%
% gait2392 contains NO wrap objects (geometry_20261103 README: 13 empty
% WrapObjectSets per model); the functional equivalent is the
% ConditionalPathPoint census with declared knee/hip ranges. This check
% verifies the conditional activation BEHAVIOR implemented in
% fr_evalRoute: the flags flip exactly at the declared ranges and the
% evaluated geometry responds (length jump quantified). The grid
% mechanics rows are fixed per the documented stage-1 simplification, so
% no torque discontinuity arises in the mechanics path; the jump is
% reported for the route-geometry layer.
ok = true;

% med_gas P2 conditional on knee_angle_r [-45, +10] deg
mga = fr_muscleroute(C.osimStock, 'med_gas_r');
P2 = mga.points{2};
ok = ok && isfield(P2, 'range');
[qAct, qIn, jumpMm] = fr_probeFlags(C, kin, mga, 'knee_angle_r', ...
    P2.range, [-44 -46]);
fprintf(['    med_gas P2 conditional knee [%.1f, %.1f] deg: ' 'active(-44deg)=%d active(-46deg)=%d, |dL| across transition ' '%.2f mm -> %s\n'], rad2deg(P2.range(1)), rad2deg(P2.range(2)), qAct, qIn, jumpMm, fr_tf(qAct ~= qIn));
ok = ok && (qAct ~= qIn);

% psoas P3 conditional on hip_flexion_r (stock model keeps the range)
pso = fr_muscleroute(C.osimStock, 'psoas_r');
found = false;
for i = 1:numel(pso.points)
    if isfield(pso.points{i}, 'range')
        P3 = pso.points{i};
        found = true;
        break
    end
end
if found
    [hAct, hIn, ~] = fr_probeFlags(C, kin, pso, 'hip_flexion_r', ...
        P3.range, [44 46]);
    fprintf(['    psoas %s conditional hip [%.1f, %.1f] deg: ' 'active(+44deg)=%d active(+46deg)=%d -> %s\n'], P3.name, rad2deg(P3.range(1)), rad2deg(P3.range(2)), hAct, hIn, fr_tf(hAct ~= hIn));
    ok = ok && (hAct ~= hIn);
else
    fprintf('    psoas stock route has no conditional point\n');
end

% the built knee case's route-length response (geometry layer)
csK = cases(strcmp({cases.name}, 'knee_bowden'));
csK = csK(1);
q0 = csK.passes{1}.qStructs{1};
q0.knee_angle_r = -46 * pi / 180;
out1 = fr_evalRoute(mga, kin, q0);
q0.knee_angle_r = -44 * pi / 180;
out2 = fr_evalRoute(mga, kin, q0);
fprintf(['    route length across the transition: %.2f mm (inactive) vs ' '%.2f mm (active), jump %.2f mm (quantified; grid mechanics rows ' 'documented fixed)\n'], out1.lengthM * 1000, out2.lengthM * 1000, abs(out2.lengthM - out1.lengthM) * 1000);
end

function [aIn, aOut, jumpMm] = fr_probeFlags(~, kin, route, coord, ...
    range, probeDeg)
qA = kin.qDefault;
qA.(coord) = probeDeg(1) * pi / 180;
qB = kin.qDefault;
qB.(coord) = probeDeg(2) * pi / 180;
oA = fr_evalRoute(route, kin, qA);
oB = fr_evalRoute(route, kin, qB);
aIn = fr_condActive(oA, range);
aOut = fr_condActive(oB, range);
jumpMm = abs(oA.lengthM - oB.lengthM) * 1000;
end

function act = fr_condActive(evalOut, range)
act = false;
for t = 1:numel(evalOut.transitions)
    if abs(evalOut.transitions(t).range(1) - range(1)) < 1e-9 && ...
            abs(evalOut.transitions(t).range(2) - range(2)) < 1e-9
        act = evalOut.transitions(t).active;
    end
end
end

function s = fr_tf(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end
