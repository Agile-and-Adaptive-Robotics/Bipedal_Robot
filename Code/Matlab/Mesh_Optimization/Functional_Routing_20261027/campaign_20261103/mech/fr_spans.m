function S = fr_spans(C, kin, pts, rowBodies, passQList, spanKind)
%FR_SPANS Required minimum tendon spans computed from geometry + hardware.
%
% master directive item 6 / target_contract.json tendonSpanPlaceholders:
% the placeholders are resolved HERE from the route's own feasible
% geometry over the full motion range plus hardware clearance - not by an
% arbitrary token length.
%
% pts: nRows x 3 current route points in their assigned frame per row
% (from fr_mechBuild's location); rowBodies: owning frame per row.
% Required span = max over the pose list of the exit-point-to-insertion
% CHORD distance (the longest span the tendon must accommodate; a shorter
% fixed tendon would jam at deep flexion) plus
%   2 x ferrule (12 mm PLACEHOLDER from the Blender library, unverified)
%   + 1 x housing-stop allowance (20 mm stand-in).
%
% spanKind: 'knee' (knee-crossing central span) or 'foot' (foot route
% minimum tendon span; no BPA hardware on the foot - the span is measured
% from the distal-tibia exit anchor to the foot insertion, which is what
% keeps the sleeve itself off the foot).
switch spanKind
    case 'foot'
        footBodies = {'talus_r', 'calcn_r', 'toes_r'};
        exitIdx = 0;
        for i = 1:numel(rowBodies)
            if ~any(strcmp(rowBodies{i}, footBodies))
                exitIdx = i;   % last non-foot body row = sleeve-side exit
            end
        end
        insIdx = numel(rowBodies);
    case 'knee'
        exitIdx = max(1, numel(rowBodies) - 1);
        insIdx = numel(rowBodies);
    otherwise
        error('fr_spans:Kind', 'unknown span kind %s', spanKind);
end
nG = numel(passQList);
span = zeros(nG, 1);
for i = 1:nG
    q = passQList{i};
    Te = fr_bodyT(kin, q, rowBodies{exitIdx});
    Ti = fr_bodyT(kin, q, rowBodies{insIdx});
    pe = Te(1:3, 1:3) * pts(exitIdx, :)' + Te(1:3, 4);
    pi_ = Ti(1:3, 1:3) * pts(insIdx, :)' + Ti(1:3, 4);
    span(i) = norm(pi_ - pe);
end
allowance = 2 * (C.bowden.ferrule_mm / 1000) + C.bowden.stopAllowance_mm / 1000;
S.spanMax_m = max(span);
S.spanAtPose_m = span;
S.allowance_m = allowance;
S.requiredSpan_m = S.spanMax_m + allowance;
S.kind = spanKind;
S.exitIdx = exitIdx;
S.note = sprintf(['required %s tendon span = max chord (%.1f mm) over ' ...
    'the case motion range + 2 ferrules (12 mm PLACEHOLDER) + stop ' ...
    'allowance (20 mm stand-in) = %.1f mm'], spanKind, ...
    S.spanMax_m * 1000, S.requiredSpan_m * 1000);
end
