function env = fr_envelope(C, caseDef, design, kin, qPoseList, xiBranch)
%FR_ENVELOPE Pose-dependent occupied envelope of one optimized group (v2).
%
% v2 (2026-11-04, Ben's corrected layering ruling): primitives are TAGGED
% so the constraint layer can split SOFT from HARD:
%   - capsules (BPA sleeve branches)        -> kind 'bpa_sleeve'  (SOFT:
%     slight BPA-BPA interference is a soft objective penalty with depth
%     and extent recorded; never a hard clearance constraint)
%   - spheres 'fitting_origin/fitting_exit' -> kind 'rigid_fitting' (HARD:
%     rigid metal hardware)
%   - sphere  'exit_row' (sleeve end bulk)  -> kind 'bpa_end'      (SOFT)
%   - sphere  'tackle_*'                    -> kind 'rigid_tackle' (HARD)
%
% The envelope is still built from the SAME solved mechanics the objective
% uses (fr_mechBuild -> the bi transmission class), evaluated at the poses
% of the group that must avoid it (physically simultaneous configurations,
% handoff item 9). Ben's layering ruling (master item 5) is retained:
% finished layers are carried into outer layers - the RIGID parts as hard
% obstacles, the BPA sleeve corridors through the soft penalty.
env.caseName = caseDef.name;
env.groupMapId = caseDef.groupMapId;
env.design = design;
env.xiBranch = xiBranch;
nG = numel(qPoseList);
qg = struct('passIdx', caseDef.envelopePassIdx, 'qList', {qPoseList});
[mech, meta] = fr_mechBuild(C, caseDef, design, kin, xiBranch, qg);
loc = meta.loc;
Drest = fr_diameter(C, caseDef.diaStr, 0);
% the SLEEVE (BPA side) spans rows 1..cross(1)-1: the capsule ends at the
% tendon EXIT row (= the last pre-crossing row), not at the first distal
% row (that would draw the envelope across the joint, which the sleeve
% never occupies - the tendon does, and tendons are not envelope bodies)
env.exitRow = caseDef.cross(1) - 1;
env.strainAtPose = mech.strain_p;
env.capsules = cell(1, nG);
env.spheres = cell(1, nG);
crossV = caseDef.crossFor{caseDef.envelopePassIdx};
for i = 1:nG
    q = qPoseList{i};
    Tw = fr_bodyT(kin, q, meta.rowBodies{1});   % prox frame -> world
    T1 = meta.T1list{i};
    T2 = meta.T2list{i};
    % per-row-frame world transform: rows < cross1 proximal (no T1), rows
    % >= cross1 transformed by T2*T1 (same rule the class uses)
    if 1 < crossV(1)
        p1 = fr_applyT(Tw, loc(1, :));
    else
        p1 = fr_applyT(Tw * T2 * T1, loc(1, :));
    end
    if env.exitRow < crossV(1)
        p2 = fr_applyT(Tw, loc(env.exitRow, :));
    else
        p2 = fr_applyT(Tw * T2 * T1, loc(env.exitRow, :));
    end
    Dinfl = fr_diameter(C, caseDef.diaStr, mech.strain_p(i));
    r = (Dinfl / 2) / 1000;   % m
    u = (p2 - p1) / max(norm(p2 - p1), 1e-9);
    u = u(:)';   % row 1x3 for the perp basis below
    if any(~isfinite(u))
        u = [1; 0; 0]';
    end
    uc = u(:);   % 3x1 column
    perp = null(uc');
    perp = perp(:, 1:min(2, size(perp, 2)));   % 3x1 or 3x2 columns
    if size(perp, 2) < 2 || any(~isfinite(perp(:)))
        % degenerate axis: fixed perpendicular frame, columns 3x1
        alt = [1; 0; 0];
        if abs(dot(alt, uc)) > 0.9
            alt = [0; 1; 0];
        end
        perp = [cross(alt, uc), cross(uc, alt)];   % 3x2, both columns perp
    end
    offs = fr_bundleOffsets(caseDef.seed.nBPA) * (Drest / 2) / 1000;
    ncap = max(1, size(offs, 1));
    env.capsules{i} = repmat(struct('p1', [], 'p2', [], 'r', 0, ...
        'branch', 0, 'kind', 'bpa_sleeve'), 1, ncap);
    for b = 1:ncap
        if b <= size(offs, 1)
            off = perp * offs(b, :)';
        else
            off = [0; 0; 0];
        end
        env.capsules{i}(b).p1 = p1 + off';
        env.capsules{i}(b).p2 = p2 + off';
        env.capsules{i}(b).r = r;
        env.capsules{i}(b).branch = b;
        env.capsules{i}(b).kind = 'bpa_sleeve';
    end
    % field-by-field (struct('p', [3x3], 'what', {cell}) would broadcast
    % into a struct array - the struct-with-array-value trap)
    sph = struct();
    sph.p = [p1; p2; p2];
    sph.r = [r, r, r];
    sph.what = {'fitting_origin', 'fitting_exit', 'exit_row'};
    sph.kind = {'rigid_fitting', 'rigid_fitting', 'bpa_end'};
    if caseDef.hasTackle
        fam = C.families(caseDef.txConfig.tackleFamily);
        sph.p(end+1, :) = p2;
        sph.r(end+1) = (fam.od_mm / 2 + 2) / 1000;
        sph.what{end+1} = ['tackle_' caseDef.txConfig.tackleFamily];
        sph.kind{end+1} = 'rigid_tackle';
    end
    env.spheres{i} = sph;
end
end

function p = fr_applyT(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end

function offs = fr_bundleOffsets(nBPA)
switch min(nBPA, 3)
    case 1
        offs = [0 0];
    case 2
        offs = [0 1; 0 -1];
    otherwise
        th = (0:2) * 2 * pi / 3;
        offs = [cos(th)', sin(th)'];
end
end
