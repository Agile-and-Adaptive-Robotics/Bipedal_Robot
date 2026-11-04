function [mech, meta] = fr_mechBuild(C, caseDef, design, kin, xiBranch, qg)
%FR_MECHBUILD Build the verified transmission mechanics for one evaluation.
%
% Mechanics core = BiPam_pulley2_20261022 (verified by
% mechanics_20261103/mech4: closure 1.1e-16 m, force balance 0 N, work
% identity 6e-14 J, chord virtual-work arm error <= 5e-7 m). ONE class
% object per evaluation, on the class's rectangular grid:
%   T = 4x4 x nG x 2, slot 1 (crossing-1 transform) at T(:,:,i,1) for
%   pose i; slot 2 (crossing-2 transform, constant across every
%   demonstrator pass - the swept joint always rides T1) replicated at
%   T(:,:,1,2) and T(:,:,2,2); the class reads gridSize = (nG, 2) and the
%   wrapper reads grid COLUMN 1 (flattened index i) as pose i.
% The multi-pose grid is REQUIRED: the tackle's travel division acts on
% deltaL = span(theta0) - span(theta) WITHIN one object's grid, so isolated
% per-pose objects would see deltaL = 0 and the tackle would divide
% nothing.
%
% xiBranch:
%   'rigid'     - compliance disabled: Xi0 = 0, Xi1 = Xi2 = Inf (the
%                 class's own rigid convention: fortzBi2 skips bracket
%                 deformation when isinf(X1) && isinf(X2); ZERO stiffness
%                 would make C_bracket singular), effectively rigid tendon
%                 (wraps 1e4), no mount/Bowden compliance, Xi3 = 0.
%   'compliant' - adopted Xi0-Xi3 (target_contract.json xiSource; Xi0-Xi3
%                 only, Xi4 waived): bracket stiffness Xi1/Xi2, tendon
%                 stiffness wraps = 6, tackle-mount / Bowden-anchor
%                 cantilever stand-in compliance (6061-T6, [25 5 40] mm,
%                 PENDING-BEN), and the Xi3 wrapping-loss factor coupled
%                 INSIDE the equilibrium by a fixed point on
%                 deltaRoute = Xi3*bendMeasure*comp^2 (the MonoPam_mult
%                 Contraction_k form, MonoPam_mult.m:449) folded through
%                 the class's Xi0 offset channel. STAGE-1 SCOPE: the class
%                 channel is scalar-per-object, so the coupled value is the
%                 max-envelope deltaRoute (conservative); the per-pose
%                 bendMeasure and deltaRoute profiles are stored on mech
%                 for reporting, and the exact per-pose coupling (a
%                 per-cell Xi3 channel in the class) is a recorded later
%                 campaign item.
%
% xi3Mode (caseDef field, OPTIONAL; default 'legacyMaxEnvelope' = the
% stage-1 path EXACTLY as above, kept switchable for A/B comparison):
%   'perPose' - the derived class BiPam_pulley2_xi3pose (campaign folder;
%     derived copy - original untouched) applies Xi3 PER POSE INSIDE the
%     force-length-deflection equilibrium: deltaRoute(i) = Xi3*bend(i)
%     *comp(i)^2 is added to the per-cell contraction demand s0(i) before
%     the class's closure solve (dissertation 30-torque.tex eq:deltaL +
%     eq:x3; MonoPam_mult.m Contraction_k lines 445-449 port; comp from
%     the Xi0-only pre-solve strain, single pass, no fixed point - comp
%     is independent of the solve). Xi0 keeps its own meaning: it is NOT
%     folded with the loss anymore (xi0Eff = the clean Xi0). The bend
%     profile is PER-CONTACT and PER-POSE: deviation arcs at BOTH crossing
%     contacts (correct parent-chain frames; the legacy fr_bendMeasure
%     single-crossing estimate is kept only for the legacy path), summed
%     per cell, wRap = C.wRap (25 mm, buildGeoExclusion lineage). The
%     compliant branch also enables the dissertation's INSERTION bracket
%     (eq:bktstiffness K = diag(chi1, chi2, chi1), two-rotation frame
%     port) as a separate, once-counted compliance alongside the origin
%     bracket [chi1 chi2 chi2], the tackle-mount stand-in, and the Bowden
%     housing stops (C.Xi.noInsBracket probe suppresses it).
%     mech.xi3PerPoseLoss_m / mech.deltaRoute carry the applied per-pose
%     loss, mech.bendArcPerContact the per-contact arcs.
%
% qg: struct with .passIdx and .qList (1xN cell of full pose structs).

pass = caseDef.passes{qg.passIdx};
qList = qg.qList;
nG = numel(qList);

% Xi3 mode: explicit A/B switch (default = stage-1 legacy path, so a
% caseDef without the field behaves byte-identically to before).
if isfield(caseDef, 'xi3Mode') && ~isempty(caseDef.xi3Mode)
    xi3Mode = caseDef.xi3Mode;
else
    xi3Mode = 'legacyMaxEnvelope';
end
if ~any(strcmp(xi3Mode, {'legacyMaxEnvelope', 'perPose'}))
    error('fr_mechBuild:Xi3Mode', ...
        ['caseDef.xi3Mode must be ''legacyMaxEnvelope'' or ' ...
         '''perPose'', got %s.'], mat2str(xi3Mode));
end

% ---- 1. location rows in their class-convention frames ------------------
[loc, rowBodies] = fr_locForPass(C, caseDef, design, kin, qg.passIdx);
caseDef = fr_caseForPass(caseDef, qg.passIdx);   % per-pass cross/rows
caseDef.rowBodiesUsed = rowBodies;

% ---- 2. sizing ------------------------------------------------------------
rest = design.L0;
kmax = campaignLminFraction(C, caseDef.diaStr) * rest;

% ---- 3. discrete transmission configuration (kept distinct) -------------
cfg = caseDef.txConfig;
% physical sheave count from the reeving gain (never inferred from travel):
% a G-line tackle block carries (G-1) moving sheaves plus one fixed
% redirect when G > 1. Recorded per crossing.
cfg.sheaveCountPhysical = tackleSheaves(cfg.gain(1)) + tackleSheaves(cfg.gain(2));
cfg.seriesSections = 1 + double(design.seriesSplit < 0.999);
cfg.sectionLengths_m = [rest * design.seriesSplit, ...
    rest * (1 - design.seriesSplit)];

% ---- 4. Xi branch constants ----------------------------------------------
switch xiBranch
    case 'rigid'
        % the class's rigid convention is INFINITE bracket stiffness
        xi0 = 0; xi1 = Inf; xi2 = Inf; wraps = 1e4;
        cfgMount = struct();
        xi3On = false;
    case 'compliant'
        xi0 = C.Xi.Xi0_flx;      % adopted flexor pick 77 (contract xiSource)
        xi1 = C.Xi.Xi1;
        xi2 = C.Xi.Xi2;
        wraps = 6;
        if isfield(C.Xi, 'wrapsOverride')
            wraps = C.Xi.wrapsOverride;   % rigid-limit recovery probe
        end
        cfgMount = struct();
        if caseDef.hasTackle
            cfgMount.mountMaterial = 'aluminum6061T6';
            cfgMount.mountSection = [0.025 0.005 0.040]; % PENDING-BEN
        end
        if any(strcmp(cfg.routingMode, 'bowden'))
            cfgMount.bowdenBracketMaterial = 'aluminum6061T6';
            cfgMount.bowdenBracketSection = [0.025 0.005 0.040];
        end
        xi3On = true;
        if strcmp(xi3Mode, 'perPose') ...
                && ~(isfield(C.Xi, 'noInsBracket') && C.Xi.noInsBracket)
            % dissertation two-bracket mechanics: the INSERTION bracket
            % (route end, tendon anchor) with the eq:bktstiffness
            % assignment diag(chi1, chi2, chi1) and the two-rotation
            % 'auto' (load-path) reference; separate hardware from the
            % tackle mount and Bowden stops, counted once (the derived
            % class carries it as its own closure coefficient).
            cfgMount.insertionBracketStiffness = [xi1, xi2, xi1];
        end
        if isfield(C.Xi, 'noCompliance') && C.Xi.noCompliance
            % rigid-limit recovery probe: suppress the tackle-mount /
            % Bowden-anchor cantilever compliance too
            cfgMount = struct();
        end
    otherwise
        error('fr_mechBuild:Branch', 'unknown xiBranch %s', xiBranch);
end

% ---- 5. transform grid: slot 1 per pose, slot 2 replicated ---------------
T = zeros(4, 4, nG, 2);
T1list = cell(1, nG);
T2list = cell(1, nG);
for i = 1:nG
    T(:,:,i, 1) = pass.T1(kin, qList{i});
    T1list{i} = T(:,:,i, 1);
    T2list{i} = pass.T2(kin, qList{i});
end
T2ref = T2list{ceil(nG / 2)};
T(:,:,1, 2) = T2ref;
if nG >= 2
    T(:,:,2, 2) = T2ref;
end

% ---- 6. construct + Xi3 coupling (explicit A/B switch) --------------------
fitn = caseFitting(C, caseDef);
deltaRoute = 0;
bendM = zeros(nG, 1);
bendCols = zeros(nG, 0);          % per-contact arcs (perPose only)
xi0Eff = xi0;
it = 1;
switch xi3Mode
    case 'legacyMaxEnvelope'
        % stage-1 path, verbatim: max-envelope deltaRoute through the
        % class's Xi0 channel (documented gap; kept for A/B comparison
        % and for the running campaign's resume compatibility).
        for it = 1:3
            xi0Eff = xi0 + deltaRoute;
            m = BiPam_pulley2_20261022(caseDef.name, loc, caseDef.cross, ...
                caseDef.diaNum, T, rest, kmax, design.tendon, fitn, ...
                C.pressure, xi0Eff, xi1, xi2, wraps, ...
                mergeCfg(cfg, cfgMount));
            if ~xi3On
                break
            end
            bendM = fr_bendMeasure(C, m, caseDef, loc, T1list, T2list);
            KMAXfrac = (rest - kmax) / rest;
            strainCol = m.strain_p(1:nG);          % grid column 1 = the poses
            comp = max(0, 1 - strainCol / KMAXfrac);
            deltaNew = max(C.Xi.Xi3 * bendM .* comp.^2);
            if it > 1 && abs(deltaNew - deltaRoute) < 1e-9
                deltaRoute = deltaNew;
                break
            end
            deltaRoute = deltaNew;
        end
    case 'perPose'
        % per-pose per-contact bend profile -> derived class; the loss is
        % applied PER CELL inside the equilibrium (see header), Xi0 stays
        % clean, single deterministic pass (comp from the Xi0-only
        % pre-solve strain; MonoPam_mult.m:445-449 port).
        if cfg.gain(2) > 1
            activeKp = 2;
        else
            activeKp = 1;
        end
        bendCols = fr_bendProfile(C, caseDef, loc, T1list, T2list, activeKp);
        bendGrid = zeros(nG, 2);
        bendGrid(:, 1) = sum(bendCols, 2);   % per-pose total over contacts
        bendGrid(:, 2) = bendGrid(ceil(nG / 2), 1);   % T2ref replica slot
        if xi3On
            xi3val = C.Xi.Xi3;
            if isfield(C.Xi, 'xi3Override') && ~isempty(C.Xi.xi3Override)
                xi3val = C.Xi.xi3Override;   % probe knob (fr_xi3_check)
            end
        else
            xi3val = 0;   % rigid branch: bit-identical base arithmetic
        end
        m = BiPam_pulley2_xi3pose(caseDef.name, loc, caseDef.cross, ...
            caseDef.diaNum, T, rest, kmax, design.tendon, fitn, ...
            C.pressure, xi0, xi1, xi2, wraps, mergeCfg(cfg, cfgMount), ...
            xi3val, bendGrid);
        deltaRoute = m.deltaRoute(1:nG);         % per-pose applied loss
        bendM = bendGrid(:, 1);
end

% ---- 7. wrapper fields (grid column 1 = the poses) ------------------------
mech = struct();
mech.caseDef = caseDef;
mech.xiBranch = xiBranch;
mech.loc = loc;
mech.rowBodies = rowBodies;
mech.rest = rest;
mech.kmax = kmax;
mech.cfg = mergeCfg(cfg, cfgMount);
mech.qg = qg;
mech.cell = m;
mech.nG = nG;
mech.strain_p = m.strain_p(1:nG);
mech.FsingleBPA = m.FsingleBPA(:, 1);
mech.Ftendon = m.Ftendon(:, 1);
mech.segTension = m.segTension(1:nG, :);
mech.deltaL = squeeze(m.deltaL(:, 1, :));
mech.PulleyInfeasible = m.PulleyInfeasible(:, 1);
mech.PulleySlack = m.PulleySlack(:, 1);
mech.PulleyTravel = m.PulleyTravel(:, 1);
mech.gama = m.gama(1:nG);
mech.ReactionFmag = m.ReactionFmag(:, 1);
mech.T1list = T1list;
mech.T2list = T2list;
mech.xi0Eff = xi0Eff;
mech.bendMeasure = bendM;
mech.xi3MaxDeltaRoute_m = max(deltaRoute);
mech.xi3Iterations = it;
mech.xi3Mode = xi3Mode;
mech.deltaRoute = deltaRoute;            % scalar (legacy max envelope) or
                                         % per-pose vector (perPose)
mech.bendArcPerContact = bendCols;       % nG x nContacts (perPose; 0 cols legacy)
if strcmp(xi3Mode, 'perPose')
    mech.xi3PerPoseLoss_m = m.deltaRoute(1:nG);  % per-pose applied loss
    mech.cbBracket = m.cbBracket;        % origin-bracket force-path projection
    mech.InsBracketDeflection = m.InsBracketDeflection;
    mech.InsBracketReaction = m.InsBracketReaction;
else
    mech.xi3PerPoseLoss_m = zeros(nG, 1);% legacy folds the loss into xi0Eff
    mech.cbBracket = [];
    mech.InsBracketDeflection = [];
    mech.InsBracketReaction = [];
end
mech.ActiveCrossing = m.ActiveCrossing;
mech.kSpr = m.kSpr;

meta = struct();
meta.loc = loc;
meta.rowBodies = rowBodies;
meta.T1list = T1list;
meta.T2list = T2list;
meta.rest = rest;
meta.kmax = kmax;
meta.cfg = mech.cfg;
meta.qg = qg;
meta.passIdx = qg.passIdx;
meta.caseDef = caseDef;
meta.xiBranch = xiBranch;
meta.xi0Eff = xi0Eff;
meta.xi3Mode = xi3Mode;
end

function c2 = fr_caseForPass(c2, passIdx)
if numel(c2.crossFor) >= passIdx
    c2.cross = c2.crossFor{passIdx};
end
mask = c2.rowUse{passIdx};
c2.nRows = sum(mask);
end

function s = tackleSheaves(g)
if g > 1
    s = (g - 1) + 1;   % moving sheaves + fixed redirect
else
    s = 0;
end
end

function c = campaignLminFraction(C, diaStr)
if strcmp(diaStr, '10')
    c = C.LminFraction.d10;
else
    c = C.LminFraction.d20;
end
end

function b = fr_bendMeasure(C, m, caseDef, loc, T1list, T2list)
% wrap-arc estimate per pose at the ACTIVE crossing:
% bendMeasure = wRap * deviation angle between the segment into the
% crossing and the segment out of it (buildGeoExclusion.m:51-53 lineage:
% "wRap*wrap angle is an arc-length estimate used only for Xi3 bend loss").
k = m.ActiveCrossing;
nG = numel(T1list);
bend = zeros(nG, 1);
rIn = caseDef.cross(k) - 1;
rOut = min(caseDef.cross(k) + 1, caseDef.nRows);
for i = 1:nG
    T1 = T1list{i};
    T2 = T2list{i};
    pIn = fr_applyT(T1, loc(rIn, :));
    pX = fr_applyT(T2 * T1, loc(caseDef.cross(k), :));
    pOut = fr_applyT(T2 * T1, loc(rOut, :));
    u1 = pX - pIn; u2 = pOut - pX;
    if norm(u1) < 1e-9 || norm(u2) < 1e-9
        bend(i) = 0;
    else
        cs = dot(u1 / norm(u1), u2 / norm(u2));
        bend(i) = C.wRap * acos(min(1, max(-1, cs)));
    end
end
b = bend;
end

function B = fr_bendProfile(C, caseDef, loc, T1list, T2list, activeK)
% Per-pose per-contact tendon bend arcs (meters): wRap * |turn angle| at
% EACH tendon-side sheave/guide contact of the routed polyline, per pose.
% Contact set (tendon-bend-only semantics - sleeves are straight pinned
% pressure vessels, only the TENDON bends over the tackle sheave and any
% downstream guides): the interior route vertices r with
% exitK <= r <= nRows-1, exitK = cross(activeK) - 1 (the tackle exit;
% the class's PulleyExitIndex default), i.e. every interior vertex
% at/after the exit lies on the tendon span exit->insertion. The turn
% angle at vertex r is measured between the segment into r and the
% segment out of it, expressed in vertex r's OWN class frame via the
% consistent parent chain - the IDENTICAL construction of the derived
% class's selfBendArcs (fr_xi3_check asserts wrapper == class). buildGeo
% Exclusion.m geo.wRap lineage: arc-length ESTIMATE used only for the Xi3
% bend loss, never added to the geometric centerline. Degenerate
% segments contribute 0. Returns nG x nContacts (per-pose rows).
crossV = caseDef.cross;
nU = size(loc, 1);
nG = numel(T1list);
exitK = crossV(activeK) - 1;
% a bend CONTACT needs both neighbors on the route (incoming and outgoing
% segments): rows r >= 2 and r <= nU-1. Short routes (e.g. the 2-point
% trunk: exit row 1, no interior tendon vertex) yield zero contacts and a
% zero-column profile - the Xi3 loss is then 0 by construction, which is
% the correct tendon-bend-only semantics for a straight 2-anchor tendon.
contacts = max(2, exitK):(nU - 1);
nC = max(numel(contacts), 0);
B = zeros(nG, nC);
for ci = 1:nC
    r = contacts(ci);
    fr = classFrameOf(r, crossV);
    for i = 1:nG
        T1 = T1list{i};
        T2 = T2list{i};
        pIn = rowToFrame(loc(r - 1, :), classFrameOf(r - 1, crossV), ...
            fr, T1, T2);
        pX = rowToFrame(loc(r, :), fr, fr, T1, T2);
        pOut = rowToFrame(loc(r + 1, :), classFrameOf(r + 1, crossV), ...
            fr, T1, T2);
        u1 = pX - pIn;
        u2 = pOut - pX;
        if norm(u1) < 1e-9 || norm(u2) < 1e-9
            B(i, ci) = 0;
        else
            cs = dot(u1 / norm(u1), u2 / norm(u2));
            B(i, ci) = C.wRap * acos(min(1, max(-1, cs)));
        end
    end
end
end

function f = classFrameOf(u, crossV)
% Class frame index (1/2/3) of per-pass used row u - BiPam rowFrame
% semantics: rows before Cross(1) proximal; duplicate crossing
% (Cross(1) == Cross(2)) has no middle frame; rows at/beyond Cross(2)
% distal.
c1 = crossV(1);
c2 = crossV(2);
if u < c1
    f = 1;
elseif c1 == c2
    f = 3;
elseif u < c2
    f = 2;
else
    f = 3;
end
end

function p = rowToFrame(p, fFrom, fTo, T1, T2)
% Express a point given in class frame fFrom in class frame fTo, chaining
% the per-pose transforms (T1 maps middle->proximal, T2 maps
% distal->middle; fr_chainT direction, gate-verified). Both directions
% handled (the bend profile maps neighboring rows ONTO the vertex frame).
while fFrom < fTo
    if fFrom == 1
        p = fr_applyT(inv(T1), p);   % proximal -> middle
        fFrom = 2;
    else
        p = fr_applyT(inv(T2), p);   % middle -> distal
        fFrom = 3;
    end
end
while fFrom > fTo
    if fFrom == 3
        p = fr_applyT(T2, p);        % distal -> middle
        fFrom = 2;
    else
        p = fr_applyT(T1, p);        % middle -> proximal
        fFrom = 1;
    end
end
end

function f = caseFitting(C, caseDef)
if isfield(caseDef, 'fittingLength') && ~isempty(caseDef.fittingLength)
    f = caseDef.fittingLength;
else
    f = C.fittingLength;
end
end

function p = fr_applyT(Tm, prow)
p = (Tm(1:3, 1:3) * prow(:) + Tm(1:3, 4))';
end

function cfg = mergeCfg(a, b)
cfg = a;
fn = fieldnames(b);
for i = 1:numel(fn)
    cfg.(fn{i}) = b.(fn{i});
end
end
