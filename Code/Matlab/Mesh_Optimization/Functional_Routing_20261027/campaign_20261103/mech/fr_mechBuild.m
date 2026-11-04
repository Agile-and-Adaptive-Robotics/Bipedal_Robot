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
% qg: struct with .passIdx and .qList (1xN cell of full pose structs).

pass = caseDef.passes{qg.passIdx};
qList = qg.qList;
nG = numel(qList);

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

% ---- 6. construct + Xi3 fixed point ---------------------------------------
fitn = caseFitting(C, caseDef);
deltaRoute = 0;
bendM = zeros(nG, 1);
for it = 1:3
    xi0Eff = xi0 + deltaRoute;
    m = BiPam_pulley2_20261022(caseDef.name, loc, caseDef.cross, ...
        caseDef.diaNum, T, rest, kmax, design.tendon, fitn, C.pressure, ...
        xi0Eff, xi1, xi2, wraps, mergeCfg(cfg, cfgMount));
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
mech.xi3MaxDeltaRoute_m = deltaRoute;
mech.xi3Iterations = it;
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
