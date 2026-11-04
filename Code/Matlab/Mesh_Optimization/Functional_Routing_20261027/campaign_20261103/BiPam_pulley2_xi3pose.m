% BiPam_pulley2_xi3pose -- PER-POSE Xi3 WRAPPING-LOSS DERIVED COPY of
% BiPam_pulley2_20261022 for the Functional_Routing_20261027 campaign
% (campaign_20261103). DERIVED COPY - ORIGINAL UNTOUCHED: Ben's
% Code/Matlab/Mesh_Optimization/BiPam_pulley2_20261022.m is not modified;
% this file carries the Xi3-per-pose campaign item recorded as the stage-1
% gap in mech/fr_mechBuild.m's header.
%
% >>> WHY (MASTER_DIRECTIVE_20261103.md item 4 / HARDWARE_AND_ZCODE_HANDOFF
% items 8-9: Xi0-Xi3 only, Xi4 waived) <<<
% The stage-1 wrapper coupled Xi3 through this class's Xi0 channel as a
% MAX-ENVELOPE scalar (one deltaRoute per object, max over poses), with the
% per-pose bendMeasure/deltaRoute profiles stored on mech for reporting
% only. The dissertation defines chi3 PER CONTACT ARC and PER POSE:
%   30-torque.tex eq:deltaL (Actuator Bending, sec:fictive_contraction):
%     Delta_l = chi3 * R * |theta_wrap - theta_k| * zeta^2   (when bent),
%     zeta = 1 - epsilon* (eq:comp), and
%     l_m = l_LMT - 2*l_fitting - chi0 - Delta_l   (eq:x3);
%   94-AppendixC.tex (tab:app_stiffness_arrays provenance paragraph):
%     "The wrapping-loss term X3 ... is evaluated as X3 R |Delta theta|
%     zeta^2 PER CONTACT ARC ... on the biomimetic configurations the bend
%     measure is w_rap|theta_turn| with w_rap = 25 mm."
%   chi3 is DIMENSIONLESS (App C tab:app_opt_bounds: chi3 = x3 in [0,1]);
%   R|Delta theta| carries the meters (an arc-length ESTIMATE:
%   buildGeoExclusion.m geo.wRap comment "wRap*wrap angle is an arc-length
%   estimate used only for Xi3 bend loss. It is not a torus arc added to
%   the Location centerline or MuscleLength"). The campaign models sleeves
%   as straight pinned pressure vessels, so ONLY THE TENDON bends over the
%   crossing guides/sheaves: the bend arcs are evaluated at the crossing
%   contacts of the routed polyline, never added to the geometric
%   centerline. The validated MonoPam implementation of record is
%   Robot_Data/MonoPam_mult.m Contraction_k (lines 435-462):
%     deltaRoute = X3 .* bendMeasure .* comp.^2   (line 449),
%     comp = max(0, 1 - strain0/KMAX) computed from the Xi0(+gama)-only
%     pre-solve strain (lines 445-447; gama = 0 pre-solve, line 290-291),
%     and Lm_adj = Lmt - tendon - 2*fit - X0 - gama - deltaRoute (line 459)
%     feeds the per-pose force balance inside fortz (line 542:
%     mL{j} = mL{j} - deltaL{j} - ...), i.e. the loss acts INSIDE the
%     force-length equilibrium as a per-pose usable-length reduction.
%
% >>> WHAT THIS DERIVED CLASS CHANGES vs BiPam_pulley2_20261022 <<<
% (1) CONSTRUCTOR: two OPTIONAL trailing arguments
%       PD = BiPam_pulley2_xi3pose(..., pulleyConfig, xi3, bendMeasure)
%     arg 16 xi3 (scalar, dimensionless, default 0) and arg 17 bendMeasure
%     (per-cell bend arc profile, default [] = self-compute). All original
%     arguments keep their order and meaning (Kmax remains the 7th).
% (2) PER-CELL Xi3 INSIDE THE EQUILIBRIUM (updateBiPulleyGeometry): after
%     the base s0 flatten, ONE new step computes the per-cell loss
%       deltaRoute(i) = Xi3 * bend_i * comp(i)^2,
%     comp from the Xi0-only pre-solve strain strainXi0f (MonoPam port;
%     single pass, deterministic - comp is mathematically independent of
%     the solve, so no fixed point is needed or possible), and adds it to
%     the per-cell contraction demand s0f BEFORE LokBi2/fortzBi2. This is
%     exactly the eq:x3 semantics (Delta_l subtracted from usable l_m
%     inside the equilibrium) but PER GRID CELL instead of per object.
%     Xi0 keeps its own meaning: it stays the constant per-object offset
%     term; the loss enters as a SEPARATE additive per-cell term. With
%     Xi3 = 0 the added step is an exact no-op (s0f + 0) and the class
%     executes the base arithmetic bit-for-bit.
% (3) BEND PROFILES: bendMeasure given as an N-vector (flattened grid,
%     ii fastest: lin = ii + (iii-1)*N1) or N1xN2 grid, in METERS
%     (wRap * turn angle per contact, summed over contacts); when empty,
%     the class self-computes the per-contact deviation arcs at the
%     TENDON-SIDE sheave/guide contacts (the interior route vertices
%     exit..npts-1 of the ACTIVE crossing - tackle exit plus downstream
%     guides; sleeves are straight pinned vessels so sleeve-side vias do
%     not bend the tendon) from the rigid geometry via the consistent
%     parent chain, bend = WrapRadius * |turn angle| per contact, 0 at
%     degenerate segments. Stored per contact on obj.BendArc
%     (N x nContacts) for audits. The campaign wrapper passes the
%     equivalent profile; fr_xi3_check asserts both agree.
% (4) INSERTION-BRACKET COMPLIANCE CHANNEL (two-bracket mechanics): the
%     base class carries ONE principal bracket (route row 1, K =
%     [X1 X2 X2], the dissertation's "Biomimetic flexor" row of
%     tab:app_stiffness_arrays). The dissertation identifies TWO principal
%     brackets (ch.30: origin bracket K = [chi1 chi1 chi2], insertion
%     bracket K = diag(chi1, chi2, chi1) = eq:bktstiffness; both
%     two-rotation yaw-pitch frames, eq:app_bracketframe). The derived
%     class adds the INSERTION bracket as a NEW, SEPARATE compliance:
%       config 'insertionBracketStiffness' (3-vector N/m; default the
%         dissertation insertion assignment [chi1 chi2 chi1] when the
%         caller enables it; zeros(3) = OFF = base-identical),
%       or 'insertionBracketCompliance' (explicit 3x3 m/N, wins),
%       'insertionBracketRefPoint' ('auto' default = the per-cell exit
%         row point, so the frame's x follows the load path per the
%         dissertation's "x axis follows the nominal load path"; or a
%         fixed 1x3 m point in the distal frame).
%     It is projected like the origin bracket (per-cell yaw-pitch bracket
%     frame port, force transformed into the frame, compliance projected
%     on the force path) and enters the closure ONCE as a span-growth
%     give at the tendon-line tension scale nBPA/G:
%       cIns(i) = (nBPA/G) * u_t' * C_ins * u_t   (pre-solve rigid u_t),
%     added to cAdd (positive sign = demand-side give, same role as the
%     mount cmT / bowden cmBw terms). Post-solve it is applied
%     geometrically (last route row shifted by C_ins*(F_t*u_t), rotation-
%     only mapping; row nPts already lives in the distal frame), stored
%     on InsBracketDeflection / InsBracketReaction. Counted ONCE: it is
%     distinct hardware from the tackle-block mount (Cmount) and the
%     Bowden housing stops (Cbowden), each of which keeps its own
%     compliance contribution and reaction term; DeflClosure_cIns is
%     stored so the audit can verify cAdd = -G*nBPA*cmB + nBPA*cmT
%     + cmBw(1) + cmBw(2) + cIns with no term counted twice.
% (5) ORIGIN-BRACKET STIFFNESS ASSIGNMENT made explicit: config
%     'bracketStiffness1' (3-vector N/m) overrides the base's hard-coded
%     K = [X1 X2 X2] (fortzBi2). Default = [X1 X2 X2] so the arithmetic is
%     unchanged; the assignment provenance (App C "Biomimetic flexor"
%     row) is now stated instead of implicit. Isotropic Xi1 is NOT
%     introduced anywhere: the tackle-mount and Bowden housing-stop
%     compliances keep their section-based cantilever model (no
%     dissertation assignment exists for them - PENDING-BEN stand-in),
%     which is the "isotropic Xi1 only where no directional assignment is
%     prescribed" rule satisfied by keeping the hardware stand-ins
%     section-derived and documented, not by pretending anisotropy.
% (6) NEW OUTPUTS (all additive; base outputs unchanged):
%       deltaRoute        N x 1 per-cell applied wrapping loss (m)
%       BendArc           N x 2 per-contact self-computed bend arcs (m)
%       cbBracket         N x 1 origin-bracket force-path projection cb
%       Bracket1Deflection / Bracket1Reaction   N x 3 bracket frame
%       InsBracketDeflection / InsBracketReaction N x 3 distal frame
%       DeflClosure_cIns  N x 1 insertion-bracket closure coefficient
%     strain_p now adds deltaRoute (the fictive loss the force law saw,
%     MonoPam_mult strain_f analog); with Xi3 = 0 this subtracts exact 0
%     and is bit-identical.
% (7) The legacy LengthCheck getter (hard-coded 0.25 contraction, ~line
%     552 of the original) is carried VERBATIM and remains unconsumed on
%     the BiPulley path (the campaign sizes with its own LminFraction
%     limits); not patched, flagged here for honesty.
% (8) DERIVED-CLASS FIX, LokBi2 rigid escape (isinf(X1) && isinf(X2) &&
%     isinf(kSpr)): the base returns slack/infeasible as DOUBLE zeros,
%     which crashes the downstream logical indexing
%     (obj.Ftendon(obj.PulleySlack) = 0) whenever the tendon rate is
%     infinite (wraps = Inf - never used by the campaign, which passes
%     wraps = 1e4, but required by fr_xi3_check's fully-rigid probes).
%     The derived class returns logicals; every logical use is identical.
% GENERALIZATIONS (stated): the dissertation's insertion-bracket stiffness
%     assignment and two-rotation frame were identified on the pinned-knee
%     flexor; applying them to the campaign's per-route insertion brackets
%     is a generalization of record (same pair chi1/chi2, insertion
%     assignment [chi1 chi2 chi1], auto reference = load-path alignment).
%     The per-contact arc set is generalized from the dissertation's
%     fixed knee-wrap radii to the routed polyline's crossing contacts
%     with the campaign's wRap = 25 mm (buildGeoExclusion.m geo.wRap; App
%     C w_rap = 25 mm biomimetic lineage).
%
% >>> END DERIVED-CLASS ADDITIONS; everything below this block is the
% original header of the dated variant, kept verbatim. <<<
%
% Pam Data -- BIARTICULAR reverse-pulley (block-and-tackle) BPA transmission
% VARIANT 2 (dated copy of BiPam_pulley.m, 2026-10-22 -- DO NOT EDIT BEN'S
% BiPam_pulley.m; this file carries the transmission-mechics upgrades of
% PULLEY_SWEEP_DIRECTIVE_20261021.md sections 2 and 3; plan items 3/4).
% Author: Ben Bolen (base class); variant 2026-10-22 per Ben's directive.
% Date: 2026-10-22
%
% >>> WHAT VARIANT 2 ADDS (everything else is verbatim BiPam_pulley) <<<
%
% (a) MOVING-VIA TACKLE-BLOCK MOUNT DEFLECTION (directive section 2).
%   The base class reports the mount reaction at the moving-via tackle
%   block, R = nBPA*F*u_bpa + F_t*u_t, as an OUTPUT only -- no compliance
%   anywhere consumes it. Variant 2 gives that block mount a 3x3 compliance
%   Cmount (config fields mountCompliance, or mountMaterial + mountSection)
%   expressed in the ACTIVE crossing's child frame, and carries the
%   deflection THROUGH THE EXISTING EQUILIBRIUM: with w = u_bpa + u_t/G
%   (the R direction per unit single-BPA force, u_bpa taken from the BPA
%   muscle origin exactly as the row-1 bracket does -- Funit is built from
%   Location(1,:) in computeForceVectorBi), the block deflection
%   delta = Cmount * (nBPA*F*w) contributes two travel-absorption terms,
%   precisely the role cb*F plays for the row-1 (Xi1/Xi2) bracket:
%     d_bpa = nBPA*F * (u_bpa' * Cmount * w)   block gives along the BPA
%                                               line -> BPA-side travel
%                                               absorbed: G*(... - d_bpa)
%     d_tdn = nBPA*F * (u_t'   * Cmount * w)   block gives along the
%                                               tendon line -> span grows:
%                                               demand side (... - d_tdn)
%   so the solved closure becomes (fortzBi2, cAdd = -G*nBPA*cmB + nBPA*cmT):
%     G*(s - s0ref - cb*F - d_bpa) = DeltaL + nBPA*F*cSpr/G - d_tdn
%   Both terms are LINEAR in F(s), so with Cmount = 0 (default) cAdd is
%   EXACTLY zero and psiFun is bit-for-bit the base psiFun -- the recorded
%   batch (nBPA = 1, G = [1 1], moving_via) reproduces at 0.000 N*m.
%   After the solve the deflection is also applied GEOMETRICALLY: the exit
%   row of L_p is shifted by delta (rotation-only mapping into the row's
%   own frame), so u_t, mA_ins, and the deformed segment lengths honor the
%   deflected moving via (first-order: closure directions from the rigid
%   geometry, exactly like cb uses the pre-force direction). Slack cells
%   deflect under the BPA-side term only (Ftendon zeroed at slack, the
%   base ReactionF convention); infeasible cells keep their deflection
%   (their torque is NaN'd downstream regardless).
%
% (b) BOWDEN ANCHOR-BRACKET DEFLECTION AT BICYCLE-GRADE STIFFNESS
%   (directive section 3). The bowden housing anchor bracket (the distal
%   stop the housing reacts against) gets a 3x3 compliance from a
%   cantilever section model, bowdenBracketCompliance or
%   bowdenBracketMaterial + bowdenBracketSection:
%       k_a  = E*A/L            (axial, A = b*h)
%       k_by = 3*E*(b*h^3/12)/L^3,  k_bz = 3*E*(h*b^3/12)/L^3
%       C = diag(1/k_a, 1/k_by, 1/k_bz)
%   with bicycle brake/shifter-grade materials 6061-T6 aluminum
%   (E = 68.9 GPa) or 4130 CroMoly (E = 205 GPa). The anchor deflects
%   under the cable tension it reacts (crossing k's own segment tension:
%   nBPA*F/G at/below the active tackle, nBPA*F above it), the give
%   entering the closure as span growth (cAdd += scale_k * u_t'Cbowden*u_t
%   per unit F) and the geometry via an exit-row shift, exactly as (a).
%   Section axes are assumed aligned with the crossing's child-frame axes;
%   hand a full 3x3 via bowdenBracketCompliance for anything else.
%   The Shimano-type packaging envelope stays: BowdenBossDia = 7 mm boss,
%   BowdenRunClearance = 5 mm run clearance, both per-crossing TUNABLE
%   config fields (bowdenBossDia / bowdenRunClearance) -- directive
%   section 3, "revisit with real hardware dimensions".
%   PENDING-BEN: if a material is given WITHOUT a section, the stand-in
%   [b h L] = [25 5 40] mm plate is used and PRINTED once per
%   construction -- never silently (campaign rule 5).
%
% REGRESSION CONTRACT: with the default config (no mountCompliance /
%   mountMaterial / bowdenBracket* fields, i.e. Cmount = Cbowden =
%   zeros(3)) this class executes the IDENTICAL code path and arithmetic
%   of BiPam_pulley -- same psiFun values bit-for-bit, same fzero
%   brackets, same outputs. The established gate (2026-10-21 doctrine
%   round): load a recorded mat, re-evaluate its xBest through this
%   class, assert the recomputed torques/objective equal the recorded
%   finalObj values to 0.000. Fixture: Results/BiPulley_opensim_*_
%   20260926_*.mat (nBPA = 1, G = [1 1], moving_via batch).
%
% >>> END VARIANT-2 ADDITIONS; the rest of this header is the base
% class's own documentation, kept verbatim. <<<
%
% Description: Biarticular spin-off of MonoPam_pulley (which itself spins
% off MonoPamDataExplicit_balance): one BPA + stiffness pipeline routed
% across TWO joints, with the route/frame bookkeeping following BiPamData's
% indexing patterns exactly (Connor Morrow's biarticular class):
%   Location    Npts x 3, NO orientation dimension. Rows before Cross(1)
%               live in the proximal frame; rows Cross(1)..Cross(2)-1 in
%               the middle frame; rows Cross(2)..end in the distal frame.
%               When Cross(1) == Cross(2) (a via-less biarticular route,
%               e.g. gastrocnemius femur->calcn) the middle frame has no
%               rows and both crossings fall on the same row.
%   Cross       1x2 [c1 c2] row indices (first row of the next frame).
%   T           4x4xN1xN2: T(:,:,ii,1) maps middle-frame -> proximal-frame
%               coordinates (joint 1 at grid angle ii), T(:,:,iii,2) maps
%               distal -> middle (joint 2 at grid angle iii). Translations
%               must be the joint pivots (child-frame origin = joint), so
%               moment arms computed about the child-frame origin are about
%               the joint -- the same convention the mono classes inherit
%               from RpToTrans-built transforms.
% Per-crossing outputs follow BiPamData's layout: N1x3xN2x2, last axis =
% crossing index; route scalars are N1xN2 grids (ii rows, iii columns).
%
% DUPLICATE-CROSSING CORRECTION (documented deviation): BiPamData's
% UnitDirection branch for C(1) == C(2) transforms pointA with
% T(:,:,ii,2)*T(:,:,iii,2) (slot/index mix) and puts the crossing-1 moment
% arm about the MIDDLE-frame origin. This class instead expresses each
% crossing's endpoints in that crossing's CHILD frame via the consistent
% parent chain (frame-3 point -> frame 2 via T2; -> frame 1 via T1*T2), so
% the crossing-1 arm is about the frame-2 origin. Distinct crossings
% reduce to BiPamData's branches exactly.
%
% Transmission model (Ben's hardware contract; see MonoPam_pulley for the
% full derivation, implemented here once with the tackle serving the
% ACTIVE crossing):
%   Tendon span  per crossing k, the TENDON segment is the whole routed
%               path from the crossing-k exit row (PulleyExitIndex(k),
%               default Cross(k)-1) to the INSERTION (the last route row)
%               -- the ask's "exit-point to insertion" -- measured as the
%               sum of segment norms in crossing k's child frame.
%               DeltaL_k = span_k(theta0_k) - span_k, with theta0_k the
%               longest-span grid cell. When the rows before the exit all
%               live in the orientation-free proximal frame (always true
%               for crossing 1; for crossing 2 whenever the pre-exit rows
%               do not cross joint 1), the rigid contraction satisfies
%               s0ref_k + DeltaL_k = s0 exactly. On a general route with
%               the tackle at the distal crossing, the pre-exit rows'
%               BPA-side rigid change p = (s0 - s0ref) - DeltaL_k is
%               nonzero; the closure carries it EXPLICITLY by referencing
%               each cell to s0 - DeltaL_k instead of s0ref (identical
%               when p = 0), so at G = 1, nPulleyBPA = 1 the solved route
%               state still equals the base-class (reference-free) solve
%               at BOTH crossings, cell for cell, for ANY two-joint
%               motion. That is the regression identity
%               Opt_sanity_BiPulley.m asserts against MonoPam_pulley
%               analogs.
%   Closure     in MonoPam_pulley's tackle-travel form, per grid cell,
%                 G*(s - s0ref_k - cb*F(s)) = DeltaL + nBPA*F(s)/(G*kSpr),
%               with s0ref_k = s0 - DeltaL_k the per-cell rigid reference
%               (= s0ref whenever the pre-exit rows are orientation-free),
%               F(s) the festo4 force law of ONE BPA, cb the bracket
%               compliance projected on the BPA-side force direction,
%               nBPA = nPulleyBPA, and G = PulleyGain (the tackle
%               travel gain, = tackleLineParts unless overridden).
%   Solve       single scalar root find in s per grid cell (fzero,
%               bracketed to [0, KMAX*Rest]); psi(0) > 0 -> slack (tendon
%               transmits nothing, F_ins = 0), psi(KMAX*Rest) < 0 ->
%               PulleyInfeasible with NaN torque.
%   Active      ONE tackle acts. PulleyActive(k) = (gain(k) > 1); the
%               DISTAL crossing governs the equilibrium. Configuring BOTH
%               gains > 1 is rejected at construction (BiPam_pulley:
%               TwoTackles): two tackles are not chained and the crossing-1
%               mount reaction assumes it is rigid.
%   Segment     the tackle divides force only at/below itself: crossing
%   tensions    k >= active carries F_t = nBPA*F_single/G(active), the
%               proximal segments carry the full bundle pull nBPA*F_single
%               (per-BPA-side line). At G = 1 both coincide with the mono
%               bookkeeping. Torque_p per crossing uses its own tension.
%   Outputs     per crossing: F_ins = Ftendon(k).*u_t(k) (slack -> 0,
%               infeasible -> NaN), Torque_ins, and the mount reaction
%               ReactionF = nBPA(k)*F_single*u_bpa(k) + F_t(k)*u_t(k)
%               ONLY where a tackle physically sits (PulleyActive);
%               ReactionFmag is NaN where no tackle sits. PulleyTravel
%               = sCon - s0ref_k - F.*cb (s0ref_k = s0 - DeltaL, the
%               per-cell rigid reference; = s0ref when the pre-exit rows
%               are orientation-free) so that
%               G*PulleyTravel = DeltaL(active) + gama holds exactly
%               (MonoPam_pulley convention).
%   Routing     'moving_via' (default): u_t recomputed per grid cell from
%               the (deformed) exit point to the insertion -- the moving
%               exit point can INCREASE the moment arm as the joint
%               rotates (dissertation note). 'bowden': the housing anchor
%               sits on the distal body, so u_t is FIXED in the crossing's
%               child frame (held at its theta0 cell); the closure math is
%               unchanged. BowdenBossDia/BowdenRunClearance carry the
%               Shimano-type min-clearance envelope (M7 x 1.0 mm barrel
%               adjuster ~7 mm boss, 5 mm housing, 5 mm bracket plate);
%               Opt_run_BiPulley echoes it.
%
% The stiffness pipeline (Xi0/Xi1/Xi2, Wraps) acts on the WHOLE routed
% length as in the mono class: the bracket deforms the mount (row 1) under
% the BPA-side force, the route sum feeds the force law, and the per-
% crossing directions/arms are recomputed from the deformed geometry. The
% base class's hard-coded frame index 92 is made robust: row 1 lives in
% the orientation-free proximal frame, so pA = Location(1,:) directly.
%
% Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef BiPam_pulley2_xi3pose < handle

    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location                    %Npts x 3 (no orientation dim; BiPamData convention)
        Cross                       %1x2, first row of the middle and distal frames
        Diameter                    %Diameter of the BPA
        TransformationMat           %4x4xN1xN2 (two-hinge orientation grid)
        RestingL                    %Resting Length of the muscle
        Kmax                        %Length of BPA at maximum contraction
        FittingLength               %Length of each end cap (center of hole to bottom port)
        TendonL                     %Length of tendon, if any
        Pressure                    %Pressure of BPA

        % --- Stiffness parameters (always provided) ---
        Xi0                         %Constant length offset
        Xi1                         %Bracket axial stiffness
        Xi2                         %Bracket bending stiffness
        Wraps                       %Number of cable wraps (affects tendon stiffness)

        % --- xi3pose additions: per-pose wrapping loss (dissertation chi3)
        Xi3                         %Dimensionless wrapping-loss factor (arg 16,
                                    %  default 0 = base-identical)
        BendMeasure                 %Per-cell bend arc profile, m (arg 17; []
                                    %  = self-compute at the crossing contacts)
        WrapRadius                  %wRap bend-arc radius, m (config 'wrapRadius',
                                    %  default 0.025: buildGeoExclusion.m geo.wRap)

        % --- Reverse-pulley transmission config (per crossing) ---
        NPulleyBPA                  %1x2 BPAs in parallel feeding crossing k's tackle
        TackleLineParts             %1x2 rope parts supporting crossing k's moving block
        PulleyGain                  %1x2 travel gain G (default = TackleLineParts)
        RoutingMode                 %1x2 cell, 'moving_via' | 'bowden' per crossing
        PulleyExitIndex             %1x2 tendon exit row per crossing (default Cross-1)
        BowdenBossDia               %Bowden barrel-adjuster boss diameter, m (TUNABLE)
        BowdenRunClearance          %Bowden housing run clearance, m (TUNABLE)

        % --- VARIANT-2 (2026-10-22): moving-via mount + bowden bracket ---
        MountCompliance             %3x3 tackle-block mount compliance, m/N (active
                                    %  crossing's child frame); zeros(3) = off
        MountMaterial               %'none' | 'aluminum6061T6' | 'chromoly4130'
        MountSection                %[b h L] m cantilever section for MountMaterial
        BowdenBracketCompliance     %3x3 bowden anchor compliance, m/N (each bowden
                                    %  crossing's child frame); zeros(3) = off
        BowdenBracketMaterial       %'none' | 'aluminum6061T6' | 'chromoly4130'
        BowdenBracketSection        %[b h L] m cantilever section (PENDING-BEN stand-in
                                    %  default [0.025 0.005 0.040], printed when used)
        CmountEffective             %resolved 3x3 mount compliance actually used
        CbowdenEffective            %resolved 3x3 bowden bracket compliance used

        % --- xi3pose additions: two-bracket mechanics ---------------------
        % Origin bracket (row 1) stiffness assignment, N/m: the base class
        % hard-codes K = [X1 X2 X2] (App C "Biomimetic flexor" row); the
        % config field 'bracketStiffness1' may override it. Default keeps
        % the base arithmetic identical.
        BracketStiffness1           %3-vector [kx ky kz] N/m (default [Xi1 Xi2 Xi2])
        % Insertion bracket (route end, tendon anchor): dissertation ch.30
        % eq:bktstiffness assignment K = diag(chi1, chi2, chi1), two-rotation
        % bracket frame (eq:app_bracketframe). Config:
        % 'insertionBracketCompliance' (3x3 m/N, wins) |
        % 'insertionBracketStiffness' (3-vector N/m) |
        % 'insertionBracketRefPoint' ('auto' = per-cell exit row point so
        % x_br follows the load path, or fixed 1x3 m in the distal frame).
        % Default OFF (zeros(3)) = base-identical.
        InsBracketCompliance        %explicit 3x3 compliance, m/N ([] = unset)
        InsBracketStiffness         %3-vector N/m assignment ([] = unset)
        InsBracketRefPoint          %'auto' | 1x3 m distal-frame point
        CinsEffective               %resolved 3x3 insertion-bracket compliance

        % --- Stiffness-aware fields (minimizeFlx-style) ---
        L_p                         %Npts x 3 x N1 x N2 deformed locations (bracket moved row 1
                                    %  [and, variant 2, the deflecting exit rows])
        Funit                       %N x 3 force direction in the proximal frame (flattened grid)
        sL_p                        %N x Nseg deformed segment lengths (flattened grid)
        Lmt_p                       %N x 1 musculotendon length with deformed geometry and Xi0
        uD_p                        %N1 x 3 x N2 x 2 per-crossing force unit direction (deformed)
        strain_p                    %N x 1 contraction with bracket deformation, tendon stretch, Xi0
        FsingleBPA                  %N1 x N2 force of ONE BPA at the solved contraction, N
        segTension                  %N x 2 per-crossing segment tension (see class header)
        F_p                         %N1 x 3 x N2 x 2 per-crossing force vector (own tension)
        mA_p                        %N1 x 3 x N2 x 2 per-crossing moment arm (deformed)
        Torque_p                    %N1 x 3 x N2 x 2 per-crossing torque (NaN where infeasible/out-of-range strain)
        gama                        %N x 1 tendon stretch delta_t = F_t/kSpr
        kSpr                        %Tendon spring rate (effective, ONE tendon line)

        % --- Pulley-solved fields (grids N1 x N2 / N1 x 3 x N2) ---
        PulleyActive                %1x2 logical, a tackle sits at crossing k (gain > 1)
        ActiveCrossing              %Crossing whose transmission solved the route state
        deltaL                      %N1 x N2 x 2 rigid exit-to-insertion span change from theta0_k
        RefIndex                    %1x2 linear grid indices of theta0_k (longest crossing-k span)
        sContraction                %N1 x N2 solved absolute BPA contraction (closure variable s)
        PulleyTravel                %N1 x N2 tackle input travel sTravel = s - s0ref - cb*F
        Ftendon                     %N1 x N2 tackle tension F_t = nBPA*F_single/G (active gain)
        u_t                         %N1 x 3 x N2 x 2 tendon unit direction, insertion toward exit
        u_bpa                       %N1 x 3 x N2 x 2 BPA force unit direction in crossing k's child frame
        F_ins                       %N1 x 3 x N2 x 2 insertion force vector per crossing
        mA_ins                      %N1 x 3 x N2 x 2 moment arm of the tendon line per crossing
        Torque_ins                  %N1 x 3 x N2 x 2 torque from the insertion force
        ReactionF                   %N1 x 3 x N2 x 2 mount reaction (NaN where no tackle sits)
        ReactionFmag                %N1 x N2 magnitude of ReactionF (NaN where no tackle sits)
        PulleyInfeasible            %N1 x N2 logical, required contraction exceeded KMAX*Rest
        PulleySlack                 %N1 x N2 logical, tendon slack (solve clamped to s = 0)

        % --- VARIANT-2 deflection outputs ---
        MountDeflection             %N1 x 3 x N2 x 2 tackle-block exit-row shift from Cmount
                                    %  (active crossing's child frame; NaN where no tackle)
        BowdenDeflection            %N1 x 3 x N2 x 2 bowden anchor exit-row shift from Cbowden
                                    %  (each bowden crossing's child frame; NaN where not bowden)

        % --- VARIANT-2 closure diagnostics (zeros when deflections off) ---
        DeflClosure_cmB             %N x 1 mount BPA-side travel coefficient per unit F
        DeflClosure_cmT             %N x 1 mount tendon-span coefficient per unit F
        DeflClosure_cmBw            %N x 2 bowden anchor give coefficient per unit F
        DeflClosure_cAdd            %N x 1 total psiFun coefficient cAdd (m/N per unit F)

        % --- xi3pose additions: per-pose loss + bracket outputs -----------
        deltaRoute                  %N x 1 per-cell applied wrapping loss (m);
                                    %  eq:x3 per pose: Xi3*bend*comp^2
        BendArc                     %N x nContacts self-computed per-contact
                                    %  bend arcs (m; NaN when bendMeasure
                                    %  supplied; contacts = tendon-side
                                    %  interior vertices exit..npts-1)
        cbBracket                   %N x 1 origin-bracket force-path projection
                                    %  cb (u_hat' * C_bracket1 * u_hat, bracket frame)
        Bracket1Deflection          %N x 3 origin-bracket deflection (bracket frame)
        Bracket1Reaction            %N x 3 origin-bracket reaction, Fmag*u_hat
                                    %  (single-BPA convention, bracket frame)
        InsBracketDeflection        %N x 3 insertion-bracket deflection (distal
                                    %  frame; zeros when off)
        InsBracketReaction          %N x 3 insertion-bracket reaction Ft*u_t
                                    %  (distal frame; zeros when off / slack)
        DeflClosure_cIns            %N x 1 insertion-bracket give coefficient
                                    %  per unit F (span-growth side)
    end

    % Dependent properties are calculated from the explicit properties.
    properties (Dependent)
        SegmentLengths              %N1 x N2 x Nseg rigid (BiPamData layout)
        LongestSegment              %N1 x N2
        MuscleLength                %N1 x N2 rigid route length (BiPamData layout)
        Contraction                 %N1 x N2 rigid
        LengthCheck
        UnitDirection               %N1 x 3 x N2 x 2 rigid
        MomentArm                   %N1 x 3 x N2 x 2 rigid
        Fmax
        Force                       %N1 x 3 x N2 x 2 rigid
        Torque                      %N1 x 3 x N2 x 2 rigid
    end


    methods
        %% ------------- Muscle Data Constructor -----------------
        % Constructor function. The first fourteen arguments mirror
        % MonoPam_pulley; the fifteenth is the pulley config struct whose
        % fields may be SCALAR (same at both crossings) or 1x2 per
        % crossing: nPulleyBPA, tackleLineParts, gain, routingMode,
        % pulleyExitIndex, bowdenBossDia, bowdenRunClearance -- all
        % optional. Omitting it configures a straight-tendon pulley-
        % disabled transmission at both crossings.
        % VARIANT-2 optional fields (single transmission-level knobs, all
        % default OFF so the default construction is bit-identical to
        % BiPam_pulley): mountCompliance (3x3), mountMaterial +
        % mountSection ([b h L] m), bowdenBracketCompliance (3x3),
        % bowdenBracketMaterial + bowdenBracketSection ([b h L] m).
        % xi3pose OPTIONAL trailing arguments (DERIVED-CLASS DIFF 1): arg 16
        % xi3 (dimensionless wrapping-loss factor, default 0) and arg 17
        % bendMeasure (per-cell bend arc profile in m; [] = self-compute).
        % New pulleyConfig fields: wrapRadius (m), bracketStiffness1
        % (3-vector N/m), insertionBracketCompliance (3x3),
        % insertionBracketStiffness (3-vector N/m), insertionBracketRefPoint
        % ('auto' | 1x3 m).
        function PD = BiPam_pulley2_xi3pose(name, location, cross, diameter, t, ...
            rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, wraps, pulleyConfig, ...
            xi3, bendMeasure)
            if nargin < 16
                xi3 = 0;
            end
            if nargin < 17
                bendMeasure = [];
            end
            if nargin == 14
                pulleyConfig = struct();
            end
            if nargin >= 14 && nargin <= 17
                if ~isstruct(pulleyConfig)
                    error('BiPam_pulley2_xi3pose:PulleyConfig', ...
                        ['The 15th input must be the pulley config struct ' ...
                         '(nPulleyBPA, tackleLineParts, gain, routingMode, ' ...
                         'pulleyExitIndex, bowdenBossDia, bowdenRunClearance; ' ...
                         'scalar or 1x2 per crossing; variant-2: mountCompliance ' ...
                         '/ mountMaterial / mountSection / bowdenBracket*; ' ...
                         'xi3pose: wrapRadius / bracketStiffness1 / ' ...
                         'insertionBracket*).'])
                end
                PD.Name = name;                   % BPA/muscle name
                PD.Location = location;           % routing-point array (Npts x 3)
                PD.Cross = cross;                 % [c1 c2] crossing rows
                PD.Diameter = diameter;           % BPA diameter, mm
                PD.TransformationMat = t;         % frame transforms (4x4xN1xN2)
                PD.RestingL = rest;               % BPA resting length, m
                PD.Kmax = kmax;                   % fully contracted length, m
                PD.TendonL = tendon;              % tendon length, m
                PD.FittingLength = fitn;          % one fitting length, m
                PD.Pressure = pres;               % BPA pressure, kPa

                PD.Xi0 = xi0;                     % constant length offset, m
                PD.Xi1 = xi1;                     % axial bracket stiffness, N/m
                PD.Xi2 = xi2;                     % bending stiffness, N/m
                PD.Wraps = wraps;                 % tendon wrap count

                % -------- xi3pose config (derived-class diffs 1-3) ------
                PD.Xi3 = xi3;                     % dimensionless loss factor
                if ~(isscalar(PD.Xi3) && isfinite(PD.Xi3) && PD.Xi3 >= 0)
                    error('BiPam_pulley2_xi3pose:Xi3', ...
                        ['xi3 (arg 16) must be a finite scalar >= 0 ' ...
                         '(dimensionless wrapping-loss factor), got %g.'], ...
                        PD.Xi3)
                end
                PD.WrapRadius = cfgField2scalar(pulleyConfig, ...
                    'wrapRadius', 0.025);   % buildGeoExclusion.m geo.wRap
                if ~(isscalar(PD.WrapRadius) && isfinite(PD.WrapRadius) ...
                        && PD.WrapRadius >= 0)
                    error('BiPam_pulley2_xi3pose:WrapRadius', ...
                        ['wrapRadius must be a finite scalar >= 0 (m), ' ...
                         'got %g.'], PD.WrapRadius)
                end
                % BendMeasure shape/length validated in xi3RouteLoss (the
                % flattened cell count N lives with the grid there).
                if isempty(bendMeasure)
                    PD.BendMeasure = [];
                else
                    PD.BendMeasure = double(bendMeasure);
                end

                % Per-crossing config: scalar fields broadcast to 1x2.
                PD.NPulleyBPA = cfgField2(pulleyConfig, 'nPulleyBPA', [1, 1]);
                PD.TackleLineParts = cfgField2(pulleyConfig, ...
                    'tackleLineParts', [1, 1]);
                PD.PulleyGain = cfgField2(pulleyConfig, 'gain', ...
                    PD.TackleLineParts);
                PD.RoutingMode = cfgFieldCell2(pulleyConfig, 'routingMode', ...
                    {'moving_via', 'moving_via'});
                PD.PulleyExitIndex = cfgField2(pulleyConfig, ...
                    'pulleyExitIndex', PD.Cross - 1);
                PD.BowdenBossDia = cfgField2(pulleyConfig, ...
                    'bowdenBossDia', 0.007);   % M7 x 1.0 mm barrel adjuster
                PD.BowdenRunClearance = cfgField2(pulleyConfig, ...
                    'bowdenRunClearance', 0.005);  % 5 mm housing + plate

                % -------- VARIANT-2 config (directive sections 2-3) ----
                % All default OFF: zeros(3) compliance reproduces the base
                % class bit-for-bit (the regression gate of record).
                PD.MountCompliance = cfgFieldMtx(pulleyConfig, ...
                    'mountCompliance', []);
                PD.MountMaterial = cfgFieldStr(pulleyConfig, ...
                    'mountMaterial', 'none');
                [PD.MountSection, mountSectionDefaulted] = ...
                    cfgFieldVec3Flag(pulleyConfig, 'mountSection', ...
                    defaultBracketSection2());
                PD.BowdenBracketCompliance = cfgFieldMtx(pulleyConfig, ...
                    'bowdenBracketCompliance', []);
                PD.BowdenBracketMaterial = cfgFieldStr(pulleyConfig, ...
                    'bowdenBracketMaterial', 'none');
                [PD.BowdenBracketSection, bowdenSectionDefaulted] = ...
                    cfgFieldVec3Flag(pulleyConfig, 'bowdenBracketSection', ...
                    defaultBracketSection2());

                PD.MountMaterial = validateBracketMaterial2(PD.MountMaterial, ...
                    'mountMaterial');
                PD.BowdenBracketMaterial = validateBracketMaterial2( ...
                    PD.BowdenBracketMaterial, 'bowdenBracketMaterial');
                PD.MountSection = validateBracketSection2(PD.MountSection, ...
                    'mountSection', PD.MountMaterial);
                PD.BowdenBracketSection = validateBracketSection2( ...
                    PD.BowdenBracketSection, 'bowdenBracketSection', ...
                    PD.BowdenBracketMaterial);

                % Resolve the effective compliances (explicit matrix wins;
                % material + section build the cantilever model; neither
                % -> zeros(3) = OFF = base-class arithmetic exactly).
                PD.CmountEffective = resolveBracketCompliance2( ...
                    PD.MountCompliance, PD.MountMaterial, PD.MountSection, ...
                    mountSectionDefaulted, 'mount');
                PD.CbowdenEffective = resolveBracketCompliance2( ...
                    PD.BowdenBracketCompliance, PD.BowdenBracketMaterial, ...
                    PD.BowdenBracketSection, bowdenSectionDefaulted, ...
                    'bowden');
                % -------- xi3pose bracket config (derived-class diffs 4-5)
                % Origin bracket (row 1) assignment: default EXACTLY the
                % base class's hard-coded K = [X1 X2 X2] so the arithmetic
                % is unchanged unless overridden.
                if isfield(pulleyConfig, 'bracketStiffness1') ...
                        && ~isempty(pulleyConfig.bracketStiffness1)
                    PD.BracketStiffness1 = ...
                        double(pulleyConfig.bracketStiffness1);
                else
                    PD.BracketStiffness1 = [PD.Xi1, PD.Xi2, PD.Xi2];
                end
                if numel(PD.BracketStiffness1) ~= 3 ...
                        || ~all((isfinite(PD.BracketStiffness1) ...
                        | isinf(PD.BracketStiffness1)) ...
                        & PD.BracketStiffness1 > 0)
                    error('BiPam_pulley2_xi3pose:BracketStiffness1', ...
                        ['bracketStiffness1 must be a 3-vector of ' ...
                         'positive stiffnesses (N/m; Inf = rigid ' ...
                         'direction, the class''s own rigid convention), ' ...
                         'got [%s].'], num2str(PD.BracketStiffness1))
                end
                % Insertion bracket (route end): dissertation ch.30
                % eq:bktstiffness assignment diag(chi1, chi2, chi1) when the
                % caller passes insertionBracketStiffness; explicit 3x3
                % compliance wins; neither -> zeros(3) = OFF.
                PD.InsBracketCompliance = cfgFieldMtx(pulleyConfig, ...
                    'insertionBracketCompliance', []);
                PD.InsBracketStiffness = [];
                if isfield(pulleyConfig, 'insertionBracketStiffness') ...
                        && ~isempty(pulleyConfig.insertionBracketStiffness)
                    PD.InsBracketStiffness = ...
                        double(pulleyConfig.insertionBracketStiffness);
                end
                if isfield(pulleyConfig, 'insertionBracketRefPoint') ...
                        && ~isempty(pulleyConfig.insertionBracketRefPoint)
                    PD.InsBracketRefPoint = ...
                        pulleyConfig.insertionBracketRefPoint;
                else
                    PD.InsBracketRefPoint = 'auto';
                end
                if ~isempty(PD.InsBracketCompliance)
                    PD.CinsEffective = PD.InsBracketCompliance;
                elseif ~isempty(PD.InsBracketStiffness)
                    if numel(PD.InsBracketStiffness) ~= 3 ...
                            || ~all((isfinite(PD.InsBracketStiffness) ...
                            | isinf(PD.InsBracketStiffness)) ...
                            & PD.InsBracketStiffness > 0)
                        error('BiPam_pulley2_xi3pose:InsBracketStiffness', ...
                            ['insertionBracketStiffness must be a ' ...
                             '3-vector of positive stiffnesses (N/m; ' ...
                             'Inf = rigid direction), got [%s].'], ...
                            num2str(PD.InsBracketStiffness))
                    end
                    PD.CinsEffective = ...
                        diag(1 ./ PD.InsBracketStiffness);
                else
                    PD.CinsEffective = zeros(3);
                end
                if ischar(PD.InsBracketRefPoint)
                    if ~strcmpi(PD.InsBracketRefPoint, 'auto')
                        error('BiPam_pulley2_xi3pose:InsBracketRefPoint', ...
                            ['insertionBracketRefPoint must be ''auto'' ' ...
                             'or a fixed 1x3 point (m, distal frame), got ' ...
                             '''%s''.'], PD.InsBracketRefPoint)
                    end
                elseif numel(PD.InsBracketRefPoint) ~= 3 ...
                        || ~all(isfinite(PD.InsBracketRefPoint))
                    error('BiPam_pulley2_xi3pose:InsBracketRefPoint', ...
                        ['insertionBracketRefPoint must be ''auto'' or ' ...
                         'a finite 1x3 point (m, distal frame).'])
                end
                % -------- END xi3pose bracket config ---------------------

                for k = 1:2
                    if ~(isscalar(PD.NPulleyBPA(k)) && isfinite(PD.NPulleyBPA(k)) ...
                            && PD.NPulleyBPA(k) >= 1 && ...
                            PD.NPulleyBPA(k) == round(PD.NPulleyBPA(k)))
                        error('BiPam_pulley2_xi3pose:NPulleyBPA', ...
                            ['nPulleyBPA(%d) must be a positive integer ' ...
                             '(BPAs in parallel feeding the tackle), got %g.'], ...
                            k, PD.NPulleyBPA(k))
                    end
                    if ~(isscalar(PD.TackleLineParts(k)) && ...
                            isfinite(PD.TackleLineParts(k)) && ...
                            PD.TackleLineParts(k) >= 1)
                        error('BiPam_pulley2_xi3pose:TackleLineParts', ...
                            ['tackleLineParts(%d) must be >= 1 (rope parts ' ...
                             'supporting the moving block = travel gain G; ' ...
                             'a four-wheel tackle rigges to at most 4), got %g.'], ...
                            k, PD.TackleLineParts(k))
                    end
                    if ~(isscalar(PD.PulleyGain(k)) && isfinite(PD.PulleyGain(k)) ...
                            && PD.PulleyGain(k) >= 1)
                        error('BiPam_pulley2_xi3pose:GainRange', ...
                            ['PulleyGain(%d) must be >= 1 (1 = pulley ' ...
                             'disabled at that crossing), got %g.'], ...
                            k, PD.PulleyGain(k))
                    end
                    PD.RoutingMode{k} = validatestring(PD.RoutingMode{k}, ...
                        {'moving_via', 'bowden'}, mfilename, 'routingMode');
                    if ~(isscalar(PD.PulleyExitIndex(k)) && ...
                            PD.PulleyExitIndex(k) >= 1 && ...
                            PD.PulleyExitIndex(k) == round(PD.PulleyExitIndex(k)) && ...
                            PD.PulleyExitIndex(k) < PD.Cross(k))
                        error('BiPam_pulley2_xi3pose:ExitIndex', ...
                            ['PulleyExitIndex(%d) must name a row before ' ...
                             'Cross(%d) (1..%d), got %d.'], ...
                            k, k, PD.Cross(k) - 1, PD.PulleyExitIndex(k))
                    end
                end
                if PD.PulleyGain(1) > 1 && PD.PulleyGain(2) > 1
                    error('BiPam_pulley2_xi3pose:TwoTackles', ...
                        ['Both crossings carry gain > 1 (%g, %g), but only ' ...
                         'ONE tackle is supported: the crossing-1 mount ' ...
                         'reaction and the shared ReactionFmag grid assume ' ...
                         'the distal crossing is the only active one. Set ' ...
                         'the other crossing''s gain to 1.'], ...
                        PD.PulleyGain(1), PD.PulleyGain(2))
                end
                if numel(cross) ~= 2 || cross(1) < 2 || cross(2) < cross(1) ...
                        || cross(2) > size(location, 1)
                    error('BiPam_pulley2_xi3pose:CrossShape', ...
                        ['Cross must be [c1 c2] with 2 <= c1 <= c2 <= ' ...
                         'npts, got [%d %d] with %d points.'], ...
                        cross(1), cross(2), size(location, 1))
                end

                % Automatically compute stiffness-aware geometry and torque
                PD = PD.updateBiPulleyGeometry();
            else
                error('BiPam_pulley2_xi3pose:ArgCount', ...
                    ['Invalid number of arguments: expected 14 to 17 ' ...
                     '(...pulleyConfig, xi3, bendMeasure), got %d.'], nargin)
            end
        end

        %% ------------- Row-frame bookkeeping --------------------
        function f = rowFrame(obj, r)
        % Frame index (1, 2, 3) of route row r. Duplicate crossing
        % (Cross(1) == Cross(2)): the middle frame is empty and rows at or
        % beyond Cross(1) live in the distal frame.
            if r < obj.Cross(1)
                f = 1;
            elseif obj.Cross(1) == obj.Cross(2)
                f = 3;
            elseif r < obj.Cross(2)
                f = 2;
            else
                f = 3;
            end
        end

        %% ------------- Segment Lengths (rigid) ------------------
        % BiPamData's SegmentLengths generalized to the duplicate-crossing
        % case with the consistent parent chain (norms are frame-invariant,
        % this is bookkeeping).
        function segLengths = get.SegmentLengths(obj)
            [N1, N2] = obj.gridSize();
            nSeg = size(obj.Location, 1) - 1;
            segLengths = zeros(N1, N2, nSeg);

            for ii = 1:N1
                for iii = 1:N2
                    for i = 1:nSeg
                        pA = obj.pointInFrame(i, 1, ii, iii);
                        pB = obj.pointInFrame(i + 1, 1, ii, iii);
                        segLengths(ii, iii, i) = norm(pA - pB);
                    end
                end
            end
        end

        %% -------------- Longest Segment -------------------------
        function longestSeg = get.LongestSegment(obj)
            segLengths = obj.SegmentLengths;
            avgSegL = squeeze(mean(mean(segLengths, 1), 2));
            [~, longestSegPointer] = max(avgSegL);
            longestSeg = segLengths(:, :, longestSegPointer);
        end

        %% ------------- Muscle Length (rigid) --------------------
        function mL = get.MuscleLength(obj)
            mL = sum(obj.SegmentLengths, 3);
        end

        %% -------------- Force Unit Direction (rigid) ------------
        % BiPamData's UnitDirection generalized: crossing k's direction is
        % expressed in crossing k's CHILD frame via the consistent chain
        % (this is where the duplicate-crossing correction lives).
        function unitD = get.UnitDirection(obj)
            [N1, N2] = obj.gridSize();
            unitD = zeros(N1, 3, N2, 2);

            for k = 1:2
                child = k + 1;
                for ii = 1:N1
                    for iii = 1:N2
                        pA = obj.pointInFrame(obj.Cross(k) - 1, child, ii, iii);
                        pB = obj.pointInFrame(obj.Cross(k), child, ii, iii);
                        d = pA - pB;
                        unitD(ii, :, iii, k) = d / norm(d);
                    end
                end
            end
        end

        %% -------------- Moment Arm (rigid) ----------------------
        % Per crossing, about the child-frame origin (= the joint).
        function mA = get.MomentArm(obj)
            [N1, N2] = obj.gridSize();
            unitD = obj.UnitDirection;
            mA = zeros(N1, 3, N2, 2);

            for k = 1:2
                child = k + 1;
                for ii = 1:N1
                    for iii = 1:N2
                        pB = obj.pointInFrame(obj.Cross(k), child, ii, iii);
                        u = reshape(unitD(ii, :, iii, k), 1, 3);
                        mA(ii, :, iii, k) = pB - u * dot(u, pB);
                    end
                end
            end
        end

        %% -------------- Contraction (rigid) ---------------------
        function contraction = get.Contraction(obj)
            rest = obj.RestingL;
            tendon = obj.TendonL;
            fitting = obj.FittingLength;

            contraction = (rest - (obj.MuscleLength - tendon - 2 .* fitting)) ./ rest;
        end

        %% -------------- Length Check ----------------------------
        function lengthCheck = get.LengthCheck(obj)
            contraction = obj.Contraction;
            maxContractPercent = 0.25;          %Contracting to 75% of length
            minContractPercent = -0.03;         %Elongating to 103% of length
            restingPamLength = obj.RestingL;

            if restingPamLength < 0
                lengthCheck = 'Unusable';
            else
                if max(contraction(:)) <= maxContractPercent
                    if min(contraction(:)) >= minContractPercent
                        lengthCheck = 'Usable';
                    else
                        lengthCheck = 'Unusable';
                    end
                else
                    lengthCheck = 'Unusable';
                end
            end
        end

        %% -------------- Maximum Force ---------------------------
        function maxF = get.Fmax(obj)
            dia = obj.Diameter;
            rest = obj.RestingL;

            if dia == 10
                maxF = maxBPAforce(rest, '10');
            elseif dia == 20
                maxF = maxBPAforce(rest, '20');
            elseif dia == 40
                maxF = maxBPAforce(rest, '40');
            else
                disp('Wrong size diameter BPA')
            end
        end

        %% -------------- Force (rigid) ---------------------------
        function F = get.Force(obj)
            unitD = obj.UnitDirection;      % N1 x 3 x N2 x 2
            strain = obj.Contraction;       % N1 x N2

            rest = obj.RestingL;
            kmax = obj.Kmax;
            KMAX = (rest - kmax) / rest;

            rel = strain ./ KMAX;

            Fn = festo4(obj.Diameter, rel, obj.Pressure);
            scalarForce = Fn .* obj.Fmax;
            scalarForce(scalarForce < 0) = 0;
            scalarForce(scalarForce > obj.Fmax) = NaN;

            F = scalarForce .* unitD;
        end

        %% -------------- Torque (rigid) --------------------------
        function tor = get.Torque(obj)
            tor = cross(obj.MomentArm, obj.Force, 2);
        end

        %% ============================================================
        % === Stiffness-aware pipeline with the reverse pulley ========
        % ============================================================
        function obj = updateBiPulleyGeometry(obj)
            [N1, N2] = obj.gridSize();
            N = N1 * N2;                    % flattened grid cells

            % Tendon spring rate: ONE physical tendon line exits the
            % tackle, so kSpr is NOT multiplied by nPulleyBPA (mono
            % convention; the parallel BPAs share the single output line).
            obj.kSpr = Spr(obj, obj.Wraps);

            % Force unit vector in the proximal frame (origin to first
            % non-duplicate point, chained into frame 1), per grid cell.
            obj.Funit = computeForceVectorBi(obj);   % N x 3 flattened

            % Which crossing's transmission solves the route state:
            % distal-first (class-header assumption).
            obj.PulleyActive = obj.PulleyGain > 1;
            if obj.PulleyActive(2)
                activeK = 2;
            else
                activeK = 1;    % proximal tackle or fully rigid route
            end
            obj.ActiveCrossing = activeK;

            % Rigid exit-to-insertion tendon spans (per crossing) and
            % their longest-span references theta0_k. When the rows
            % before the exit are orientation-free (proximal frame),
            % s0ref_k + DeltaL_k = s0 holds exactly; on general distal-
            % tackle routes the BPA-side rigid change enters the closure
            % through the per-cell reference s0 - DeltaL_k (fortzBi /
            % class header).
            nPts = size(obj.Location, 1);
            obj.deltaL = zeros(N1, N2, 2);
            obj.RefIndex = zeros(1, 2);
            for k = 1:2
                child = k + 1;
                span = zeros(N1, N2);
                for ii = 1:N1
                    for iii = 1:N2
                        tot = 0;
                        for r = obj.PulleyExitIndex(k):nPts - 1
                            pA = obj.pointInFrame(r, child, ii, iii);
                            pB = obj.pointInFrame(r + 1, child, ii, iii);
                            tot = tot + norm(pA - pB);
                        end
                        span(ii, iii) = tot;
                    end
                end
                [~, linRef] = max(span(:));
                obj.RefIndex(k) = linRef;
                obj.deltaL(:, :, k) = span(linRef) - span;
            end

            % Rigid-geometry contraction per grid cell and at theta0 of
            % the active crossing (mono s0 formula).
            mL = obj.MuscleLength;                       % N1 x N2
            s0 = obj.RestingL - mL + obj.Xi0 + obj.TendonL + 2 * obj.FittingLength;
            s0ref = s0(obj.RefIndex(activeK));
            deltaLact = obj.deltaL(:, :, activeK);

            % Contraction from the constant length offset only
            % (base-identical pre-force for the bracket direction).
            strain_Xi0 = (obj.RestingL - (mL - obj.TendonL ...
                - 2 * obj.FittingLength - obj.Xi0)) ./ obj.RestingL;

            % Flatten the grid for the scalar-equilibrium helpers.
            strainXi0f = reshape(strain_Xi0, N, 1);
            s0f = reshape(s0, N, 1);
            deltaLf = reshape(deltaLact, N, 1);

            % -------- xi3pose: PER-CELL wrapping loss INSIDE the
            % equilibrium (derived-class diff 2) --------------------------
            % deltaRoute(i) = Xi3 * bend(i) * comp(i)^2, comp from the
            % Xi0-only pre-solve strain (validated MonoPam_mult.m
            % Contraction_k form, lines 445-449: strain0 excludes gama and
            % the loss; single pass - comp is mathematically independent
            % of the solve, so no fixed point is needed). Dissertation
            % eq:x3: l_m = l_LMT - 2*l_fitting - chi0 - Delta_l, so the
            % per-cell contraction demand grows by deltaRoute(i) BEFORE
            % LokBi2/fortzBi2 solves the closure. Xi0 keeps its own
            % meaning: it stays the constant offset above and is NOT
            % folded with the loss. With Xi3 = 0 the added term is exact
            % zero (IEEE no-op) and the path is bit-identical to base.
            obj.deltaRoute = xi3RouteLoss(obj, strainXi0f, N);
            s0f = s0f + obj.deltaRoute;

            % -------- VARIANT-2: deflection closure coefficients -------
            % Pre-solve per-cell directions (RIGID geometry, first-order,
            % exactly the role cb's pre-force direction plays) for the
            % moving-via tackle mount and any bowden anchor brackets:
            %   ubPre   BPA-side force direction from the BPA muscle
            %           ORIGIN (Location(1,:) -> Funit), mapped into the
            %           active crossing's child frame.
            %   utPre   tendon exit->insertion direction per crossing in
            %           that crossing's child frame (bowden crossings
            %           pinned at their theta0 cell, base convention).
            % From these, per-unit-F travel-absorption coefficients for
            % the equilibrium (see the variant-2 header block). With the
            % default zero compliances every coefficient is EXACTLY 0 and
            % the closure below is bit-for-bit the base class's.
            Gact = obj.PulleyGain(activeK);
            nBPAact = obj.NPulleyBPA(activeK);
            [ubPre, utPre] = preDeflectionDirections2(obj, activeK);
            [cAdd, cmBf, cmTf, cmBwf, cInsf] = deflectionClosureCoeffs2( ...
                obj, ubPre, utPre, activeK);

            % Deformed geometry, tendon stretch, and the pulley
            % equilibrium (single route state; see class header). The
            % variant passes the per-cell deflection coefficient cAdd
            % into the solve (fortzBi2: psi += cAdd*F(s)).
            [L_p_f, gama_f, sCon_f, Fmag_f, infeas_f, slack_f, cb_f, ...
                dBr1_f, uHat1_f] = LokBi2( ...
                obj, obj.Xi1, obj.Xi2, obj.kSpr, obj.Funit, ...
                strainXi0f, obj.Xi0, deltaLf, s0f, s0ref, cAdd);

            % Deformed locations back onto the 4-D grid.
            obj.L_p = reshape(L_p_f, [size(obj.Location, 1), 3, N1, N2]);
            obj.gama = gama_f;
            obj.sContraction = reshape(sCon_f, N1, N2);
            obj.FsingleBPA = reshape(Fmag_f, N1, N2);

            % Tackle tension and travel at the active crossing
            % (MonoPam_pulley conventions).
            obj.Ftendon = obj.FsingleBPA * nBPAact / Gact;
            % Per-cell rigid reference s0 - DeltaL (equals s0ref exactly
            % when the pre-exit rows are orientation-free; class header).
            s0refGrid = reshape(s0f - deltaLf, N1, N2);
            obj.PulleyTravel = reshape(sCon_f, N1, N2) - s0refGrid ...
                - obj.FsingleBPA .* reshape(cb_f, N1, N2);

            obj.PulleyInfeasible = reshape(infeas_f, N1, N2);
            obj.PulleySlack = reshape(slack_f, N1, N2);
            % Slack tendon carries no tension: the solve keeps the
            % zero-strain BPA force, so without this the taut tension
            % would leak into F_ins and the mount reaction at slack
            % cells (MonoPam_pulley convention).
            obj.Ftendon(obj.PulleySlack) = 0;

            % -------- VARIANT-2: apply the deflections geometrically ---
            % delta = C * R at the solved force state; the exit rows of
            % L_p are shifted (rotation-only mapping into each row's own
            % frame, so the displacement picks up no joint translation),
            % so the downstream u_t / mA_ins / deformed segment lengths
            % honor the deflected moving via / bowden anchor. Skipped
            % entirely (zero shift) when the compliances are zero, which
            % keeps the default path bit-identical.
            obj.MountDeflection = nan(N1, 3, N2, 2);
            obj.BowdenDeflection = nan(N1, 3, N2, 2);
            obj = applyExitDeflections2(obj, ubPre, utPre, activeK);
            % -------- END VARIANT-2 deflection application -------------

            % -------- xi3pose: bracket outputs (derived-class diff 6) ---
            % Origin bracket (row 1) per-cell deflection and reaction in
            % the bracket frame (returned by LokBi2; zeros in the rigid
            % escape). Reaction = Fmag*u_hat, the single-BPA convention
            % the base closure's cb/e_bkt terms use.
            obj.cbBracket = cb_f;
            obj.Bracket1Deflection = dBr1_f;
            obj.Bracket1Reaction = dBr1_f;      % overwritten right below
            obj.Bracket1Reaction = Fmag_f .* uHat1_f;
            % Insertion bracket (route end): geometric give under the
            % tendon-line tension (separate hardware from the tackle mount
            % and the Bowden housing stops; counted ONCE in the closure
            % via DeflClosure_cIns). zeros when off.
            obj = applyInsBracketDeflection2(obj, utPre, activeK);
            % -------- END xi3pose bracket outputs -----------------------

            % Deformed segment lengths and musculotendon length.
            obj.sL_p = segBi(obj);
            Lmt_p_f = sum(obj.sL_p, 2) - obj.Xi0;
            obj.Lmt_p = Lmt_p_f;

            % Route contraction with bracket deformation, tendon stretch,
            % Xi0, and the per-cell Xi3 wrapping loss (the fictive extra
            % contraction the force law saw; MonoPam_mult strain_f analog,
            % Contraction_k with the SOLVED gama). deltaRoute = 0 keeps
            % this bit-identical to base.
            obj.strain_p = (obj.RestingL - (Lmt_p_f - obj.TendonL - gama_f ...
                - 2 * obj.FittingLength - obj.deltaRoute)) ./ obj.RestingL;

            % Per-crossing directions, arms, and segment tensions.
            obj.uD_p = zeros(N1, 3, N2, 2);
            obj.mA_p = zeros(N1, 3, N2, 2);
            obj.segTension = zeros(N, 2);
            for k = 1:2
                child = k + 1;
                for ii = 1:N1
                    for iii = 1:N2
                        pA = obj.pointInFrameDeformed(obj.Cross(k) - 1, child, ii, iii);
                        pB = obj.pointInFrameDeformed(obj.Cross(k), child, ii, iii);
                        d = pA - pB;
                        u = d / norm(d);
                        obj.uD_p(ii, :, iii, k) = u;
                        obj.mA_p(ii, :, iii, k) = pB - u * dot(u, pB);
                    end
                end
                % Tension: at/below the active tackle the single tendon
                % line carries F_t = nBPA*F_single/G; the proximal
                % segments (BPA side of the tackle) carry the full bundle
                % pull nBPA*F_single.
                if k >= activeK
                    obj.segTension(:, k) = Fmag_f * nBPAact / Gact;
                else
                    obj.segTension(:, k) = Fmag_f * nBPAact;
                end
            end

            % Per-crossing force vectors and torques. Torque NaN mirrors
            % the base Tor rule (strain below the floor) plus the pulley
            % infeasibility flag.
            strainGrid = reshape(obj.strain_p, N1, N2);
            obj.F_p = zeros(N1, 3, N2, 2);
            for k = 1:2
                obj.F_p(:, :, :, k) = ...
                    reshape(obj.segTension(:, k), N1, 1, N2) .* obj.uD_p(:, :, :, k);
            end
            obj.Torque_p = cross(obj.mA_p, obj.F_p, 2);
            nanMask2 = obj.PulleyInfeasible | (strainGrid < obj.strainFloor());
            for k = 1:2
                obj.Torque_p(:, :, :, k) = maskNaN3(obj.Torque_p(:, :, :, k), nanMask2);
            end

            % --------------------------------------------------------
            % Reverse-pulley outputs (equilibrium-honest quantities)
            % --------------------------------------------------------
            obj.u_t = zeros(N1, 3, N2, 2);
            obj.u_bpa = zeros(N1, 3, N2, 2);
            obj.F_ins = zeros(N1, 3, N2, 2);
            obj.mA_ins = zeros(N1, 3, N2, 2);
            obj.Torque_ins = zeros(N1, 3, N2, 2);
            obj.ReactionF = nan(N1, 3, N2, 2);   % NaN where no tackle sits
            obj.ReactionFmag = nan(N1, N2);

            for k = 1:2
                child = k + 1;

                % Tendon-side unit direction from the (deformed) exit
                % point to the INSERTION (route end), in crossing k's
                % child frame, recomputed per grid cell. The force on the
                % distal body pulls the insertion toward the exit (base
                % UD sign convention: proximal-side point minus distal
                % point).
                u_t_i = zeros(N1, 3, N2);
                for ii = 1:N1
                    for iii = 1:N2
                        pE = obj.pointInFrameDeformed( ...
                            obj.PulleyExitIndex(k), child, ii, iii);
                        pI = obj.pointInFrameDeformed(nPts, child, ii, iii);
                        u_t_i(ii, :, iii) = pE - pI;
                    end
                end
                u_t_i = normalize3(u_t_i);
                if strcmp(obj.RoutingMode{k}, 'bowden')
                    % Housing anchor on the distal body: the line
                    % direction is fixed in the child frame (theta0).
                    % RefIndex is a LINEAR index over the (ii, iii) grid
                    % with ii fastest; the three components of that cell
                    % inside the N1 x 3 x N2 array live at
                    % ii0 + (c-1)*N1 + (iii0-1)*3*N1 (c = 1..3) -- use
                    % sub2ind so no stride is missed (consecutive linear
                    % indices, or RefIndex + (c-1)*N1 alone, mix
                    % DIFFERENT cells' components).
                    ii0 = mod(obj.RefIndex(k) - 1, N1) + 1;
                    iii0 = floor((obj.RefIndex(k) - 1) / N1) + 1;
                    cellIdx = sub2ind([N1, 3, N2], ii0, 1:3, iii0);
                    vref = reshape(u_t_i(cellIdx), 1, 3);
                    u_t_i = repmat(reshape(vref, 1, 3, 1), [N1, 1, N2]);
                end
                obj.u_t(:, :, :, k) = u_t_i;

                % Insertion force = the crossing's own segment tension
                % times u_t; slack transmits nothing; infeasible is NaN
                % (MonoPam_pulley conventions).
                F_ins_i = reshape(obj.segTension(:, k), N1, 1, N2) .* u_t_i;
                F_ins_i = maskZero3(F_ins_i, obj.PulleySlack);
                F_ins_i = maskNaN3(F_ins_i, obj.PulleyInfeasible);
                obj.F_ins(:, :, :, k) = F_ins_i;

                % Moment arm of the tendon line about the joint (the
                % perpendicular foot is line-invariant; the insertion
                % point plays the base class's crossing-point role).
                mA_ins_i = zeros(N1, 3, N2);
                for ii = 1:N1
                    for iii = 1:N2
                        pI = obj.pointInFrameDeformed(nPts, child, ii, iii);
                        u = reshape(u_t_i(ii, :, iii), 1, 3);
                        mA_ins_i(ii, :, iii) = pI - u * dot(u, pI);
                    end
                end
                obj.mA_ins(:, :, :, k) = mA_ins_i;

                Torque_ins_i = cross(mA_ins_i, F_ins_i, 2);
                Torque_ins_i = maskNaN3(Torque_ins_i, obj.PulleyInfeasible);
                obj.Torque_ins(:, :, :, k) = Torque_ins_i;

                if ~obj.PulleyActive(k)
                    continue;   % no tackle at this crossing: no reaction
                end

                % BPA-side force direction in crossing k's child frame so
                % the reaction sum is frame-consistent. DIRECTIONS rotate
                % between frames but must NOT pick up the joint
                % translations, so this uses the rotation-only twin
                % mapDirToFrame (a point transform would corrupt the
                % direction by the joint offset).
                u_bpa_i = zeros(N1, 3, N2);
                for ii = 1:N1
                    for iii = 1:N2
                        u1 = reshape(obj.Funit(ii + (iii - 1) * N1, :), 1, 3);
                        u_bpa_i(ii, :, iii) = obj.mapDirToFrame( ...
                            u1, child, ii, iii, 1);
                    end
                end
                u_bpa_i = normalize3(u_bpa_i);
                obj.u_bpa(:, :, :, k) = u_bpa_i;

                % Mount reaction at the pulley block: the nPulleyBPA
                % parallel BPA lines plus the tendon line (NaN'd where
                % infeasible; MonoPam_pulley convention). The grid
                % scalars are reshaped to N1 x 1 x N2 for implicit
                % expansion against the N1 x 3 x N2 direction arrays.
                Fbpa3 = reshape(obj.FsingleBPA, N1, 1, N2);
                Ft3 = reshape(obj.Ftendon, N1, 1, N2);
                R = (nBPAact * Fbpa3) .* u_bpa_i + Ft3 .* u_t_i;
                R = maskNaN3(R, obj.PulleyInfeasible);
                obj.ReactionF(:, :, :, k) = R;
                obj.ReactionFmag = reshape(vecnorm(R, 2, 2), N1, N2);
            end

            % -------- VARIANT-2 diagnostics stash (closure coefficients)
            % Kept on the object for audits/tests; zeros when off.
            obj.DeflClosure_cmB = cmBf;
            obj.DeflClosure_cmT = cmTf;
            obj.DeflClosure_cmBw = cmBwf;
            obj.DeflClosure_cIns = cInsf;
            obj.DeflClosure_cAdd = cAdd;
            % -------- END VARIANT-2 diagnostics ------------------------
        end

        %% -------------- Grid size -------------------------------
        function [N1, N2] = gridSize(obj)
            N1 = size(obj.TransformationMat, 3);
            N2 = size(obj.TransformationMat, 4);
        end

        %% -------------- Strain floor (base Tor rule) ------------
        function ss = strainFloor(obj)
            switch obj.Diameter
                case 20
                    ss = -.03;       %maximum allowable strain
                case 10
                    ss = -.02;
                otherwise
                    ss = -.02;
            end
        end

        %% ------ Point/coordinate frame mapping (rigid) ----------
        function p = pointInFrame(obj, r, frame, ii, iii)
        % Route row r expressed in the requested frame at grid (ii, iii).
            p = obj.mapVectorToFrame(obj.Location(r, :), frame, ii, iii, ...
                obj.rowFrame(r));
        end

        function p = pointInFrameDeformed(obj, r, frame, ii, iii)
        % Same, but reading the DEFORMED location grid L_p.
            p = obj.mapVectorToFrame( ...
                reshape(obj.L_p(r, :, ii, iii), 1, 3), frame, ii, iii, ...
                obj.rowFrame(r));
        end

        function v = mapVectorToFrame(obj, v, frame, ii, iii, sourceFrame)
        % Express a vector/point given in sourceFrame in the requested
        % frame, chaining T1/T2 consistently (parent chain downward, inverse
        % chain upward). POINT transform: translations included.
            T1 = obj.TransformationMat(:, :, ii, 1);   % frame2 -> frame1
            T2 = obj.TransformationMat(:, :, iii, 2);  % frame3 -> frame2
            f = sourceFrame;
            while f > frame
                if f == 3
                    v = RowVecTrans(T2, v);
                    f = 2;
                else
                    v = RowVecTrans(T1, v);
                    f = 1;
                end
            end
            while f < frame
                if f == 1
                    v = RowVecTrans(inv(T1), v);
                    f = 2;
                else
                    v = RowVecTrans(inv(T2), v);
                    f = 3;
                end
            end
        end

        function v = mapDirToFrame(obj, v, frame, ii, iii, sourceFrame)
        % Rotation-only twin of mapVectorToFrame for DIRECTIONS: a
        % direction rotates between frames but must not pick up the
        % joint translations (mapping a direction with the point
        % transform corrupts it by the joint offset).
            T1 = obj.TransformationMat(:, :, ii, 1);   % frame2 -> frame1
            T2 = obj.TransformationMat(:, :, iii, 2);  % frame3 -> frame2
            f = sourceFrame;
            while f > frame
                if f == 3
                    v = RowVecDir(T2, v);
                    f = 2;
                else
                    v = RowVecDir(T1, v);
                    f = 1;
                end
            end
            while f < frame
                if f == 1
                    v = RowVecDir(inv(T1), v);
                    f = 2;
                else
                    v = RowVecDir(inv(T2), v);
                    f = 3;
                end
            end
        end

    end % methods

end % classdef

%% =====================================================================
%% Helper functions (derived from the minimizer calculations; flattened
%% N = N1*N2 grid-cell formulation)
%% =====================================================================

%% ------------- Per-crossing config field (scalar or 1x2) ---------
function val = cfgField2(cfg, name, default)
if isfield(cfg, name)
    val = cfg.(name);
    if isscalar(val)
        val = [val, val];       %scalar config: same at both crossings
    end
    val = reshape(val, 1, 2);
else
    val = default;
end
end

%% ------------- Per-crossing config field (cell, scalar or 1x2) ---
function val = cfgFieldCell2(cfg, name, default)
if isfield(cfg, name)
    v = cfg.(name);
    if ischar(v)
        val = {char(v), char(v)};
    elseif isscalar(v)
        val = {char(v), char(v)};
    else
        val = {char(v{1}), char(v{2})};
    end
else
    val = default;
end
end

%% ------------- VARIANT-2 config helpers ---------------------------
function val = cfgFieldStr(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = char(cfg.(name));
else
    val = default;
end
end

function [val, defaulted] = cfgFieldVec3Flag(cfg, name, default)
defaulted = false;
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = double(cfg.(name));
    val = reshape(val, 1, 3);
else
    val = default;
    defaulted = true;
end
end

function val = cfgFieldMtx(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = double(cfg.(name));
    if ~isequal(size(val), [3, 3])
        error('BiPam_pulley2_xi3pose:%sShape', ...
            '%s must be a 3x3 compliance matrix (m/N), got %dx%d.', ...
            name, size(val, 1), size(val, 2))
    end
    if ~all(isfinite(val(:)))
        error('BiPam_pulley2_xi3pose:%sFinite', ...
            '%s must be all-finite.', name)
    end
else
    val = default;
end
end

function m = validateBracketMaterial2(m, fieldName)
m = validatestring(lower(m), ...
    {'none', 'aluminum6061t6', 'chromoly4130'}, mfilename, fieldName);
end

function s = validateBracketSection2(s, fieldName, material)
% The section only matters when a material is active; validate shape
% whenever a real (non-default-empty) value is present.
if ~strcmp(material, 'none') || ~isequal(s, [])
    if numel(s) ~= 3 || ~all(isfinite(s)) || any(s <= 0)
        error('BiPam_pulley2_xi3pose:%sShape', ...
            ['%s must be [b h L] in meters, all finite and positive ' ...
             '(cantilever bracket section), got [%s].'], ...
            fieldName, num2str(s))
    end
end
end

function C = resolveBracketCompliance2(Cexplicit, material, section, ...
    sectionDefaulted, which)
% Explicit 3x3 wins; otherwise material + cantilever section; otherwise
% zeros(3) = OFF. The PENDING-BEN default section is PRINTED whenever it
% is actually used (campaign rule 5: never guessed silently).
if ~isempty(Cexplicit)
    C = Cexplicit;
elseif ~strcmp(material, 'none')
    if sectionDefaulted
        fprintf(['BiPam_pulley2_xi3pose: %s bracket section ' ...
            'PENDING-BEN stand-in [b h L] = [%g %g %g] mm (%s, E = ' ...
            '%.1f GPa) -- set %sSection for real hardware dims.\n'], ...
            which, 1000 * section(1), 1000 * section(2), ...
            1000 * section(3), material, materialE2(material) / 1e9, which)
    end
    C = bracketCompliance2(material, section);
else
    C = zeros(3);
end
end

function E = materialE2(material)
% Bicycle brake/shifter bracket materials (directive section 3):
% 6061-T6 aluminum and 4130 CroMoly (chromoly) steel.
switch lower(material)
    case 'aluminum6061t6'
        E = 68.9e9;     %Pa
    case 'chromoly4130'
        E = 205e9;      %Pa
    otherwise
        error('BiPam_pulley2_xi3pose:MaterialE', ...
            'No Young''s modulus for material ''%s''.', material)
end
end

function C = bracketCompliance2(material, section)
% Cantilever bracket compliance from material + rectangular section
% [b h L] (m), section axes aligned with the compliance frame axes:
%   axial   k_a  = E*A/L,             A  = b*h
%   bending k_by = 3*E*Iy/L^3,        Iy = b*h^3/12
%   bending k_bz = 3*E*Iz/L^3,        Iz = h*b^3/12
% C = diag(1/k_a, 1/k_by, 1/k_bz). For anything else (angled mounts,
% coupled sections) hand the full 3x3 in via the *Compliance config.
E = materialE2(material);
b = section(1);
h = section(2);
L = section(3);
A = b * h;
Iy = b * h^3 / 12;
Iz = h * b^3 / 12;
ka = E * A / L;
kby = 3 * E * Iy / L^3;
kbz = 3 * E * Iz / L^3;
C = diag([1 / ka, 1 / kby, 1 / kbz]);
end

function s = defaultBracketSection2()
% PENDING-BEN stand-in plate section [b h L] = [25 5 40] mm (directive
% section 3: "revisit with real hardware dimensions"; printed when used).
s = [0.025, 0.005, 0.040];
end

%% ------------- VARIANT-2 deflection helpers -----------------------
function [ubPre, utPre] = preDeflectionDirections2(obj, activeK)
% Pre-solve per-cell directions on the RIGID geometry (first-order, the
% same role cb's pre-force direction plays in the base closure):
%   ubPre  N1x3xN2    BPA-side force direction, from the BPA muscle
%                    ORIGIN (Location(1,:) -- the same origin the row-1
%                    Xi1/Xi2 bracket references; Funit is built from it
%                    in computeForceVectorBi), mapped rotation-only into
%                    the ACTIVE crossing's child frame.
%   utPre  N1x3xN2x2  per-crossing rigid exit->insertion tendon
%                    direction in that crossing's child frame; bowden
%                    crossings pinned at their theta0 cell (base u_t
%                    convention, incl. the sub2ind stride fix).
[N1, N2] = obj.gridSize();
nPts = size(obj.Location, 1);

ubPre = zeros(N1, 3, N2);
for ii = 1:N1
    for iii = 1:N2
        u1 = reshape(obj.Funit(ii + (iii - 1) * N1, :), 1, 3);
        ubPre(ii, :, iii) = obj.mapDirToFrame(u1, activeK + 1, ii, iii, 1);
    end
end
ubPre = normalize3(ubPre);

utPre = zeros(N1, 3, N2, 2);
for k = 1:2
    child = k + 1;
    utk = zeros(N1, 3, N2);
    for ii = 1:N1
        for iii = 1:N2
            pE = obj.pointInFrame(obj.PulleyExitIndex(k), child, ii, iii);
            pI = obj.pointInFrame(nPts, child, ii, iii);
            utk(ii, :, iii) = pE - pI;
        end
    end
    utk = normalize3(utk);
    if strcmp(obj.RoutingMode{k}, 'bowden')
        ii0 = mod(obj.RefIndex(k) - 1, N1) + 1;
        iii0 = floor((obj.RefIndex(k) - 1) / N1) + 1;
        cellIdx = sub2ind([N1, 3, N2], ii0, 1:3, iii0);
        vref = reshape(utk(cellIdx), 1, 3);
        utk = repmat(reshape(vref, 1, 3, 1), [N1, 1, N2]);
    end
    utPre(:, :, :, k) = utk;
end
end

function [cAdd, cmBf, cmTf, cmBwf, cInsf] = deflectionClosureCoeffs2(obj, ...
    ubPre, utPre, activeK)
% Per-cell travel-absorption coefficients per unit single-BPA force F
% (see the variant-2 header block):
%   tackle mount (active crossing only, needs a tackle):
%     w   = u_bpa + u_t/G           reaction direction per unit F
%     cmB = u_bpa' * Cmount * w     BPA-side travel absorbed per unit F
%     cmT = u_t'   * Cmount * w     tendon span growth per unit F
%   bowden anchors (every bowden crossing):
%     cmBw(k) = u_t' * Cbowden * u_t   anchor give per unit tension
%   INSERTION bracket (xi3pose, derived-class diff 4): separate hardware
%     from the mount/bowden stops; give at the route end under the
%     tendon-line tension (crossing 2's segment tension, scale nBPA/G
%     always since activeK <= 2):
%       cIns = (nBPA/G) * u_t' * Cins * u_t
%     projected through the per-cell two-rotation bracket frame port
%     (yaw-pitch x_br along the load path; eq:app_bracketframe).
% with tension scale nBPA/G at/below the active tackle, nBPA above it:
%   cAdd = -G*nBPA*cmB + nBPA*cmT + sum_k scale_k * cmBw(k) + cIns
% psiFun2 = psiFun_base + cAdd*F(s); every coefficient is EXACTLY zero
% when the compliances are zero, so the default is bit-identical.
[N1, N2] = obj.gridSize();
N = N1 * N2;
Gact = obj.PulleyGain(activeK);
nBPAact = obj.NPulleyBPA(activeK);
Cm = obj.CmountEffective;
Cb = obj.CbowdenEffective;

cmBf = zeros(N, 1);
cmTf = zeros(N, 1);
cmBwf = zeros(N, 2);

if any(obj.PulleyActive) && any(Cm(:))
    utA = utPre(:, :, :, activeK);
    for ii = 1:N1
        for iii = 1:N2
            lin = ii + (iii - 1) * N1;
            uB = reshape(ubPre(ii, :, iii), 3, 1);
            uT = reshape(utA(ii, :, iii), 3, 1);
            w = uB + uT / Gact;
            cmBf(lin) = uB.' * (Cm * w);
            cmTf(lin) = uT.' * (Cm * w);
        end
    end
end

for k = 1:2
    if ~strcmp(obj.RoutingMode{k}, 'bowden') || ~any(Cb(:))
        continue;
    end
    if k >= activeK
        scale_k = nBPAact / Gact;
    else
        scale_k = nBPAact;
    end
    utk = utPre(:, :, :, k);
    for ii = 1:N1
        for iii = 1:N2
            lin = ii + (iii - 1) * N1;
            uT = reshape(utk(ii, :, iii), 3, 1);
            cmBwf(lin, k) = scale_k * (uT.' * (Cb * uT));
        end
    end
end

cInsf = insBracketClosureCoeffs2(obj, utPre);
cAdd = -Gact * nBPAact * cmBf + nBPAact * cmTf + sum(cmBwf, 2) + cInsf;
end

%% ------------- xi3pose: insertion-bracket closure coefficient -----
function cInsf = insBracketClosureCoeffs2(obj, utPre)
% Per-cell give coefficient of the INSERTION bracket (route end) under
% the tendon-line tension, projected through a per-cell port of the
% dissertation's two-rotation bracket frame (eq:app_bracketframe: yaw
% about the frame z, pitch about y_br, x_br aligned with the vector from
% the bracket reference point to the attachment point; "the frame's x
% axis follows the nominal load path"). The ACTUAL load direction (the
% pre-solve rigid exit->insertion unit vector utPre crossing 2, distal
% frame) is transformed into the per-cell bracket frame and projected on
% diag compliance -- the same force-into-frame projection the origin
% bracket's LokBi2 does. With the default 'auto' reference the x axis
% aligns with the load path, so the projection reduces to the axial 1/kx
% of the dissertation assignment diag(chi1, chi2, chi1). Tension scale =
% the tendon line at the insertion = nBPA/G (crossing 2 >= activeK
% always; the caller's PulleyGain/NPulleyBPA at the ACTIVE crossing is
% the scale used throughout the variant-2 coefficients).
[N1, N2] = obj.gridSize();
N = N1 * N2;
Gact = obj.PulleyGain(obj.ActiveCrossing);
nBPAact = obj.NPulleyBPA(obj.ActiveCrossing);
Cins = obj.CinsEffective;
cInsf = zeros(N, 1);
if ~any(Cins(:))
    return
end
nPts = size(obj.Location, 1);
exR = obj.PulleyExitIndex(2);
useRef = ~ischar(obj.InsBracketRefPoint);
utk = utPre(:, :, :, 2);
for ii = 1:N1
    for iii = 1:N2
        lin = ii + (iii - 1) * N1;
        pI = obj.pointInFrame(nPts, 3, ii, iii);
        if useRef
            pRef = obj.InsBracketRefPoint;
        else
            pRef = obj.pointInFrame(exR, 3, ii, iii);
        end
        Rhbr = insBracketFrameRot2(pI, pRef);
        uT = reshape(utk(ii, :, iii), 3, 1);
        uHat = (Rhbr' * uT).';   %load direction in bracket-frame coords
        cInsf(lin) = (nBPAact / Gact) * (uHat * (Cins * uHat.'));
    end
end
end

%% ------------- xi3pose: insertion-bracket frame rotation ----------
function Rhbr = insBracketFrameRot2(pI, pRef)
% Yaw-pitch bracket frame rotation (port of the origin bracket's
% construction in LokBi2, per cell): yaw about the fixed z axis aligns
% the sagittal projection of the reference->attachment vector; pitch
% about the current y_br aligns the full vector with x_br. Rhbr maps
% bracket coords to distal coords.
phbr = pI(:) - pRef(:);
thetabrA = atan2(phbr(2), phbr(1));
RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...
    sin(thetabrA) cos(thetabrA) 0; ...
    0    0   1];
pbrhA = RhbrZ' * phbr;
thetaY = atan2(pbrhA(3), pbrhA(1));
Ry = [cos(thetaY)  0  sin(thetaY);
    0            1  0;
    -sin(thetaY) 0  cos(thetaY)];
Rhbr = RhbrZ * Ry';
end

%% ------------- xi3pose: insertion-bracket deflection --------------
function obj = applyInsBracketDeflection2(obj, utPre, ~)
% Post-solve geometric application (derived-class diff 4): the route-end
% row shifts by Cins * (Ftendon * u_t) under the tendon-line tension
% (slack already zeroed in Ftendon: slack cells give nothing; infeasible
% cells keep their deflection, their torque is NaN'd downstream). The
% last row lives in the distal frame = crossing 2's child frame, so the
% rotation-only mapping into the row's own frame is the identity. The
% direction is utPre(:, :, :, 2) -- the pre-solve rigid exit->insertion
% direction (first-order, exactly the role the variant-2 mount/bowden
% terms play; obj.u_t does not exist yet at this pipeline point).
% Skipped entirely (zero shift) when the compliance is zeros(3), keeping
% the default path bit-identical.
Cins = obj.CinsEffective;
obj.InsBracketDeflection = zeros(numel(obj.Ftendon), 3);
obj.InsBracketReaction = zeros(numel(obj.Ftendon), 3);
if ~any(Cins(:))
    return
end
[N1, N2] = obj.gridSize();
nPts = size(obj.Location, 1);
utIns = utPre(:, :, :, 2);
dI = zeros(N1, 3, N2);
for ii = 1:N1
    for iii = 1:N2
        uT = reshape(utIns(ii, :, iii), 3, 1);
        Rv = obj.Ftendon(ii, iii) * uT;
        dI(ii, :, iii) = reshape(Cins * Rv, 1, 3);
        obj.InsBracketReaction(ii + (iii - 1) * N1, :) = reshape(Rv, 1, 3);
    end
end
obj.InsBracketDeflection = reshape(permute(dI, [1, 3, 2]), N1 * N2, 3);
obj.L_p(nPts, :, :, :) = obj.L_p(nPts, :, :, :) ...
    + reshape(permute(dI, [2, 1, 3]), [1, 3, N1, N2]);
end

%% ------------- xi3pose: per-cell wrapping loss (eq:x3) ------------
function dR = xi3RouteLoss(obj, strainXi0f, N)
% deltaRoute(i) = Xi3 * bend(i) * comp(i)^2 (MonoPam_mult.m Contraction_k
% line 449 port; dissertation 30-torque.tex eq:deltaL + eq:x3). comp uses
% the Xi0-only pre-solve strain (lines 445-447: comp = max(0, 1 -
% strain0/KMAX)); the loss grows the per-cell contraction demand when
% folded into s0 by the caller. bend(i) is the per-cell total bend arc
% (m): the supplied BendMeasure profile, or the self-computed sum of
% WrapRadius*|turn angle| over the crossing contacts (tendon-bend-only
% semantics: sleeves are straight pinned vessels, the tendon bends over
% the crossing guides - buildGeoExclusion.m geo.wRap lineage; App C
% "w_rap|theta_turn|, w_rap = 25 mm"; the arc is NOT added to the
% geometric centerline).
dR = zeros(N, 1);
obj.BendArc = nan(N, 2);
if ~(isscalar(obj.Xi3) && obj.Xi3 > 0)
    return
end
KMAX = (obj.RestingL - obj.Kmax) / obj.RestingL;
comp = max(0, 1 - strainXi0f / KMAX);
if ~isempty(obj.BendMeasure)
    bend = obj.BendMeasure;
    if isequal(size(bend), [obj.gridSize()])
        bend = reshape(bend, [], 1);
    elseif numel(bend) ~= N
        error('BiPam_pulley2_xi3pose:BendMeasureSize', ...
            ['bendMeasure (arg 17) must be [] for self-compute, an ' ...
             'N1xN2 grid, or an N = N1*N2 vector (ii fastest), got ' ...
             '%dx%d with N = %d.'], size(bend, 1), size(bend, 2), N)
    else
        bend = bend(:);
    end
    if ~all(isfinite(bend)) || any(bend < 0)
        error('BiPam_pulley2_xi3pose:BendMeasureValues', ...
            'bendMeasure must be finite and non-negative (m).')
    end
else
    bend = selfBendArcs(obj);   % N x 1 total over contacts
end
dR = obj.Xi3 .* bend .* comp.^2;
end

function bend = selfBendArcs(obj)
% Per-contact deviation arcs at the TENDON-side sheave/guide contacts of
% the routed polyline, per grid cell. Contact set (tendon-bend-only
% semantics - the sleeves are straight pinned pressure vessels, only the
% tendon bends over the tackle sheave and any downstream guides): the
% interior route vertices r with exitK <= r <= npts-1, where exitK =
% PulleyExitIndex of the ACTIVE crossing (the tackle exit; every vertex
% at/after it lies on the tendon span exit->insertion). The turn angle at
% vertex r is between the segment into r and the segment out of it,
% measured in vertex r's OWN class frame via the consistent parent chain
% (norms are frame-invariant; the wrapper's fr_bendProfile implements the
% identical construction and fr_xi3_check asserts they agree). The arcs
% are an ESTIMATE for the Xi3 bend loss only - never added to the
% geometric centerline (buildGeoExclusion.m geo.wRap lineage; App C
% "w_rap|theta_turn|, w_rap = 25 mm"). Degenerate segments contribute 0.
[N1, N2] = obj.gridSize();
N = N1 * N2;
bend = zeros(N, 1);
npts = size(obj.Location, 1);
exitK = obj.PulleyExitIndex(obj.ActiveCrossing);
contacts = exitK:(npts - 1);
nC = numel(contacts);
obj.BendArc = zeros(N, max(nC, 0));
if nC == 0
    return
end
for ci = 1:nC
    r = contacts(ci);
    arcs = zeros(N, 1);
    for ii = 1:N1
        for iii = 1:N2
            fr = obj.rowFrame(r);
            pIn = obj.pointInFrame(r - 1, fr, ii, iii);
            pX = obj.pointInFrame(r, fr, ii, iii);
            pOut = obj.pointInFrame(r + 1, fr, ii, iii);
            u1 = pX - pIn;
            u2 = pOut - pX;
            if norm(u1) < 1e-9 || norm(u2) < 1e-9
                arcs(ii + (iii - 1) * N1) = 0;
            else
                cs = dot(u1 / norm(u1), u2 / norm(u2));
                arcs(ii + (iii - 1) * N1) = ...
                    obj.WrapRadius * acos(min(1, max(-1, cs)));
            end
        end
    end
    obj.BendArc(:, ci) = arcs;
    bend = bend + arcs;
end
end

%% ------------- xi3pose: scalar config helper ----------------------
function val = cfgField2scalar(cfg, name, default)
if isfield(cfg, name) && ~isempty(cfg.(name))
    val = double(cfg.(name));
    if ~isscalar(val)
        val = val(1);
    end
else
    val = default;
end
end

function obj = applyExitDeflections2(obj, ubPre, utPre, activeK)
% Post-solve geometric application (variant 2): deflect the exit rows of
% L_p under the solved force state and store the deflection grids.
%   tackle mount (active crossing):  delta = Cmount * (nBPA*F*u_bpa +
%     Ftendon*u_t) -- the SAME reaction the base class reports as
%     ReactionF (Ftendon already slack-zeroed), i.e. the tackle-increased
%     internal force, with u_bpa referenced to the BPA muscle origin.
%   bowden anchor (each bowden crossing): delta = Cbowden * (tension*u_t)
%     with the crossing's own segment tension (nBPA*F/G at/below the
%     active tackle, nBPA*F above it).
% Displacements are mapped rotation-only from the crossing's child frame
% into the exit row's own frame (a displacement vector picks up no joint
% translation). Infeasible cells keep their deflection (their torque is
% NaN'd downstream); slack cells carry the BPA-side term only.
[N1, N2] = obj.gridSize();
nBPAact = obj.NPulleyBPA(activeK);
Gact = obj.PulleyGain(activeK);
Cm = obj.CmountEffective;
Cb = obj.CbowdenEffective;

% --- tackle-block mount deflection ---------------------------------
if any(obj.PulleyActive) && any(Cm(:))
    utA = utPre(:, :, :, activeK);
    dM = zeros(N1, 3, N2);
    for ii = 1:N1
        for iii = 1:N2
            Rv = (nBPAact * obj.FsingleBPA(ii,iii)) * ...
                reshape(ubPre(ii, :, iii), 3, 1) ...
                + obj.Ftendon(ii,iii) * reshape(utA(ii, :, iii), 3, 1);
            dM(ii, :, iii) = reshape(Cm * Rv, 1, 3);
        end
    end
    obj.MountDeflection(:, :, :, activeK) = dM;
    exR = obj.PulleyExitIndex(activeK);
    frEx = obj.rowFrame(exR);
    for ii = 1:N1
        for iii = 1:N2
            dLoc = obj.mapDirToFrame(reshape(dM(ii, :, iii), 1, 3), ...
                frEx, ii, iii, activeK + 1);
            obj.L_p(exR, :, ii, iii) = ...
                reshape(obj.L_p(exR, :, ii, iii), 1, 3) + dLoc;
        end
    end
end

% --- bowden anchor deflection --------------------------------------
for k = 1:2
    if ~strcmp(obj.RoutingMode{k}, 'bowden') || ~any(Cb(:))
        continue;
    end
    if k >= activeK
        tk = obj.Ftendon;               %single tendon line tension
    else
        tk = nBPAact * obj.FsingleBPA;  %full bundle pull above the tackle
    end
    utk = utPre(:, :, :, k);
    dB = zeros(N1, 3, N2);
    for ii = 1:N1
        for iii = 1:N2
            Rv = tk(ii,iii) * reshape(utk(ii, :, iii), 3, 1);
            dB(ii, :, iii) = reshape(Cb * Rv, 1, 3);
        end
    end
    obj.BowdenDeflection(:, :, :, k) = dB;
    exR = obj.PulleyExitIndex(k);
    frEx = obj.rowFrame(exR);
    for ii = 1:N1
        for iii = 1:N2
            dLoc = obj.mapDirToFrame(reshape(dB(ii, :, iii), 1, 3), ...
                frEx, ii, iii, k + 1);
            obj.L_p(exR, :, ii, iii) = ...
                reshape(obj.L_p(exR, :, ii, iii), 1, 3) + dLoc;
        end
    end
end
end

%% -------------Force unit direction ---------------
function F_unit = computeForceVectorBi(obj)
% Force unit direction from the origin (row 1, proximal frame) to the next
% non-duplicate point, chained into the proximal frame, per grid cell.
% Generalizes the base class's computeForceVector (its k == C branch).

L = obj.Location;
[N1, N2] = obj.gridSize();
N = N1 * N2;

F_unit = zeros(N, 3);
for ii = 1:N1
    for iii = 1:N2
        lin = ii + (iii - 1) * N1;
        pt1 = L(1, :);
        pt2 = pt1;
        found = false;
        for r = 2:size(L, 1)
            d = norm(L(r, :) - L(1, :));
            if d > 1e-6
                pt2 = obj.pointInFrame(r, 1, ii, iii);
                found = true;
                break;
            end
        end
        if ~found
            warning("Grid cell (%d,%d): No valid second point, using pt1=pt2", ii, iii);
        end
        F_unit(lin, :) = pt2 - pt1;
    end
end
F_unit = normalizeRows(F_unit);

end

%% ------------- Location / equilibrium ------------------------
function [LOCf, gema, sCon, Fmag, infeasible, slack, cbOut, dBr1, uHat1] = LokBi2(obj, X1, X2, kSpr, Funit, strain_predef, X0, DeltaL, s0, s0ref, cAdd) %#ok<INUSL,X0>
% VARIANT-2 twin of BiPam_pulley's LokBi (verbatim except the cAdd
% pass-through): the bracket (row 1) is deformed per grid cell under the
% BPA-side force and the per-cell scalar pulley equilibrium is solved,
% now WITH the mount/bowden deflection travel-absorption coefficient
% cAdd (m/N per unit single-BPA force; exactly zero when the variant-2
% compliances are off).
% xi3pose additions (derived-class diffs 5-6): the per-cell bracket
% deflection dBr1 (N x 3, bracket frame) and force direction uHat1
% (N x 3, bracket frame) are returned for the Bracket1* outputs; the
% bracket stiffness assignment comes from obj.BracketStiffness1
% (default [X1 X2 X2] = the base class's hard-coded K, so the default
% arithmetic is unchanged).

L = obj.Location;
rest = obj.RestingL;
Fm = obj.Fmax;
P = obj.Pressure;
D = obj.Diameter;
KMAX = (rest - obj.Kmax) / rest;
N = numel(strain_predef);

% Force magnitude at the Xi0-only pre-contraction, along the route
% direction (proximal frame).
relstrain = strain_predef / KMAX;
FF = festo4(D, relstrain, P) * Fm;
FF(FF < 0) = 0;
F = FF .* Funit;    % N x 3, proximal frame

% Bracket frame from the (orientation-free) origin point. Frame-count
% robustness: the base class hard-codes L(1,:,92); row 1 here lives in the
% proximal frame with no orientation dependence, so pA = L(1,:) directly.
pA = L(1, :);
switch obj.Diameter
    case 20
        Pbr = [9.48  -36.21   30.86] / 1000;   %bracket bolt pattern centroid
    case 10
        Pbr = [-19 22 27.6] / 1000;            %centroid of bracket cantilever
    otherwise
        Pbr = [0 0 0];
end

phbrA = pA - Pbr;
thetabrA = atan2(phbrA(2), phbrA(1));
RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...
    sin(thetabrA) cos(thetabrA) 0; ...
    0    0   1];
pbrhA = RhbrZ' * phbrA';
thetaY = atan2(pbrhA(3), pbrhA(1));
Ry = [cos(thetaY)  0  sin(thetaY);
    0            1  0;
    -sin(thetaY) 0  cos(thetaY)];
Rhbr = RhbrZ * Ry';
Thbr = RpToTrans(Rhbr, Pbr');

Fbrh = zeros(N, 3);
for lin = 1:N
    Fbrh(lin, :) = RowVecTrans(Thbr \ eye(4), F(lin, :));
end

if isinf(X1) && isinf(X2) && isinf(kSpr)
    [epsilon, delta, beta, gema, sCon, Fmag] = deal(zeros(N, 1));
    % DERIVED-CLASS FIX (diff 7): the base escape returns slack/infeasible
    % as DOUBLE zeros, which breaks the logical indexing downstream
    % (obj.Ftendon(obj.PulleySlack) = 0) - the base class was never
    % exercised with an infinite tendon rate (wraps = Inf) on this path.
    % Logical false is identical for every logical use and enables the
    % fully-rigid-tendon probes of fr_xi3_check.
    infeasible = false(N, 1);
    slack = false(N, 1);
    sCon = s0;      %rigid escape: no elastic absorption anywhere
    cbOut = zeros(N, 1);
    dBr1 = zeros(N, 3);
    uHat1 = zeros(N, 3);
else
    [epsilon, delta, beta, gema, sCon, Fmag, infeasible, slack, cbOut, ...
        uHat1] = ...
        fortzBi2(obj, Fbrh, X1, X2, kSpr, DeltaL, s0, s0ref, cAdd);
    dBr1 = [epsilon, delta, beta];
end
deflection = [epsilon, delta, beta];
pbrAnew = [norm(pbrhA), 0, 0] + deflection;

% Rebuild the deformed route as Npts x 3 x N (flattened), row 1 moved:
LOCf = zeros(size(L, 1), 3, N);
for lin = 1:N
    LOCf(:, :, lin) = L;
    LOCf(1, :, lin) = RowVecTrans(Thbr, pbrAnew(lin, :));
end

end

%% Pulley equilibrium per grid cell (variant-2 twin of fortzBi)
function [e_axial, e_bendY, e_bendZ, e_cable, sCon, Fmag, infeasible, slack, cbOut, uHat1] = fortzBi2(obj, Fbr, X1, X2, kSpr, DeltaL, s0, s0ref, cAdd) %#ok<INUSL>
% Same math as BiPam_pulley.fortzBi, generalized to the two-joint grid,
% plus the variant-2 deflection terms: solve s from
%   G*(s - (s0(i) - DeltaL(i)) - cb*F(s) - d_bpa(s))
%       = DeltaL + nBPA*F(s)/(G*kSpr) - d_tdn(s)
% with the mount/bowden deflection travel terms (LINEAR in F):
%   d_bpa + d_tdn-equivalent  =  cAdd(i) * F(s)   (see class header);
% i.e. psiFun = base psiFun + cAdd(i)*F(s). cAdd is EXACTLY zero when the
% variant-2 compliances are off, making psiFun bit-for-bit the base
% psiFun (IEEE: adding an exact 0 term changes no value). Slack and
% infeasibility clamps are unchanged. xi3pose: returns uHat1 (the
% per-cell force direction in the bracket frame) for the Bracket1*
% outputs, and the bracket stiffness uses obj.BracketStiffness1
% (default [X1 X2 X2] = the base hard-coded K, bit-identical).

N = size(Fbr, 1);
[e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N, 1));
[sCon, Fmag] = deal(zeros(N, 1));
infeasible = false(N, 1);
slack = false(N, 1);
cbOut = zeros(N, 1);
uHat1 = zeros(N, 3);

D = obj.Diameter;
rest = obj.RestingL;
mif = obj.Fmax;         %maximum force of ONE BPA
kmax = obj.Kmax;
KMAX = (rest - kmax) / rest;
P = obj.Pressure;
G = obj.PulleyGain(obj.ActiveCrossing);      %tackle gain
nBPA = obj.NPulleyBPA(obj.ActiveCrossing);   %parallel BPAs feeding the tackle

norms = vecnorm(Fbr, 2, 2);
valid = norms > 1e-3 & all(~isnan(Fbr), 2);
u_hat_all = normalizeRows(Fbr);

K = obj.BracketStiffness1;   %xi3pose: default [X1 X2 X2] = base hard-coded
C_bracket = diag([1 / K(1), 1 / K(2), 1 / K(3)]);
cSpr = 1 / kSpr;             %compliance of the ONE tendon line

% Force law at absolute contraction s for ONE BPA (festo4 returns >= 0 and
% 0 for relative strain above 1)
Ffun = @(s) festo4(D, (s / rest) / KMAX, P) * mif;

for i = 1:N
    if ~valid(i)
        sCon(i) = s0(i);    %base-identical escape: no deformation solved
        continue;
    end

    u_hat_i = u_hat_all(i, :);
    cb_i = (u_hat_i * C_bracket) * u_hat_i.';
    cbOut(i) = cb_i;
    uHat1(i, :) = u_hat_i;   %xi3pose: bracket-frame direction for outputs

    % Scalar equilibrium (see class header): the rigid reference is the
    % per-cell s0(i) - DeltaL(i), NOT s0ref -- on a general distal-tackle
    % route the pre-exit rows' BPA-side rigid change makes s0ref +
    % DeltaL(i) differ from s0(i), and referencing every cell to
    % s0ref would drop that change (the G = 1 limit would miss the
    % base-class solve). The two references coincide when the pre-exit
    % rows are orientation-free:
    %   G*(s - (s0(i) - DeltaL(i)) - cb*F(s) - d_bpa(s))
    %       = DeltaL(i) + nBPA*F(s)*cSpr/G - d_tdn(s)
    % with d_bpa/d_tdn folded into the linear coefficient cAdd(i):
    s0ref_i = s0(i) - DeltaL(i);
    psiFun = @(s) G * s - G * cb_i * Ffun(s) - nBPA * Ffun(s) * cSpr / G ...
        - DeltaL(i) - G * s0ref_i + cAdd(i) * Ffun(s);

    sLo = 0;                %contract per the contract's clamp
    sHi = KMAX * rest;      %fully contracted
    psiLo = psiFun(sLo);
    psiHi = psiFun(sHi);

    if psiLo > 0
        % Even at zero contraction the span demands less takeup than the
        % tackle would feed: the tendon goes slack. Clamp to s = 0.
        sCon(i) = sLo;
        Fmag(i) = Ffun(sLo);
        e_bkt = C_bracket * (Fmag(i) * u_hat_i.');
        e_axial(i) = e_bkt(1);
        e_bendY(i) = e_bkt(2);
        e_bendZ(i) = e_bkt(3);
        e_cable(i) = 0;     %slack tendon carries no stretch
        slack(i) = true;
    elseif psiHi < 0
        % Required contraction exceeds KMAX*Rest: infeasible cell. Clamp
        % to full contraction (force law gives ~0 there) and flag; the
        % torque is NaN'd by the caller.
        sCon(i) = sHi;
        Fmag(i) = Ffun(sHi);
        e_bkt = C_bracket * (Fmag(i) * u_hat_i.');
        e_axial(i) = e_bkt(1);
        e_bendY(i) = e_bkt(2);
        e_bendZ(i) = e_bkt(3);
        e_cable(i) = nBPA * Fmag(i) * cSpr / G;
        infeasible(i) = true;
    else
        try
            sCon(i) = fzero(psiFun, [sLo, sHi]);
        catch
            sCon(i) = sLo;   %mirror the base class's r = 0 escape
            slack(i) = true;
        end
        sCon(i) = max(sCon(i), 0); %guard against s being slightly negative.
        Fmag(i) = Ffun(sCon(i));

        if sCon(i) == 0 || Fmag(i) == 0
            continue;
        elseif isinf(X1) && isinf(X2)
            % Rigid body: no bracket deformation
            e_axial(i) = 0;
            e_bendY(i) = 0;
            e_bendZ(i) = 0;
            e_cable(i) = nBPA * Fmag(i) * cSpr / G;  %tendon-side stretch
            continue;
        else
            % Bracket displacement under the BPA-side force
            e_bkt = C_bracket * (Fmag(i) * u_hat_i.');
            e_axial(i) = e_bkt(1);
            e_bendY(i) = e_bkt(2);
            e_bendZ(i) = e_bkt(3);
            % Tendon stretch: the tackle divides the summed force by G
            % (base: F/kSpr).
            e_cable(i) = nBPA * Fmag(i) * cSpr / G;
        end
    end
end
end

%% tendon springrate
function springrate = Spr(obj, wraps)
% wraps: number of cable wraps around the post (affects effective stiffness)
switch obj.Diameter
    case 10
        mult = 2;
    case 20
        mult = 6;
    otherwise
        mult = 2;
end

if nargin < 2 || isempty(wraps)
    wraps = mult;
end

Aeff = 1.51*10^-6;%Effective area for 19-strand cable
E = 193*10^9;       %Young's Modulus
L = obj.TendonL;      %tendon length

springrate = wraps*Aeff*E/L;
end

%% ------------- Segment lengths (deformed, flattened) ---------
function SL = segBi(obj)
[N1, N2] = obj.gridSize();
N = N1 * N2;
M = size(obj.Location, 1);
SL = zeros(N, M - 1);

for ii = 1:N1
    for iii = 1:N2
        lin = ii + (iii - 1) * N1;
        for i = 1:M - 1
            pA = obj.pointInFrameDeformed(i, 1, ii, iii);
            pB = obj.pointInFrameDeformed(i + 1, 1, ii, iii);
            SL(lin, i) = norm(pA - pB);
        end
    end
end
end

%% Grid-cell masking helpers (N1 x 3 x N2 arrays, N1 x N2 masks)
function out = maskNaN3(in, mask)
m3 = repmat(reshape(logical(mask), size(in, 1), 1, []), [1, size(in, 2), 1]);
out = in;
out(m3) = NaN;
end

function out = maskZero3(in, mask)
m3 = repmat(reshape(logical(mask), size(in, 1), 1, []), [1, size(in, 2), 1]);
out = in;
out(m3) = 0;
end

%% Row normalization over trailing 3-D (N1 x 3 x N2 -> unit rows)
function vhat = normalize3(v)
% Per-cell normalization of N1 x 3 x N2 direction arrays. The 3-wide
% component axis is dim 2, so the per-cell norm is vecnorm along dim 2
% (an N1 x 1 x N2 grid) and the division uses implicit expansion.
% IMPORTANT: never reshape(v, [], 3) here -- dim 1 (grid cells) is the
% fastest-varying axis of an N1 x 3 x N2 array, so consecutive linear
% elements are DIFFERENT cells' x-components, not one cell's vector;
% a reshape-based normalization silently scrambles the components.
n = vecnorm(v, 2, 2);
valid = n > 1e-9 & all(~isnan(v), 2);
nsafe = n;
nsafe(~valid) = 1;
vhat = v ./ nsafe;
bad = repmat(~valid, [1, size(v, 2), 1]);
vhat(bad) = 0;
end

%% Row normalization for N x 3 matrices
function vhat = normalizeRows(v)
norms = vecnorm(v, 2, 2);
valid = norms > 1e-3 & all(~isnan(v), 2);
vhat = zeros(size(v));
vhat(valid, :) = v(valid, :) ./ norms(valid);
end

%% Row-vector DIRECTION transform (rotation block only)
function v = RowVecDir(T, v)
% Rotates a row direction vector between frames: applies ONLY the
% rotation block of the homogeneous transform, never its translation
% (RowVecTrans is the point-transform counterpart).
R = T(1:3, 1:3);
v = (R * v.').';
end
