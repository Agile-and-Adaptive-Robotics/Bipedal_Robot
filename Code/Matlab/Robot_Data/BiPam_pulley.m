% Pam Data -- BIARTICULAR reverse-pulley (block-and-tackle) BPA transmission
% Author: Ben Bolen
% Date: 2026-09-24
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

classdef BiPam_pulley < handle

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

        % --- Reverse-pulley transmission config (per crossing) ---
        NPulleyBPA                  %1x2 BPAs in parallel feeding crossing k's tackle
        TackleLineParts             %1x2 rope parts supporting crossing k's moving block
        PulleyGain                  %1x2 travel gain G (default = TackleLineParts)
        RoutingMode                 %1x2 cell, 'moving_via' | 'bowden' per crossing
        PulleyExitIndex             %1x2 tendon exit row per crossing (default Cross-1)
        BowdenBossDia               %Bowden barrel-adjuster boss diameter, m
        BowdenRunClearance          %Bowden housing run clearance, m

        % --- Stiffness-aware fields (minimizeFlx-style) ---
        L_p                         %Npts x 3 x N1 x N2 deformed locations (bracket moved row 1)
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
        function PD = BiPam_pulley(name, location, cross, diameter, t, ...
            rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, wraps, pulleyConfig)
            if nargin == 14
                pulleyConfig = struct();
            end
            if nargin == 14 || nargin == 15
                if ~isstruct(pulleyConfig)
                    error('BiPam_pulley:PulleyConfig', ...
                        ['The 15th input must be the pulley config struct ' ...
                         '(nPulleyBPA, tackleLineParts, gain, routingMode, ' ...
                         'pulleyExitIndex, bowdenBossDia, bowdenRunClearance; ' ...
                         'scalar or 1x2 per crossing).'])
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

                for k = 1:2
                    if ~(isscalar(PD.NPulleyBPA(k)) && isfinite(PD.NPulleyBPA(k)) ...
                            && PD.NPulleyBPA(k) >= 1 && ...
                            PD.NPulleyBPA(k) == round(PD.NPulleyBPA(k)))
                        error('BiPam_pulley:NPulleyBPA', ...
                            ['nPulleyBPA(%d) must be a positive integer ' ...
                             '(BPAs in parallel feeding the tackle), got %g.'], ...
                            k, PD.NPulleyBPA(k))
                    end
                    if ~(isscalar(PD.TackleLineParts(k)) && ...
                            isfinite(PD.TackleLineParts(k)) && ...
                            PD.TackleLineParts(k) >= 1)
                        error('BiPam_pulley:TackleLineParts', ...
                            ['tackleLineParts(%d) must be >= 1 (rope parts ' ...
                             'supporting the moving block = travel gain G; ' ...
                             'a four-wheel tackle rigges to at most 4), got %g.'], ...
                            k, PD.TackleLineParts(k))
                    end
                    if ~(isscalar(PD.PulleyGain(k)) && isfinite(PD.PulleyGain(k)) ...
                            && PD.PulleyGain(k) >= 1)
                        error('BiPam_pulley:GainRange', ...
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
                        error('BiPam_pulley:ExitIndex', ...
                            ['PulleyExitIndex(%d) must name a row before ' ...
                             'Cross(%d) (1..%d), got %d.'], ...
                            k, k, PD.Cross(k) - 1, PD.PulleyExitIndex(k))
                    end
                end
                if PD.PulleyGain(1) > 1 && PD.PulleyGain(2) > 1
                    error('BiPam_pulley:TwoTackles', ...
                        ['Both crossings carry gain > 1 (%g, %g), but only ' ...
                         'ONE tackle is supported: the crossing-1 mount ' ...
                         'reaction and the shared ReactionFmag grid assume ' ...
                         'the distal crossing is the only active one. Set ' ...
                         'the other crossing''s gain to 1.'], ...
                        PD.PulleyGain(1), PD.PulleyGain(2))
                end
                if numel(cross) ~= 2 || cross(1) < 2 || cross(2) < cross(1) ...
                        || cross(2) > size(location, 1)
                    error('BiPam_pulley:CrossShape', ...
                        ['Cross must be [c1 c2] with 2 <= c1 <= c2 <= ' ...
                         'npts, got [%d %d] with %d points.'], ...
                        cross(1), cross(2), size(location, 1))
                end

                % Automatically compute stiffness-aware geometry and torque
                PD = PD.updateBiPulleyGeometry();
            else
                fprintf('Invalid number of arguments\n')
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

            % Deformed geometry, tendon stretch, and the pulley
            % equilibrium (single route state; see class header).
            [L_p_f, gama_f, sCon_f, Fmag_f, infeas_f, slack_f, cb_f] = LokBi( ...
                obj, obj.Xi1, obj.Xi2, obj.kSpr, obj.Funit, ...
                strainXi0f, obj.Xi0, deltaLf, s0f, s0ref);

            % Deformed locations back onto the 4-D grid.
            obj.L_p = reshape(L_p_f, [size(obj.Location, 1), 3, N1, N2]);
            obj.gama = gama_f;
            obj.sContraction = reshape(sCon_f, N1, N2);
            obj.FsingleBPA = reshape(Fmag_f, N1, N2);

            % Tackle tension and travel at the active crossing
            % (MonoPam_pulley conventions).
            Gact = obj.PulleyGain(activeK);
            nBPAact = obj.NPulleyBPA(activeK);
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

            % Deformed segment lengths and musculotendon length.
            obj.sL_p = segBi(obj);
            Lmt_p_f = sum(obj.sL_p, 2) - obj.Xi0;
            obj.Lmt_p = Lmt_p_f;

            % Route contraction with bracket deformation, tendon stretch,
            % and Xi0 (base-identical bookkeeping).
            obj.strain_p = (obj.RestingL - (Lmt_p_f - obj.TendonL - gama_f ...
                - 2 * obj.FittingLength)) ./ obj.RestingL;

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
function [LOCf, gema, sCon, Fmag, infeasible, slack, cbOut] = LokBi(obj, X1, X2, kSpr, Funit, strain_predef, X0, DeltaL, s0, s0ref) %#ok<INUSL,X0>
% Flattened-grid version of MonoPam_pulley's Lok: the bracket (row 1) is
% deformed per grid cell under the BPA-side force and the per-cell scalar
% pulley equilibrium is solved. Row 1 is orientation-free (proximal frame),
% so the bracket-frame construction is done once.

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
    -sin(thetaY) 0   cos(thetaY)];
Rhbr = RhbrZ * Ry';
Thbr = RpToTrans(Rhbr, Pbr');

Fbrh = zeros(N, 3);
for lin = 1:N
    Fbrh(lin, :) = RowVecTrans(Thbr \ eye(4), F(lin, :));
end

if isinf(X1) && isinf(X2) && isinf(kSpr)
    [epsilon, delta, beta, gema, sCon, Fmag, infeasible, slack] = ...
        deal(zeros(N, 1));
    sCon = s0;      %rigid escape: no elastic absorption anywhere
    cbOut = zeros(N, 1);
else
    [epsilon, delta, beta, gema, sCon, Fmag, infeasible, slack, cbOut] = ...
        fortzBi(obj, Fbrh, X1, X2, kSpr, DeltaL, s0, s0ref);
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

%% Pulley equilibrium per grid cell (MonoPam_pulley's fortz, flattened)
function [e_axial, e_bendY, e_bendZ, e_cable, sCon, Fmag, infeasible, slack, cbOut] = fortzBi(obj, Fbr, X1, X2, kSpr, DeltaL, s0, s0ref) %#ok<INUSL>
% Same math as MonoPam_pulley.fortz, generalized to the two-joint grid:
% solve s from
%   G*(s - (s0(i) - DeltaL(i)) - cb*F(s)) = DeltaL + nBPA*F(s)/(G*kSpr)
% per grid cell, with slack/infeasible clamps. The per-cell reference
% s0(i) - DeltaL(i) equals s0ref exactly when the pre-exit rows are
% orientation-free (then this is verbatim MonoPam_pulley's solve) and
% carries the BPA-side rigid change otherwise (class header). G and nBPA
% are the ACTIVE crossing's config
% (obj.PulleyGain/obj.NPulleyBPA(obj.ActiveCrossing)).

N = size(Fbr, 1);
[e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N, 1));
[sCon, Fmag] = deal(zeros(N, 1));
infeasible = false(N, 1);
slack = false(N, 1);
cbOut = zeros(N, 1);

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

K = [X1, X2, X2];
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

    % Scalar equilibrium (see class header): the rigid reference is the
    % per-cell s0(i) - DeltaL(i), NOT s0ref -- on a general distal-tackle
    % route the pre-exit rows' BPA-side rigid change makes s0ref +
    % DeltaL(i) differ from s0(i), and referencing every cell to
    % s0ref would drop that change (the G = 1 limit would miss the
    % base-class solve). The two references coincide when the pre-exit
    % rows are orientation-free:
    %   G*(s - (s0(i) - DeltaL(i)) - cb*F(s))
    %       = DeltaL(i) + nBPA*F(s)*cSpr/G
    s0ref_i = s0(i) - DeltaL(i);
    psiFun = @(s) G * s - G * cb_i * Ffun(s) - nBPA * Ffun(s) * cSpr / G ...
        - DeltaL(i) - G * s0ref_i;

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
