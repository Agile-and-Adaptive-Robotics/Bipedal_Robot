% Pam Data -- reverse-pulley (block-and-tackle) BPA transmission
% Author: Ben Bolen
% Date: 2026-09-24
% Description: Spin-off of MonoPamDataExplicit_balance for Ben's
% block-and-tackle rig (2026-09-24): nPulleyBPA 20 mm BPAs in parallel
% (the physical rig used 2) drive a four-wheel tackle mounted with them on
% the proximal (foot) body; the artificial tendon exits the tackle, crosses
% the joint, and inserts on the distal (tibia) body. Effect: insertion
% travel is multiplied by the tackle gain G at the cost of a proportional
% loss of force transfer (insertion force = total BPA force / G), and the
% reaction force on the foot mount grows (approximately F_BPA + F_tendon).
% G = 1 with nPulleyBPA = 1 means "pulley disabled": the class then
% reproduces MonoPamDataExplicit_balance to ~1e-8 relative on F_p, mA_p,
% and Torque_p for a geometry whose crossed segment is the tendon line
% (the regression identity asserted by Opt_sanity_pulley.m).
%
% Transmission model (Ben's hardware contract, implemented here):
%   theta0     reference orientation = the pose with the LONGEST rigid
%              tendon span (the rope's unstretched/installation pose), so
%              DeltaL >= 0 at every other pose for a flexor route.
%   s          solved ABSOLUTE BPA contraction travel from the rest length
%              (s/RestingL is the contraction fraction fed to festo4);
%              stored as sContraction.
%   DeltaL(i)  Lspan(theta0) - Lspan(i): rigid span change of the TENDON
%              segment (row PulleyExitIndex -> row Cross), computed from
%              the undeformed routing points.
%   delta_t    tendon stretch, delta_t = F_t/kSpr with F_t the tension.
%   Closure    in tackle-travel form (the contract's literal form)
%                 G*sTravel = DeltaL + delta_t,
%              where sTravel = s - s0ref - cb*F is the BPA-side travel the
%              tackle actually receives: s0ref is the rigid-geometry
%              contraction already present at the installation pose
%              theta0 (the rope is rigged with the BPA pre-tensioned
%              there), and cb*F is the share of the stroke the BPA-side
%              bracket compliance absorbs before the tackle sees it.
%              In terms of the absolute contraction the closure is
%                 G*(s - s0ref - cb*F(s)) = DeltaL + F_total(s)/(G*kSpr)
%              with F_total = nPulleyBPA*F_single(s) the summed parallel-
%              BPA force at contraction s. The bracket share and the
%              installation offset are REQUIRED by the G = 1 regression
%              identity: at G = 1, nPulleyBPA = 1 the equilibrium reduces
%              algebraically to MonoPamDataExplicit_balance's solve
%              r = s - s0(i) = (cb + 1/kSpr)*F, because for any route
%              whose non-tendon rows are orientation-constant
%              DeltaL(i) = s0(i) - s0ref. Opt_sanity_pulley.m asserts both
%              the identity and the literal closure residual.
%   nPulleyBPA BPAs IN PARALLEL feeding the tackle: total pull
%              nPulleyBPA*F_single, travel gain unchanged (G only),
%              F_t = nPulleyBPA*F_single/G, and the mount reaction sums
%              the nPulleyBPA BPA lines plus the tendon line. Each BPA
%              keeps its own bracket, so the per-BPA stroke share cb*F
%              (not nPulleyBPA*cb*F) enters the closure. The parallel BPAs
%              are assumed identical, matching the
%              biPulleySpecsFromOpenSim placement convention: bundle
%              symmetric about the original OpenSim line so the first-
%              order line of action and moment arm are preserved --
%              ASYMMETRIC placement shifts the line of action.
%   Solve      single scalar root find in s per orientation (fzero,
%              bracketed to [0, KMAX*Rest]; psi is monotone because the
%              festo4 force law does not increase with contraction). If
%              the required contraction exceeds KMAX*Rest the frame is
%              INFEASIBLE: PulleyInfeasible is set, s clamps to KMAX*Rest,
%              and the torque is NaN, as the base class does for strain
%              out of range. NOTE the direction: by the closure the
%              required contraction (DeltaL + delta_t)/G DECREASES with G,
%              so infeasibility is reached when DeltaL exceeds
%              G*(KMAX*Rest - s0ref), not by making G large.
%   Bracket    the Xi0/Xi1/Xi2 mount compliance (Lok/fortz machinery)
%              deforms the mount under the BPA-side force exactly as
%              MonoPamDataExplicit_balance does: the deformed origin
%              point A (row 1) drives the deformed route bookkeeping
%              (strain_p, F_p, mA_p, Torque_p), which is kept IDENTICAL
%              to the base class so the G = 1 regression is structural.
%   Routing    'moving_via' (default): the exit point is body-fixed to
%              the proximal body; the tendon-side unit direction u_t is
%              recomputed per orientation from the (deformed) exit point
%              to the insertion, so the moving exit/wrap point is honored
%              and the moment arm can INCREASE as the joint rotates (the
%              dissertation's "more daunting calculation").
%              'bowden': the tendon runs in a housing whose anchor
%              bracket sits on the distal (tibia) body, so u_t is FIXED
%              in the tibia frame (held at its theta0 value) and no
%              longer rotates with the joint; the closure/force-balance
%              math is unchanged. Packaging cost is the explicit
%              min-clearance envelope carried in BowdenBossDia /
%              BowdenRunClearance (standard Shimano-type parts: M7 x 1.0
%              mm barrel adjuster ~7 mm boss, 5 mm brake housing, 5 mm
%              bracket plate); Opt_run_pulley reports it.
%
% Outputs exposed per the contract: PulleyGain, F_ins, Torque_ins,
% ReactionF (with ReactionFmag), PulleyInfeasible -- plus the solved
% equilibrium state (sContraction, PulleyTravel, Ftendon, FsingleBPA,
% gama = delta_t, deltaL, RefIndex, PulleySlack) and the insertion-side
% moment arm mA_ins. Torque about the joint follows the base class
% moment-arm formula using the tendon line (exit point as the effective
% via point).
%
% Frame-count robustness: the base class hard-codes pA = L(1,:,92) in its
% Lok helper; this class derives the reference frame from
% size(Location,3) (iRef = min(92, N), which is exactly 92 for the base
% class's 92-orientation sweeps and the last frame for shorter sweeps).
%
%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPam_pulley < handle

    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which row corresponds with a location where the muscle crosses into a new reference frame
        Diameter                    %Diameter of the BPA
        TransformationMat           %Contains a transformation matrix to change the
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

        % --- Reverse-pulley transmission config (Ben, 2026-09-24) ---
        NPulleyBPA                  %BPAs in parallel feeding the tackle (rig used 2)
        TackleLineParts             %Rope parts supporting the moving block (integer, hardware description)
        PulleyGain                  %Travel gain G used in the math (= TackleLineParts unless overridden)
        RoutingMode                 %'moving_via' (exit rotates with the joint) | 'bowden' (u_t fixed in tibia frame)
        PulleyExitIndex             %Row the tendon leaves the tackle from (proximal body)
        BowdenBossDia               %Bowden barrel-adjuster boss diameter, m (min-clearance envelope)
        BowdenRunClearance          %Bowden housing run clearance, m (min-clearance envelope)

        % --- Stiffness-aware fields (minimizeFlx-style, base-identical) ---
        L_p                         %Deformed location matrix (updated attachment points)
        Funit                       %Force direction in hip frame
        sL_p                        %Segment lengths
        Lmt_p                       %Musculotendon length with deformed geometry and constant length offset
        uD_p                        %Force unit direction in tibia frame with deformed geometry
        strain_p                    %Contraction with bracket deformation, tendon stretch, and constant length offset
        F_p                         %Force vector with stiffness effects
        mA_p                        %Moment arm with stiffness effects
        Torque_p                    %Torque with stiffness effects
        gama                        %Tendon stretch delta_t = F_t/kSpr (cable elongation)
        kSpr                        %Tendon spring rate (effective, ONE tendon line)

        % --- Pulley-solved fields ---
        deltaL                      %Nx1 rigid tendon-span change from theta0
        RefIndex                    %theta0 frame index (longest rigid tendon span)
        sContraction                %Nx1 solved absolute BPA contraction travel from rest, m
        PulleyTravel                %Nx1 tackle input travel sTravel = s - s0ref - cb*F, m
        FsingleBPA                  %Nx1 force of ONE BPA at the solved contraction, N
        Ftendon                     %Nx1 tendon tension F_t = nPulleyBPA*F_single/G, N
        u_t                         %Nx3 tendon unit direction at the insertion, tibia frame
        u_bpa                       %Nx3 total BPA force unit direction, tibia frame
        F_ins                       %Nx3 insertion force vector = F_t * u_t (tibia frame)
        mA_ins                      %Nx3 moment arm of the tendon line about the joint
        Torque_ins                  %Nx3 torque from the insertion force (NaN where infeasible)
        ReactionF                   %Nx3 mount reaction = nBPA*F_BPA*u_bpa + F_t*u_t (tibia frame)
        ReactionFmag                %Nx1 magnitude of ReactionF
        PulleyInfeasible            %Nx1 logical, required contraction exceeded KMAX*Rest
        PulleySlack                 %Nx1 logical, tendon slack (solve clamped to s = 0)
    end

    % Dependent properties are calculated from the explicit properties.
    % MATLAB evaluates each get method when the property is queried.
    properties (Dependent)
        SegmentLengths
        LongestSegment
        MuscleLength
        Contraction
        LengthCheck
        UnitDirection
        MomentArm
        Fmax
        Force
        Torque
    end


    methods
        %% ------------- Muscle Data Constructor -----------------
        % Constructor function. The first fourteen arguments mirror
        % MonoPamDataExplicit_balance; the fifteenth is the pulley config
        % struct (fields nPulleyBPA, tackleLineParts, gain, routingMode,
        % pulleyExitIndex, bowdenBossDia, bowdenRunClearance -- all
        % optional). Omitting it configures a straight-tendon pulley-
        % disabled transmission (nPulleyBPA = 1, G = 1, 'moving_via',
        % exit row = Cross-1, standard Shimano-type Bowden envelope).
        function PD = MonoPam_pulley(name, location, cross, diameter, t, rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, wraps, pulleyConfig)
            if nargin == 14
                pulleyConfig = struct();
            end
            if nargin == 14 || nargin == 15
                if ~isstruct(pulleyConfig)
                    error('MonoPam_pulley:PulleyConfig', ...
                        ['The 15th input must be the pulley config struct ' ...
                         '(nPulleyBPA, tackleLineParts, gain, routingMode, ' ...
                         'pulleyExitIndex, bowdenBossDia, bowdenRunClearance).'])
                end
                PD.Name = name;                   % BPA/muscle name
                PD.Location = location;           % routing-point array
                PD.Cross = cross;                 % first row in the next frame
                PD.Diameter = diameter;           % BPA diameter, mm
                PD.TransformationMat = t;         % frame transforms
                PD.RestingL = rest;               % BPA resting length, m
                PD.Kmax = kmax;                   % fully contracted length, m
                PD.TendonL = tendon;              % tendon length, m
                PD.FittingLength = fitn;           % one fitting length, m
                PD.Pressure = pres;               % BPA pressure, kPa

                PD.Xi0 = xi0;                     % constant length offset, m
                PD.Xi1 = xi1;                     % axial bracket stiffness, N/m
                PD.Xi2 = xi2;                     % bending stiffness, N/m
                PD.Wraps = wraps;                 % tendon wrap count

                % Pulley config with defaults (1 = pulley disabled).
                PD.NPulleyBPA = cfgField(pulleyConfig, 'nPulleyBPA', 1);
                PD.TackleLineParts = cfgField(pulleyConfig, ...
                    'tackleLineParts', 1);
                PD.PulleyGain = cfgField(pulleyConfig, ...
                    'gain', PD.TackleLineParts);
                PD.RoutingMode = cfgField(pulleyConfig, ...
                    'routingMode', 'moving_via');
                PD.PulleyExitIndex = cfgField(pulleyConfig, ...
                    'pulleyExitIndex', cross - 1);
                PD.BowdenBossDia = cfgField(pulleyConfig, ...
                    'bowdenBossDia', 0.007);   % M7 x 1.0 mm barrel adjuster
                PD.BowdenRunClearance = cfgField(pulleyConfig, ...
                    'bowdenRunClearance', 0.005);  % 5 mm housing + plate

                if ~(isscalar(PD.NPulleyBPA) && isfinite(PD.NPulleyBPA) ...
                        && PD.NPulleyBPA >= 1 && PD.NPulleyBPA == round(PD.NPulleyBPA))
                    error('MonoPam_pulley:NPulleyBPA', ...
                        ['nPulleyBPA must be a positive integer (BPAs in ' ...
                         'parallel feeding the tackle), got %g.'], PD.NPulleyBPA)
                end
                if ~(isscalar(PD.TackleLineParts) && isfinite(PD.TackleLineParts) ...
                        && PD.TackleLineParts >= 1)
                    error('MonoPam_pulley:TackleLineParts', ...
                        ['tackleLineParts must be >= 1 (rope parts ' ...
                         'supporting the moving block = travel gain G; ' ...
                         'a four-wheel tackle rigges to at most 4), got %g.'], ...
                        PD.TackleLineParts)
                end
                if ~(isscalar(PD.PulleyGain) && isfinite(PD.PulleyGain) ...
                        && PD.PulleyGain >= 1)
                    error('MonoPam_pulley:GainRange', ...
                        'PulleyGain must be >= 1 (1 = pulley disabled), got %g.', ...
                        PD.PulleyGain)
                end
                PD.RoutingMode = validatestring(PD.RoutingMode, ...
                    {'moving_via', 'bowden'}, mfilename, 'routingMode');
                if ~(isscalar(PD.PulleyExitIndex) && PD.PulleyExitIndex >= 1 ...
                        && PD.PulleyExitIndex == round(PD.PulleyExitIndex) ...
                        && PD.PulleyExitIndex < PD.Cross)
                    error('MonoPam_pulley:ExitIndex', ...
                        ['PulleyExitIndex must name a proximal-body row ' ...
                         'before Cross (1..Cross-1), got %d with Cross = %d.'], ...
                        PD.PulleyExitIndex, PD.Cross)
                end

                % Automatically compute stiffness-aware geometry and torque
                PD = PD.updatePulleyGeometry();
            else
                fprintf('Invalid number of arguments\n')
            end
        end


        %% ------------- Segment Lengths ------------------------
        function segLengths = get.SegmentLengths(obj)
            L = obj.Location;
            C = obj.Cross;
            T = obj.TransformationMat;
            segLengths = zeros(size(T, 3), size(L, 1) - 1);

            for ii = 1:size(T, 3)                          %Repeat for each orientation
                for i = 1:size(L, 1)-1                      %Repeat for all muscle segments
                    pointA = L(i, :,ii);
                    pointB = L(i+1, :,ii);
                    if i+1 == C
                        pointB = RowVecTrans(T(:, :, ii), pointB);
                    end
                    segLengths(ii, i) = norm(pointA - pointB);
                end
            end
        end

        %% -------------- Longest Segment Calculation ----------------
        function longestSeg = get.LongestSegment(obj)
            L = obj.Location;
            segLengths = obj.SegmentLengths;

            % Calculate which muscle segment is the longest on average.
            % This will be where the Pam resides.
            avgSegL = zeros(size(L, 1) - 1);
            for i = 1:size(segLengths, 2)
                avgSegL(i) = mean(segLengths(:, i));
            end

            longestSegPointer = 1;
            if size(avgSegL, 1) > 1
                for i = 1:size(avgSegL, 1) - 1
                    if avgSegL(i + 1) > avgSegL(i)
                        longestSegPointer = i + 1;
                    end
                end
            end
            longestSeg = segLengths(:, longestSegPointer);
        end

        %% ------------- Muscle Length ------------------------
        %Function that calculates the muscle length, based
        function mL = get.MuscleLength(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            mL = zeros(size(T, 3), 1);
            segLengths = obj.SegmentLengths;

            for ii = 1:size(mL, 1)                          %Repeat for each orientation
                for i = 1:size(L, 1)-1                      %Repeat for all muscle segments
                    mL(ii, 1) = mL(ii, 1) + segLengths(ii, i);
                end
            end
        end

        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
        function unitD = get.UnitDirection(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            C = obj.Cross;
            direction = zeros(size(T, 3), 3);
            unitD = zeros(size(direction));

            for i = 1:size(T, 3)
                pointA = L(C-1, :, i);
                pointB = L(C, :, i);
                direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
                unitD(i, :) = direction(i, :)/norm(direction(i, :));
            end
        end

        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
        %For every ViaPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
        function mA = get.MomentArm(obj)
            T = obj.TransformationMat;
            L = obj.Location;
            C = obj.Cross;
            unitD = obj.UnitDirection;
            mA = zeros(size(T, 3), 3);

            for i = 1:size(T, 3)
                pointB = L(C, :, i);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
                %mA(i, :) = cross(pointB, unitD(i, :));
            end
        end

        %% -------------- Contraction of the PAM --------------------------
        function contraction = get.Contraction(obj)
            mL = obj.MuscleLength;
            rest = obj.RestingL;
            tendon = obj.TendonL;
            fitting = obj.FittingLength;

            contraction = (rest-(mL-tendon-2.*fitting))./rest;
        end

        %% -------------- Length Check --------------------------
        function lengthCheck = get.LengthCheck(obj)
            contraction = obj.Contraction;
            maxContractPercent = 0.25;          %Contracting to 75% of length
            minContractPercent = -0.03;          %Elongating to 103% of length
            restingPamLength = obj.RestingL;

            if restingPamLength < 0
                lengthCheck = 'Unusable';
            else
                if max(contraction) <= maxContractPercent
                    if min(contraction) >= minContractPercent
                        lengthCheck = 'Usable';
                    else
                        lengthCheck = 'Unusable';
                    end
                else
                    lengthCheck = 'Unusable';
                end
            end
        end

        %% -------------- Maximum Force --------------------------
        %Calculate the direction of the forced applied by the muscle
        function maxF = get.Fmax(obj)
        %Inputs:
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %Outputs:
        %maxF == Maximum Force, N, produced by BPA at 0% contraction and
        %           620 kPa

            dia = obj.Diameter;
            rest = obj.RestingL;

            if dia == 10
                maxF = maxBPAforce(rest,'10');
            elseif dia ==20
                maxF = maxBPAforce(rest,'20');
            elseif dia ==40
                maxF = maxBPAforce(rest,'40');
            else
                disp('Wrong size diameter BPA')
            end
        end

%% -------------- Force --------------------------
% Calculate PAM force using rigid/original geometry.
% Stiffness-adjusted force is handled separately by obj.F_p / obj.Torque_p.
function F = get.Force(obj)

    unitD = obj.UnitDirection;      % Nx3
    strain = obj.Contraction;       % Nx1

    rest = obj.RestingL;
    kmax = obj.Kmax;
    KMAX = (rest - kmax)/rest;      % max contraction fraction

    rel = strain ./ KMAX;           % relative strain

    Fn = festo4(obj.Diameter, rel, obj.Pressure);  % normalized force

    scalarForce = Fn .* obj.Fmax;   % N

    scalarForce(scalarForce < 0) = 0;
    scalarForce(scalarForce > obj.Fmax) = NaN;

    F = scalarForce .* unitD;       % Nx3

end
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
        function tor = get.Torque(obj)
            mA = obj.MomentArm;
            F = obj.Force;
            tor = cross(mA,F,2);
        end

        % ============================================================
        % === Stiffness-aware pipeline with the reverse pulley =======
        % ============================================================

        % Xi0: constant length offset
        % Xi1, Xi2: bracket stiffness components
        % Wraps: number of cable wraps (affects tendon spring rate)
        % nPulleyBPA/tackleLineParts: the block-and-tackle config
        function obj = updatePulleyGeometry(obj)
            % Tendon spring rate: ONE physical tendon line exits the
            % tackle, so kSpr is NOT multiplied by nPulleyBPA (the parallel
            % BPAs share the single output line).
            obj.kSpr = Spr(obj, obj.Wraps);

            % Force unit vector in hip frame (origin to first non-duplicate point)
            Funit_i = computeForceVector(obj);
            obj.Funit = Funit_i;

            % Rigid tendon span per orientation: the crossed segment, exit
            % row -> insertion row, measured in ONE frame (the tibia frame
            % for the exit; norm is frame-invariant).
            N = size(obj.Location, 3);
            E = obj.PulleyExitIndex;
            C = obj.Cross;
            span = zeros(N, 1);
            for ii = 1:N
                exitT1 = RowVecTrans(obj.TransformationMat(:,:,ii)\eye(4), ...
                    obj.Location(E,:,ii));
                span(ii) = norm(exitT1 - obj.Location(C,:,ii));
            end

            % Reference orientation theta0 = longest rigid tendon span (the
            % rope's unstretched/installation pose), so DeltaL >= 0 at the
            % other poses for a flexor route.
            [~, obj.RefIndex] = max(span);
            obj.deltaL = span(obj.RefIndex) - span;

            % Rigid-geometry contraction per frame and at theta0
            s0 = obj.RestingL - obj.MuscleLength + obj.Xi0 ...
                + obj.TendonL + 2*obj.FittingLength;

            % Contraction from constant length offset only (base-identical
            % pre-force used for the bracket-frame direction mask)
            strain_Xi0 = Contraction_k(obj, [], [], obj.Xi0);

            % Deformed geometry, tendon stretch, and the pulley equilibrium
            [L_p_i, gama_i, sCon_i, Fmag_i, infeas_i, slack_i, cb_i] = Lok( ...
                obj, obj.Xi1, obj.Xi2, obj.kSpr, ...
                obj.Funit, strain_Xi0, obj.Xi0, ...
                obj.deltaL, s0, s0(obj.RefIndex));
            obj.L_p = L_p_i;
            obj.gama = gama_i;
            obj.sContraction = sCon_i;
            obj.FsingleBPA = Fmag_i;

            % Unit direction with deformed geometry
            uD_p_i = UD(obj, obj.L_p);
            obj.uD_p = uD_p_i;

            % Segment lengths with deformed geometry
            sL_p_i = seg(obj, obj.L_p);
            obj.sL_p = sL_p_i;

            % Musculotendon length with deformed geometry and Xi0
            Lmt_p_i = LMT(obj.sL_p, obj.Xi0);
            obj.Lmt_p = Lmt_p_i;

            % Contraction with bracket deformation, tendon stretch, and Xi0
            % (base-identical route bookkeeping)
            strain_p_i = Contraction_k( ...
                obj, Lmt_p_i, obj.gama, []);
            obj.strain_p = strain_p_i;

            % Force with stiffness effects
            F_p_i = Force_p(obj, obj.uD_p, obj.strain_p);
            obj.F_p = F_p_i;

            % Moment arm with stiffness effects
            mA_p_i = Mom(obj, obj.L_p, obj.uD_p);
            obj.mA_p = mA_p_i;

            % Torque with stiffness effects; infeasible frames are NaN'd
            % per the pulley contract on top of the base strain rule.
            Torque_p_i = Tor(obj, obj.mA_p, obj.F_p, obj.strain_p);
            Torque_p_i(infeas_i, :) = NaN;
            obj.Torque_p = Torque_p_i;

            % --------------------------------------------------------
            % Reverse-pulley outputs (equilibrium-honest quantities)
            % --------------------------------------------------------
            G = obj.PulleyGain;
            obj.Ftendon = obj.NPulleyBPA * Fmag_i / G;
            % Slack tendon carries no tension: the solve clamps to s = 0
            % but keeps the zero-strain BPA force in Fmag, so without
            % this the taut tension would leak into F_ins and the mount
            % reaction at slack frames.
            obj.Ftendon(slack_i) = 0;

            % Tackle input travel: solved contraction minus the
            % installation offset and the bracket's share of the stroke,
            % so that G*PulleyTravel = DeltaL + delta_t holds exactly.
            obj.PulleyTravel = sCon_i - s0(obj.RefIndex) - Fmag_i .* cb_i;

            % Tendon-side unit direction at the insertion, in the tibia
            % frame. The force on the tibia pulls the insertion toward the
            % exit (base UD sign convention: proximal-side point minus
            % distal point). moving_via: recomputed per orientation from
            % the DEFORMED exit point, so the moving exit/wrap point is
            % honored and the moment arm can grow with joint angle.
            % bowden: the housing anchor sits on the tibia, so the line
            % direction is fixed in the tibia frame (held at theta0).
            u_t_i = zeros(N, 3);
            for ii = 1:N
                exitT1 = RowVecTrans(obj.TransformationMat(:,:,ii)\eye(4), ...
                    obj.L_p(E,:,ii));
                u_t_i(ii, :) = exitT1 - obj.L_p(C,:,ii);
            end
            u_t_i = normalize(u_t_i);
            if strcmp(obj.RoutingMode, 'bowden')
                u_t_i = repmat(u_t_i(obj.RefIndex, :), N, 1);
            end
            obj.u_t = u_t_i;

            % BPA-side total force direction expressed in the tibia frame
            % so the reaction sum is frame-consistent. DIRECTIONS rotate
            % but do not translate: apply only the rotation block of
            % T^-1 (a full RowVecTrans point transform would add the
            % joint offset and corrupt the direction).
            u_bpa_i = zeros(N, 3);
            for ii = 1:N
                Tinv = obj.TransformationMat(:,:,ii)\eye(4);
                u_bpa_i(ii, :) = (Tinv(1:3, 1:3) * Funit_i(ii, :).').';
            end
            u_bpa_i = normalize(u_bpa_i);
            obj.u_bpa = u_bpa_i;

            % Insertion force, moment arm, and torque about the joint
            F_ins_i = obj.Ftendon .* u_t_i;
            F_ins_i(slack_i, :) = 0;        % slack tendon transmits nothing
            F_ins_i(infeas_i, :) = NaN;
            obj.F_ins = F_ins_i;

            mA_ins_i = zeros(N, 3);
            for ii = 1:N
                pointB = obj.L_p(C,:,ii);
                mA_ins_i(ii, :) = pointB - u_t_i(ii, :)*dot(u_t_i(ii, :), pointB);
            end
            obj.mA_ins = mA_ins_i;

            Torque_ins_i = zeros(N, 3);
            for ii = 1:N
                Torque_ins_i(ii, :) = cross(mA_ins_i(ii, :), F_ins_i(ii, :));
            end
            Torque_ins_i(infeas_i, :) = NaN;
            obj.Torque_ins = Torque_ins_i;

            % Mount reaction at the pulley block: the nPulleyBPA parallel
            % BPA lines plus the tendon line (NaN'd where infeasible; at
            % slack the tendon term vanishes on its own).
            Reaction_i = (obj.NPulleyBPA * Fmag_i) .* u_bpa_i ...
                + obj.Ftendon .* u_t_i;
            Reaction_i(infeas_i, :) = NaN;
            obj.ReactionF = Reaction_i;
            obj.ReactionFmag = vecnorm(obj.ReactionF, 2, 2);

            obj.PulleyInfeasible = infeas_i;
            obj.PulleySlack = slack_i;
        end

    end % methods

end % classdef

%% =====================================================================
%% Helper functions (derived from the minimizer calculations)
%% =====================================================================

%% ------------- Pulley config field ----------------
function val = cfgField(cfg, name, default)
% Pull one optional field out of the pulley config struct.
    if isfield(cfg, name) && ~isempty(cfg.(name))
        val = cfg.(name);
    else
        val = default;
    end
end

%% -------------Force unit direction ---------------
function F_unit = computeForceVector(obj)
%Calculate the force unit direction from muscle origin (hip frame) to the next
%real point. This takes into account if there are any additional via
%points between muscle origin and muscle insertion. It also takes into
%account if a homogenous transform+ation matrix needs to be used to
%convert the second point into the first points frame.

L = obj.Location;      %Location (wrapping, attachment points)
C = obj.Cross;         %Cross point (moves from one frame to another)
T = obj.TransformationMat;       %Transformation matrix

% Step 1: Detect the first valid segment (non-repeated)
N = size(L, 3);      % Number of samples/frames
pt1 = squeeze(L(1,:,:))';  % Origin point, 3*N -> N*3
pt2 = NaN(N, 3);

for i = 1:N
    pt1(i,:) = L(1,:,i);  % muscle origin
    found = false;

    for k = 2:size(L,1)
        d = norm(L(k,:,i) - L(1,:,i));
        if d > 1e-6
            if k == C
                pt2(i,:) = RowVecTrans(T(:,:,i), L(k,:,i));
            else
                pt2(i,:) = L(k,:,i);
            end
            found = true;
            break;
        end
    end

    if ~found
        warning("Frame %d: No valid second point, using pt1=pt2", i);
        pt2(i,:) = pt1(i,:);  % fallback
    end
end
% Force Direction vector (hip frame)
F_vec = pt2 - pt1;
F_unit = normalize(F_vec);

end

%% -------------- Contraction of the PAM --------------------------
function contraction = Contraction_k(obj,Lmt,gema,X0)
rest = obj.RestingL;      %resting length
tendon = obj.TendonL;     %artificial tendon length
fitting = obj.FittingLength;   %fitting length

if isempty(Lmt)
    Lmt = obj.MuscleLength;
end

if isempty(gema)
    gema = 0;
end

if isempty(X0)
    X0 = 0;
end

Lm = Lmt-tendon-gema-2*fitting-X0;  %active BPA muscle length
contraction = (rest-Lm)/rest;    %contracted percent of original
end

%% ------------- Location  ------------------------
function [LOC, gema, sCon, Fmag, infeasible, slack, cbOut] = Lok(obj,X1,X2,kSpr,Funit,strain_predef,X0,DeltaL,s0,s0ref)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   Funit, force unit direction in the hip frame
%   strain_predef - N*1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
%   DeltaL, N*1 rigid tendon-span change from theta0 (pulley closure)
%   s0, N*1 rigid-geometry contraction; s0ref, its value at theta0
L = obj.Location;          %Location of wrapping and attachment points
rest = obj.RestingL;      %resting length
Fm = obj.Fmax;          %maximum isometric force (ONE BPA)
P = obj.Pressure;            %BPA pressure
D = obj.Diameter;         %BPA diameter
KMAX = (rest - obj.Kmax)/rest; %maximum contracted length (meters)
N = size(L,3);

% Compute Force
relstrain = strain_predef / KMAX;  %Relative strain
FF = festo4(D, relstrain, P) * Fm; %Force magnitude
FF (FF < 0) = 0;
F = FF.*Funit;  % N*3, already in hip frame

% Frame-count reference: the base class hard-codes frame 92; derive it
% from size(Location,3) (identical behavior at the base's 92-orientation
% sweeps, last frame otherwise).
iRef = min(92, N);
pA = L(1,:,iRef);                                  %Distance from hip origin to muscle insertion
switch obj.Diameter
    case 20
%       Pbr = [-0.8100  -20.222   31.66]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
        Pbr = [9.48  -36.21   30.86]/1000;       %from hip origin to bracket bolt pattern centroid
    case 10
        Pbr = [-19 22 27.6]/1000;       %from hip origin centroid of bracket cantilever
%         Pbr = [-21.33  -79   6.94]/1000;       %from centroid of bracket bolts.
    otherwise
        Pbr = [0 0 0];
end

phbrA = pA-Pbr;                                  %vector from bracket to point A (in the hip frame)
thetabrA = atan2(phbrA(2),phbrA(1));             %angle between pbrA and x axis
RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...     %Rotation matrix
       sin(thetabrA) cos(thetabrA) 0; ...
       0    0   1];
pbrhA = RhbrZ'*phbrA';       %Vector in the bracket frame
% Now calculate angle from x-axis to this vector
thetaY = atan2(pbrhA(3), pbrhA(1));  % z vs x (in bracket frame)

% Rotation matrix about y-axis (local frame adjustment)
Ry = [cos(thetaY)  0  sin(thetaY);
      0            1  0;
     -sin(thetaY) 0   cos(thetaY)];
Rhbr = RhbrZ*Ry';            %Rotate about y-axis in body frame
Thbr = RpToTrans(Rhbr, Pbr');    %Transformation matrix, represent bracket frame in hip frame

Fbrh = zeros(N,3);
pAnew = zeros(N,3);     %New point A, in the hip frame
for ii = 1:N                          %Repeat for each orientation
    Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),F(ii,:));            %Force vector in the hip frame represented in the bracket frame
end
if isinf(X1) && isinf(X2) && isinf(kSpr)
    [epsilon, delta, beta, gema, sCon, Fmag, cbOut] = ...
        deal(zeros(N,1));
    sCon = s0;      %rigid escape: no elastic absorption anywhere
    infeasible = false(N,1);
    slack = false(N,1);
else
    [epsilon, delta, beta, gema, sCon, Fmag, infeasible, slack, cbOut] = ...
        fortz(obj,Fbrh,X1,X2,kSpr,DeltaL,s0,s0ref);  %pulley equilibrium
end
deflection = [epsilon, delta, beta];    %bracket movement
pbrAnew = [norm(pbrhA),0,0]+deflection; %New point A, represented in the bracket frame
% pbrAnew = [norm(pbrhA(1:2)),0,pbrhA(3)]+deflection; %New point A, represented in the bracket frame
LOC = L;
for ii = 1:N                          %Repeat for each orientation
    pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
    LOC(1,:,ii) = pAnew(ii,:);      %Update location matrix
end

end

%% Pulley equilibrium: solve s from G*(s - s0ref - cb*F) = DeltaL + delta_t
%% with delta_t = nPulleyBPA*F/(G*kSpr) (tackle equilibrium F_t = F_BPA/G)
function [e_axial, e_bendY, e_bendZ, e_cable, sCon, Fmag, infeasible, slack, cbOut] = fortz(obj,Fbr,X1,X2,kSpr,DeltaL,s0,s0ref)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch (= delta_t, the stretch beyond the
%          blocked reference state; the ABSOLUTE tension stays
%          kSpr*e_cable = F_t)
% sCon, solved absolute BPA contraction; Fmag, its single-BPA force
% infeasible/slack, N*1 logical flags
% cbOut, per-frame bracket compliance along the SAME direction u_hat the
% equilibrium used (RowVecTrans is affine, so re-deriving the direction
% elsewhere would not match; PulleyTravel needs this exact cb)

N = size(Fbr,1);
% Initialize outputs
[e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N,1));
[sCon, Fmag] = deal(zeros(N,1));
cbOut = zeros(N,1);
infeasible = false(N,1);
slack = false(N,1);

D = obj.Diameter;         %BPA diameter
rest = obj.RestingL;      %resting length
mif = obj.Fmax;         %maximum force of ONE BPA
kmax = obj.Kmax;      %maximum contracted length
KMAX = (rest-kmax)/rest; %turn it into a percentage
P = obj.Pressure;            %pressure
G = obj.PulleyGain;          %tackle gain
nBPA = obj.NPulleyBPA;       %parallel BPAs feeding the tackle

% Normalize force vectors safely
norms = vecnorm(Fbr, 2, 2);
valid = norms > 1e-3 & all(~isnan(Fbr), 2);
u_hat_all = normalize(Fbr);

% Vectorized k_b computation
K = [X1, X2, X2];   % bracket stiffness array
C_bracket = diag([1/K(1), 1/K(2), 1/K(3)]);             % compliance matrix
cSpr = 1/kSpr;              % compliance of the ONE tendon line

% Force law at absolute contraction s for ONE BPA (festo4 returns >= 0,
% and 0 for relative strain above 1)
Ffun = @(s) festo4(D, (s/rest)/KMAX, P) * mif;

for i = 1:N
    if ~valid(i)
        sCon(i) = s0(i);    %base-identical escape: no deformation solved
        continue;
    end

    u_hat_i = u_hat_all(i, :);
    cb_i = (u_hat_i*C_bracket)*u_hat_i.';   %bracket compliance along u
    cbOut(i) = cb_i;

    % Scalar equilibrium (see class header): the tackle multiplies the
    % BPA-side travel it receives by G, so in absolute-contraction form
    %   G*(s - s0ref - cb*F(s)) = DeltaL + nBPA*F(s)/(G*kSpr).
    % At G = 1, nBPA = 1 this is algebraically the base class solve
    % r = s - s0(i) = (cb + cSpr)*F(s), because s0ref + DeltaL = s0(i)
    % for any route whose non-tendon rows are orientation-constant.
    psiFun = @(s) G*s - G*cb_i*Ffun(s) - nBPA*Ffun(s)*cSpr/G ...
        - DeltaL(i) - G*s0ref;

    sLo = 0;                %contract per the contract's clamp
    sHi = KMAX*rest;        %fully contracted
    psiLo = psiFun(sLo);
    psiHi = psiFun(sHi);

    if psiLo > 0
        % Even at zero contraction the span demands less takeup than the
        % tackle would feed: the tendon goes slack. Clamp to s = 0; the
        % mount still deflects under the zero-strain BPA force.
        sCon(i) = sLo;
        Fmag(i) = Ffun(sLo);
        e_bkt = C_bracket*(Fmag(i)*u_hat_i.');
        e_axial(i) = e_bkt(1);
        e_bendY(i) = e_bkt(2);
        e_bendZ(i) = e_bkt(3);
        e_cable(i) = 0;     %slack tendon carries no stretch
        slack(i) = true;
    elseif psiHi < 0
        % Required contraction exceeds KMAX*Rest: infeasible frame. Clamp
        % to full contraction (the force law gives ~0 there) and flag;
        % the torque is NaN'd by the caller.
        sCon(i) = sHi;
        Fmag(i) = Ffun(sHi);
        e_bkt = C_bracket*(Fmag(i)*u_hat_i.');
        e_axial(i) = e_bkt(1);
        e_bendY(i) = e_bkt(2);
        e_bendZ(i) = e_bkt(3);
        e_cable(i) = nBPA*Fmag(i)*cSpr/G;
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
            e_cable(i) = nBPA*Fmag(i)*cSpr/G;  %tendon-side stretch
            continue;
        else
            % Bracket displacement under the BPA-side force
            e_bkt = C_bracket * (Fmag(i) * u_hat_i.');

            e_axial(i) = e_bkt(1);
            e_bendY(i) = e_bkt(2);
            e_bendZ(i) = e_bkt(3);

            % Tendon stretch: the tackle divides the force by G, so the
            % stretch is nBPA*F/(G*kSpr) (base: F/kSpr).
            e_cable(i) = nBPA*Fmag(i)*cSpr/G;
        end
    end
end
end

%% ------------- Segment Lengths ------------------------
function SL = seg(obj, L)
C = obj.Cross;
T = obj.TransformationMat;
N = size(T, 3);
M = size(L, 1);
SL = zeros(N, M-1);

for ii = 1:N                    %Repeat for each orientation
    for i = 1:M-1               %Calculate all segments
        pointA = L(i,:,ii);
        pointB = L(i+1,:,ii);
        if i+1 == C
            pointB = RowVecTrans(T(:,:,ii), pointB);
        end
        SL(ii,i) = norm(pointA - pointB);
    end
end
end

%% ------------- Muscle Length ------------------------
%Function that calculates the musclutendon length
function Lmt = LMT(sL, X0)
% Compute muscle-tendon length from segment lengths and offset Xi0
Lmt = sum(sL, 2);  % Nx1, sum across segments

if ~isempty(X0)
    Lmt = Lmt - X0;
end
end

%% -------------- Force Unit Direction ----------------
%Calculate the unit direction of the muscle force about the joint.
function unitD = UD(obj, L_p)
T = obj.TransformationMat;
C = obj.Cross;
direction = zeros(size(T, 3), 3);
unitD = zeros(size(direction));

for i = 1:size(T, 3)
    pointA = L_p(C-1, :, i);
    pointB = L_p(C, :, i);
    direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
    unitD(i, :) = direction(i, :)/norm(direction(i, :));
end
end

%% -------------- Moment Arm --------------------------
%Calculate the moment arm about the joint
%For every ViaPoint, calculate the moment arm of the muscle about
%the joint it crosses over
function mA = Mom(obj, L_p, unitD_p)
T = obj.TransformationMat;
C = obj.Cross;
mA = zeros(size(T, 3), 3);

for i = 1:size(T, 3)
    pointB = L_p(C, :, i);
    mA(i, :) = pointB - unitD_p(i, :)*dot(unitD_p(i, :), pointB);
end
end

%% -------------- Force --------------------------
%Calculate the direction of the forced applied by the muscle
function F = Force_p(obj, unitD_p, strain)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%pres == measured pressure
%kmax == maximum contraction length
%Outputs:
%F == Force, N
rest = obj.RestingL;
kmax = obj.Kmax;
KMAX = (rest-kmax)/rest; %turn it into a percentage

rel = strain./KMAX;                    %relative strain

Fn = festo4(obj.Diameter,rel,obj.Pressure);

scalarForce = Fn.*obj.Fmax;
scalarForce(scalarForce < 0) = 0;

F = scalarForce.*unitD_p;

end

%% ---------------------- Torque --------------
%Calculate torque by multiplying the the force along the
%Useful information
% i -> Index for Crossing Points/Joints
% ii -> Index for every degree of motion
% iii -> Index for axes of interest to observe Torque about
function Mz = Tor(obj, mA_p, F_p, strain_p)
Mz = zeros(size(F_p));

switch obj.Diameter
    case 20
        ss = -.03;       %maximum allowable strain
    case 10
        ss = -.02;
    otherwise
        ss = -.02;
end

for i = 1:size(F_p, 1)
    if strain_p(i,:) < ss
        Mz(i,:) = NaN;
    else
        Mz(i, :) = cross(mA_p(i, :), F_p(i, :));
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

%% Subfunctions
function vhat = normalize(v)
N = size(v,1);
norms = vecnorm(v,2,2);
valid = norms > 1e-3 & all(~isnan(v), 2);
vhat = zeros(N, 3);
vhat(valid, :) = v(valid, :) ./ norms(valid);
end
