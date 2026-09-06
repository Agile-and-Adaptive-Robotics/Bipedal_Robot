% Pam Data
% Author: Connor Morrow & Ben Bolen
% Date: 6/2026
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification
% The stiffness-aware pipeline was developed from the Xi0-Xi3 minimizer
% calculations. Standard Dependent properties are calculated when queried;
% updateStiffnessGeometry precomputes and stores the coupled deformation,
% tendon-stretch, force, moment-arm, and torque results during construction.
% The current Lok implementation deforms a femur/hip-frame bracket point.
%
%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPamDataExplicit_balanceX3 < handle
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location                    %BPA routing point locations
        Cross                       %Designates which Location row corresponds with a location where the muscle crosses into a new reference frame
        Diameter                    %Diameter of the BPA
        TransformationMat           %Contains a transformation matrix to change the 
        RestingL                    %Resting Length of the muscle
        Kmax                        %Length of BPA at maximum contraction
        FittingLength               %Length of each end cap (center of hole to bottom port)
        TendonL                     %Length of tendon, if any
        Pressure                    %Pressure of BPA
        AngleD                      %Angle vector in degrees
        
        % --- Stiffness parameters (always provided) ---
        Xi0                         %Constant length offset
        Xi1                         %Bracket axial stiffness
        Xi2                         %Bracket bending stiffness
        Xi3                         %Factor for loss of usable length
        BendMeasure                 % Nx1 geometric sum(R*alpha), m
        Wraps                       %Number of cable wraps (affects tendon stiffness)
        BPAcount                    %Number of BPAs in parallel

        % --- Stiffness-aware fields (minimizeExt-style) ---
        L_p                         %Deformed location matrix (updated attachment points)
        Lmt_p                       %Length of musculotendon including constant length offset Xi0
        delta_L                     %Loss of usable length due to bending
        strain_p                    %Contraction with bracket deformation, tendon stretch, Xi0, but not delta_L
        strain_f                    %Contraction with bracket deformation, tendon stretch, Xi0 and delta_L
        F_p                         %Force vector with stiffness effects
        mA_p                        %Moment arm with stiffness effects
        Torque_p                    %Torque with stiffness effects
        gama                        %Tendon stretch (cable elongation)
        kSpr                        %Tendon spring rate (effective)
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
        Fbal
        Torque
    end
    
   
    methods
        %% ------------- Muscle Data Constructor -----------------
        % Constructor function. The stored stiffness-aware fields are
        % populated once by updateStiffnessGeometry below.
        function PD = MonoPamDataExplicit_balanceX3(name, location, cross, diameter, t, rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, xi3, wraps, angleD, bpaCount, bendMeasure)
            if nargin == 17
                % Backwards compatibility with older callers.
                bendMeasure = [];

            elseif nargin ~= 18
                error('MonoPamDataExplicit_balanceX3:BadInputCount', ...
                    'Expected 17 or 18 inputs, got %d.', nargin);
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
                PD.Xi3 = xi3;                     % bend-loss scale factor
                PD.Wraps = wraps;                 % tendon wrap count
                PD.AngleD = angleD(:);            % joint-angle vector, degrees
                PD.BPAcount = bpaCount;            % equivalent parallel BPAs
                PD.BendMeasure = bendMeasure;      % R*angle by joint position, m
                % Automatically compute stiffness-aware geometry and torque.
                PD = PD.updateStiffnessGeometry();
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
%             contraction = zeros(length(mL), 1);
%             for i = 1:length(mL)
%                 contraction(i) = (rest-(mL(i,1)-tendon-2*fitting))/rest;
% %                   contraction(i) = ((rest+tendon+2*fitting)-mL(i,1))/(rest+tendon+2*fitting);
%             end
        end
        
        %% -------------- Length Check --------------------------
        function lengthCheck = get.LengthCheck(obj)
            contraction = obj.Contraction;
            maxContractPercent = 0.25;          %Contracting to 75% of length
            minContractPercent = -0.1;          %Elongating to 110% of length
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
                maxF = maxBPAforce(rest,620);
            elseif dia ==20
                maxF = maxBPAforce(rest,'20',620);
            elseif dia ==40
                maxF = 6000;
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
            
%             for i = 1:size(mA, 1)
%                 tor(i, :) = cross(mA(i, :), F(i, :));
%             end
        end    
        
        % ============================================================
        % === Stiffness-aware pipeline (minimizeExt-style) ===========
        % ============================================================
        
        % Xi0: constant length offset
        % Xi1, Xi2: bracket stiffness components
        % Xi3: factor for loss of usable length
        % Wraps: number of cable wraps (affects tendon spring rate)
        function obj = updateStiffnessGeometry(obj)

        % Tendon spring rate.
        % If two BPAs use two parallel tendon paths, multiply this by BPAcount.
        % If they share one tendon path, do not multiply it.
        obj.kSpr = obj.BPAcount * Spr(obj, obj.Wraps);
    
        % Force unit vector from origin
        Funit = computeForceVector(obj);
    
        % First pass:
        % Include Xi0 and Xi3 to estimate force for deformation calculation.
        [strain_Xi3, delta_L_i] = Contraction_k( ...
            obj, [], obj.Xi0, [], obj.Xi3);
    
        obj.delta_L = delta_L_i;
    
        % Bracket deformation and tendon stretch.
        % This matches minimizeExt:
        %   Lok(..., strain_Xi3, delta_L + Xi0)
        [L_p_i, gama_i] = Lok( ...
            obj, ...
            obj.Xi1, ...
            obj.Xi2, ...
            obj.kSpr, ...
            Funit, ...
            strain_Xi3, ...
            delta_L_i + obj.Xi0);
    
        obj.L_p = L_p_i;
        obj.gama = gama_i;
    
        % Segment lengths after deformation.
        sL_p = seg(obj, obj.L_p);
    
        % sL_p is the deformed path lengths without Xi0.
        % Therefore Lmt_p = sL_p - Xi0.
        Lmt_p_i = LMT( ...
            sL_p, obj.Xi0);
        obj.Lmt_p = Lmt_p_i;
    
        % Force strain:
        % includes Xi3 and is used to calculate BPA force.
        [strain_f_i, ~] = Contraction_k( ...
            obj, ...
            Lmt_p_i, ...
            [], ...
            obj.gama, ...
            obj.Xi3);
    
        obj.strain_f = strain_f_i;
    
        % Unit direction after deformation.
        unitD_p = UD(obj, obj.L_p);
    
        % Force_p(redicted) uses strain_f(ictive) to get accurate torque
        F_p_i = Force_p(obj, unitD_p, obj.strain_f);
        obj.F_p = F_p_i;
    
        % Moment arm after deformation.
        mA_p_i = Mom(obj, obj.L_p, unitD_p);
        obj.mA_p = mA_p_i;
    
        % Measured-comparison strain:
        % excludes Xi3. Use this to check whether the real BPA is stretching.
        [strain_p_i, ~] = Contraction_k( ...
            obj, ...
            Lmt_p_i, ...
            [], ...
            obj.gama, ...
            []);
    
        obj.strain_p = strain_p_i;
    
        % Torque uses F_p from strain_f, but Tor checks strain_p for stretching.
        obj.Torque_p = Tor(obj, obj.mA_p, obj.F_p, obj.strain_p);

end
        
    end % methods
    
end % classdef

%% =====================================================================
%% Helper functions (derived from the minimizer calculations)
%% =====================================================================

%% -------------Force unit direction ---------------
function F_unit = computeForceVector(obj)
%Calculate the force unit direction from muscle origin (hip frame) to the next
%real point. This takes into account if there are any additional via
%points between muscle origin and muscle insertion. It also takes into
%account if a homogenous transform+ation matrix needs to be used to
%convert the second point into the first points frame.

L = obj.Location;
C = obj.Cross;         %Cross point (moves from one frame to another)
T = obj.TransformationMat;       %Transformation matrix
    
% Step 1: Detect the first valid segment (non-repeated)
N = size(L, 3);      % Number of samples/frames
pt1 = squeeze(L(1,:,:))';  % Origin point, 3×N → N×3
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
function [contraction, delta_L] = Contraction_k(obj, Lmt, X0, gama, X3)

rest   = obj.RestingL;
tendon = obj.TendonL;
fitn   = obj.FittingLength;
theta_k = obj.AngleD(:);     % degrees
N = numel(theta_k);

KMAX = (rest - obj.Kmax)/rest;

if isempty(gama)
    gama = zeros(N,1);
end

if isempty(X0)
    X0 = 0;
end

if isempty(Lmt)
    Lmt = obj.MuscleLength;
end

delta_L = zeros(N,1);

if ~isempty(X3)

    Lm0 = Lmt - tendon - 2*fitn - X0 - gama;

    strain0 = (rest - Lm0) / rest;

    relstrain0 = strain0 / KMAX;

    comp = 1 - relstrain0;
    comp = max(0, comp);


    if ~isempty(obj.BendMeasure)

        bendMeasure = obj.BendMeasure(:);

        if numel(bendMeasure) ~= N
            error('MonoPamDataExplicit_balanceX3:BendMeasureSize', ...
                'BendMeasure must contain one value per knee position.')
        end

        % Geometric Xi3 model:
        %
        %   delta_L = Xi3 * sum(R*alpha) * comp^2
        %
        % Xi3 itself is unchanged.

        delta_L = ...
            X3 .* bendMeasure .* comp.^2;

    else

        % ---------------------------------------------------------
        % Backwards-compatible original pinned-knee formulation.
        % ---------------------------------------------------------

        ang = -9.19;

        angleRad = deg2rad((ang - theta_k)*80/(ang + 120));

        idx = angleRad > 0;

        R1 = 0.022;
        R2 = 0.176;

        delta_L1 = X3*R1*deg2rad(28).*comp.^2;

        delta_L2 = zeros(N,1);

        delta_L2(idx) = X3*R2 .*angleRad(idx).*comp(idx).^2;

        delta_L = delta_L1 + delta_L2;

    end

end

Lm_adj = Lmt - tendon - 2*fitn - X0 - gama - delta_L;
contraction = (rest - Lm_adj) / rest;

end

%% ------------- Location  ------------------------
function [LOC, gema] = Lok(obj,X1,X2,kSpr,Funit,strain_predef,deltaL)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   Funit, force unit direction in the hip frame
%   strain_predef – N×1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
L = obj.Location;          %Location of wrapping and attachment points
rest = obj.RestingL;      %resting length
Fm = obj.Fmax;          %maximum isometric force
P = obj.Pressure;            %BPA pressure
D = obj.Diameter;         %BPA diameter
KMAX = (rest - obj.Kmax)/rest; %maximum contracted length (meters)
N = size(L,3);

% Compute Force
relstrain = strain_predef / KMAX;  %Relative strain
FF = festo4(D, relstrain, P) .* Fm; %Force magnitude, single BPA
FF(FF < 0) = 0;

FF = obj.BPAcount .* FF;        %Now make it total force
F = FF .* Funit;                % Force vector N×3, already in hip frame

%Bracket transform
pA = L(1,:,92);               %Distance from hip origin to muscle insertion
switch obj.Diameter
    case 20
%       Pbr = [-0.8100  -20.222   31.66]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
        Pbr = [9.48  -33.38   30.86]/1000;       %from hip origin to bracket bolt pattern centroid
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
    [epsilon, delta, beta, gema] = deal(zeros(N,1));
else
    [epsilon, delta, beta, gema] = fortz(obj,Fbrh,X1,X2,kSpr,deltaL);  %shared length change from force balance
end
deflection = [epsilon, delta, beta];    %bracket movement
pbrAnew = [norm(pbrhA),0,0]+deflection; %New point A, represented in the bracket frame
% pbrAnew = [norm(pbrhA(1:2)),0,pbrhA(3)]+deflection; %New point A, represented in the bracket frame
LOC = L;

for ii = 1:N
    pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:));
    LOC(1,:,ii) = pAnew(ii,:);

    % Rows eliminated from the beginning of the femur-side route repeat
    % the original p1. Move those repeated rows with the deformed p1 so
    % they remain zero-length bookkeeping segments.
    for jj = 2:obj.Cross-1
        if norm(L(jj,:,ii) - L(1,:,ii)) < 1e-10
            LOC(jj,:,ii) = pAnew(ii,:);
        else
            break
        end
    end
end

end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(obj,Fbr,X1,X2,kSpr,deltaL)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change
    
N = size(Fbr,1);
% Initialize outputs
[e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N,1));

D = obj.Diameter;         %BPA diameter
if isempty(deltaL)
    deltaL = 0;
end
rest = obj.RestingL;      %resting length
tendon = obj.TendonL;      %tendon length
fitn = obj.FittingLength;    %fitting length
mL = obj.MuscleLength - deltaL - tendon - 2*fitn;   %musculotendon length
mif = obj.Fmax;         %maximum force
kmax = obj.Kmax;      %maximum contracted length
KMAX = (rest-kmax)/rest; %turn it into a percentage
P = obj.Pressure;            %pressure

% Normalize force vectors safely
norms = vecnorm(Fbr, 2, 2);
valid = norms > 1e-3 & all(~isnan(Fbr), 2);
u_hat_all = normalize(Fbr);

% Vectorized k_b computation
K = [X1, X2, X2]; %bracket stiffness array
K_bracket = diag(K);       %bracket stiffness matrix
C_bracket = diag([1/K(1), 1/K(2), 1/K(3)]);       %bracket compliance matrix
u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN] %reshape u_hat vector for all knee angles
C_rep = repmat(C_bracket, [1, 1, N]);   % [3x3xN] %repeat compliance bracket
c_b = pagemtimes(pagemtimes(u_hat, C_rep), permute(u_hat, [2, 1, 3])); %project bracket compliance onto force direction
c_b = reshape(c_b, [N, 1]);             %reshape bracket compliance                
cSpr = 1/kSpr;  %tendon compliance
c_eff = c_b+cSpr; %effective compliance
k_eff = 1 ./ c_eff;  % effective stiffness

parfor i = 1:N
    if ~valid(i)
        continue;
    end
    
    keff = k_eff(i);
    unit_vec = u_hat_all(i, :);
    Lm = mL(i);
    
    contraction0    = ( rest - Lm ) / rest;
    relstrain0      = contraction0 / KMAX;  %relative strain
    if relstrain0 >= 1
        r = 0;
    else
        relfun = @(r) ...
            festo4( D, ...
                (rest - (Lm -  r)) / rest / KMAX, ...
                P  ...
            ) * mif - keff * r;

        try
            r = fzero(relfun, [0, Lm-kmax]);
        catch
            r = 0;
        end
        r = max(r,0); %guard against r being slightly negative.
    end

    if r == 0
        continue;
    elseif isinf(X1) && isinf(X2)
        % Rigid body: no bracket deformation
        e_axial(i) = 0;
        e_bendY(i) = 0;
        e_bendZ(i) = 0;
        e_cable(i) = r;  % All elongation goes to cable
        continue;
    else
        % Final force magnitude
        contraction = (rest - (Lm - r)) / rest;
        relstrain = contraction / KMAX;
        F_mag = festo4(D, relstrain, P) * mif;

        % Bracket displacement
        e_bkt = K_bracket \ (F_mag * unit_vec');

        e_axial(i) = e_bkt(1);
        e_bendY(i) = e_bkt(2);
        e_bendZ(i) = e_bkt(3);
        % Cable elongation
        if tendon > 0
            r_bracket = unit_vec * e_bkt;
            r_cable = r-r_bracket;
            e_cable(i) = F_mag/kSpr;

            if abs(r_cable - e_cable(i)) > 1e-6
                warning('fortz:LengthBalanceMismatch', ...
                    'Frame %d: e_cable = %.9g, r_cable = %.9g, diff = %.9g', ...
                    i, e_cable(i), r_cable, r_cable - e_cable(i));
            end
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
%Calculate the moment arm about a joint
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

scalarForceSingle = Fn.*obj.Fmax;
scalarForceSingle(scalarForceSingle < 0) = 0;            

% Total force from parallel BPAs.
scalarForce = obj.BPAcount .* scalarForceSingle;

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

if nargin >= 2 && ~isempty(wraps)
    mult = mult * wraps;
end

Aeff = 1.51*10^-6;%Effective area for 19-strand cable
E = 193*10^9;       %Young's Modulus
L = obj.TendonL;      %tendon length

springrate = mult*Aeff*E/L;        
end

%% Subfunctions
function vhat = normalize(v)
N = size(v,1);
norms = vecnorm(v,2,2);
valid = norms > 1e-3 & all(~isnan(v), 2);
vhat = zeros(N, 3);
vhat(valid, :) = v(valid, :) ./ norms(valid);
end
