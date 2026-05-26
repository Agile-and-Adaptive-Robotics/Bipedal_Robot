% Pam Data
% Author: Connor Morrow & Ben Bolen
% Date: 6/2026
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification
%
%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPamDataExplicit_balance < handle
    
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
        
        % --- Stiffness-aware fields (minimizeFlx-style) ---
        L_p                         %Deformed location matrix (updated attachment points)
        strain_p                    %Contraction with bracket deformation, tendon stretch, and constant length offset
        F_p                         %Force vector with stiffness effects
        mA_p                        %Moment arm with stiffness effects
        Torque_p                    %Torque with stiffness effects
        gama                        %Tendon stretch (cable elongation)
        kSpr                        %Tendon spring rate (effective)
    end
    
    %Dependent properties are those that are calculated by the explicit
    %properties. Matlab will not calculate these until it is queried in the
    %main script
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
        %Constructor Function. By calling 'MuscleData' and entering the
        %muscle information, we construct an object for that muscle.
        function PD = MonoPamDataExplicit_balance(name, location, cross, diameter, t, rest, kmax, tendon, fitn, pres, xi0, xi1, xi2, wraps)
            if nargin == 14
                PD.Name = name;
                PD.Location = location;
                PD.Cross = cross;
                PD.Diameter = diameter;
                PD.TransformationMat = t;
                PD.RestingL = rest;
                PD.Kmax = kmax;
                PD.TendonL = tendon;
                PD.FittingLength = fitn;
                PD.Pressure = pres;
                
                PD.Xi0 = xi0;
                PD.Xi1 = xi1;
                PD.Xi2 = xi2;
                PD.Wraps = wraps;
                
                % Automatically compute stiffness-aware geometry and torque
                PD = PD.updateStiffnessGeometry();
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
            avgSegL = zeros(size(L, 1, 1) - 1);
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
                for i = 1:size(L, 1, 1)-1                      %Repeat for all muscle segments
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
                maxF = 1500;
            elseif dia ==40
                maxF = 6000;
            else
                disp('Wrong size diameter BPA')
            end
        end

        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
        function F = get.Force(obj)
        %Inputs:
        %Lmt == muscle-tendon length, scalar
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %pres == measured pressure
        %kmax == maximum contraction length
        %Outputs:
        %F == Force, N           
            dia = obj.Diameter;
            unitD = obj.UnitDirection;
            contract = obj.Contraction;
            mL = obj.MuscleLength;
            rest = obj.RestingL;
            fitting = obj.FittingLength;
            pres = obj.Pressure;
            kmax = obj.Kmax;  
            KMAX = (rest-kmax)/rest; %turn it into a percentage 
            maxF = obj.Fmax;
            stiff = 0;               % no longer used for balance here
            tendon =  obj.TendonL;   %Length of artificial tendon and air fittings

           sF = zeros(size(obj.UnitDirection));
            
            %This function will solve for a new muscle length that
            %balances tendon force (fT) and the muscle force (fM) 
            function balanceF = muscleF(nmL)
                k = (rest-(mL-nmL-tendon-2.*fitting))./rest;  %strain 
                rel = k./KMAX;                             %relative strain  
                
               if dia == 10
                Fz = cell(1,4);
                [Fz{1}, Fz{2}, Fz{3}, Fz{4}] = normF10(rel, pres);
                Fn =Fz{3};
                fM = Fn.*maxF;
               elseif dia ~= 10
                fM = festo4(dia, rel, pres); 
               end

                fT = stiff*nmL;

                %Balance of the muscle forces. Solving to find when it
                %becomes equal to 0
                balanceF = fM - fT;
            end
            
          %Repeat the force calculation for every rotation of the joint
           for i = 1:size(unitD, 1)
                ML = obj.MuscleLength(i);              %Muscle-Tendon Length, which is the full calculated length between points in OpenSim

                %Set Function solver parameters
                options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);

                %Determine the normalized muscle length that solves the force
                %equations
                snmL = fsolve(@muscleF, ML, options);

                %With the solved muscle length, we can determine the scalar muscle
                %force by plugging it back into one of the force equations (fT)
                sF = stiff.*snmL;
                
           end

            for i = 1:size(unitD, 1)
                if sF(i) < 0
                    sF(i) = 0;
                end
                if sF(i) > maxF
                    sF(i) = NaN;
                end
            end
            
            SF = diag(sF);
            F = SF*unitD;

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
        % === Stiffness-aware pipeline (minimizeFlx-style) ===========
        % ============================================================
        
        % Xi0: constant length offset
        % Xi1, Xi2: bracket stiffness components
        % Wraps: number of cable wraps (affects tendon spring rate)
        function obj = updateStiffnessGeometry(obj)
            % Tendon spring rate
            obj.kSpr = Spr_local(obj, obj.Wraps);
            
            % Force unit vector in hip frame (origin to first non-duplicate point)
            Funit = computeForceVector_local(obj);
            
            % Contraction from constant length offset only
            strain_Xi0 = Contraction_local(obj, [], [], obj.Xi0);
            
            % Deformed geometry and tendon stretch
            [L_p_local, gama_local] = Lok_local(obj, obj.Xi1, obj.Xi2, obj.kSpr, Funit, strain_Xi0, obj.Xi0);
            obj.L_p = L_p_local;
            obj.gama = gama_local;
            
            % Unit direction with deformed geometry
            unitD_p = UD_local(obj, obj.L_p);
            
            % Segment lengths with deformed geometry
            sL_p = seg_local(obj, obj.L_p);
            
            % Musculotendon length with deformed geometry and Xi0
            Lmt_p_local = LMT_local(sL_p, obj.Xi0);
            
            % Contraction with bracket deformation, tendon stretch, and Xi0
            strain_p_local = Contraction_local(obj, Lmt_p_local, obj.gama, []);
            obj.strain_p = strain_p_local;
            
            % Force with stiffness effects
            F_p_local = Force_local(obj, unitD_p, obj.strain_p);
            obj.F_p = F_p_local;
            
            % Moment arm with stiffness effects
            mA_p_local = Mom_local(obj, obj.L_p, unitD_p);
            obj.mA_p = mA_p_local;
            
            % Torque with stiffness effects
            obj.Torque_p = Tor_local(obj, obj.mA_p, obj.F_p, obj.strain_p);
        end
        
    end % methods
    
end % classdef

%% =====================================================================
%% Local helper functions (mirroring minimizeFlx nested functions)
%% =====================================================================

%% -------------Force unit direction ---------------
function F_unit = computeForceVector_local(klass)
%Calculate the force unit direction from muscle origin (hip frame) to the next
%real point. This takes into account if there are any additional via
%points between muscle origin and muscle insertion. It also takes into
%account if a homogenous transform+ation matrix needs to be used to
%convert the second point into the first points frame.

L = klass.Location;      %Location (wrapping, attachment points)
C = klass.Cross;         %Cross point (moves from one frame to another)
T = klass.TransformationMat;       %Transformation matrix
    
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
F_unit = normalize_local(F_vec);

end

%% -------------- Contraction of the PAM --------------------------
function contraction = Contraction_local(klass,Lmt,gema,X0)
rest = klass.RestingL;      %resting length
tendon = klass.TendonL;     %artificial tendon length
fitting = klass.FittingLength;   %fitting length

if isempty(Lmt)
    Lmt = klass.MuscleLength;
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
function [LOC, gema] = Lok_local(klass,X1,X2,kSpr,Funit,strain_predef,X0)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   Funit, force unit direction in the hip frame
%   strain_predef – N×1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
L = klass.Location;          %Location of wrapping and attachment points
rest = klass.RestingL;      %resting length
Fm = klass.Fmax;          %maximum isometric force
P = klass.Pressure;            %BPA pressure
D = klass.Diameter;         %BPA diameter
KMAX = (rest - klass.Kmax)/rest; %maximum contracted length (meters)
N = size(L,3);

% Compute Force
relstrain = strain_predef / KMAX;  %Relative strain
FF = festo4(D, relstrain, P) * Fm; %Force magnitude
FF (FF < 0) = 0;
F = FF.*Funit;  % N×3, already in hip frame

pA = L(1,:,1);                                  %Distance from hip origin to muscle insertion
switch klass.Diameter
    case 20
%       Pbr = [-0.8100  -20.222   31.66]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
        Pbr = [9.48  -36.2   30.76]/1000;       %from hip origin to bracket bolt pattern centroid
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
Thbr = RpToTrans(RhbrZ, Pbr');    %Transformation matrix, represent bracket frame in hip frame              

Fbrh = zeros(N,3);
pAnew = zeros(N,3);     %New point A, in the hip frame
for ii = 1:N                          %Repeat for each orientation
    Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),F(ii,:));            %Force vector in the hip frame represented in the bracket frame
end
if isinf(X1) && isinf(X2) && isinf(kSpr)
    [epsilon, delta, beta, gema] = deal(zeros(N,1));
else
    [epsilon, delta, beta, gema] = fortz_local(klass,Fbrh,X1,X2,kSpr,X0);  %strain from force divided by tensile stiffness
end
deflection = [epsilon, delta, beta];    %bracket movement
% pbrAnew = [norm(pbrhA),0,0]+deflection; %New point A, represented in the bracket frame
pbrAnew = [norm(pbrhA(1:2)),0,pbrhA(3)]+deflection; %New point A, represented in the bracket frame
LOC = L;
for ii = 1:N                          %Repeat for each orientation 
    pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
    LOC(1,:,ii) = pAnew(ii,:);      %Update location matrix
end

end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz_local(klass,Fbr,X1,X2,kSpr,X0)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change
    
N = size(Fbr,1);
% Initialize outputs
[e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N,1));

D = klass.Diameter;         %BPA diameter
if isempty(X0)
    X0 = 0;
end
rest = klass.RestingL;      %resting length
tendon = klass.TendonL;      %tendon length
fitn = klass.FittingLength;    %fitting length
mL = klass.MuscleLength - X0 - tendon - 2*fitn;   %musculotendon length
mif = klass.Fmax;         %maximum force
kmax = klass.Kmax;      %maximum contracted length
KMAX = (rest-kmax)/rest; %turn it into a percentage
P = klass.Pressure;            %pressure

% Normalize force vectors safely
norms = vecnorm(Fbr, 2, 2);
valid = norms > 1e-3 & all(~isnan(Fbr), 2);
u_hat_all = normalize_local(Fbr);

% Vectorized k_b computation
K_bracket = diag([X1, X2, X1]);       %project bracket stiffness onto force direction
u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
K_rep = repmat(K_bracket, [1, 1, N]);   % [3x3xN]
k_b = pagemtimes(pagemtimes(u_hat, K_rep), permute(u_hat, [2, 1, 3]));
k_b = reshape(k_b, [N, 1]);
k_eff = 1 ./ (1 ./ k_b + 1 ./ kSpr);  % Nx1

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
            r = fzero(relfun, [0, .5]);
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
        e_cable(i) = F_mag/kSpr;
    end

end
end

%% ------------- Segment Lengths ------------------------
function SL = seg_local(klass, L)
C = klass.Cross;
T = klass.TransformationMat;
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
function Lmt = LMT_local(sL, X0)
% Compute muscle-tendon length from segment lengths and offset Xi0
Lmt = sum(sL, 2);  % Nx1, sum across segments

if ~isempty(X0)
    Lmt = Lmt - X0;
end
end          

%% -------------- Force Unit Direction ----------------
%Calculate the unit direction of the muscle force about the joint.
function unitD = UD_local(klass, L_p)
T = klass.TransformationMat;
C = klass.Cross;
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
function mA = Mom_local(klass, L_p, unitD_p)
T = klass.TransformationMat;
C = klass.Cross;
mA = zeros(size(T, 3), 3);

for i = 1:size(T, 3)
    pointB = L_p(C, :, i);
    mA(i, :) = pointB - unitD_p(i, :)*dot(unitD_p(i, :), pointB);
end
end        
        
%% -------------- Force --------------------------
%Calculate the direction of the forced applied by the muscle
function F = Force_local(klass, unitD_p, strain)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%pres == measured pressure
%kmax == maximum contraction length
%Outputs:
%F == Force, N           
rest = klass.RestingL;
kmax = klass.Kmax;  
KMAX = (rest-kmax)/rest; %turn it into a percentage 

rel = strain./KMAX;                    %relative strain        

Fn = festo4(klass.Diameter,rel,klass.Pressure);

scalarForce = Fn.*klass.Fmax;
scalarForce(scalarForce < 0) = 0;            

F = scalarForce.*unitD_p;

end
        
%% ---------------------- Torque --------------
%Calculate torque by multiplying the the force along the 
%Useful information
% i -> Index for Crossing Points/Joints
% ii -> Index for every degree of motion
% iii -> Index for axes of interest to observe Torque about
function Mz = Tor_local(klass, mA_p, F_p, strain_p)  
Mz = zeros(size(F_p));

switch klass.Diameter
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
function springrate = Spr_local(klass, wraps)
% wraps: number of cable wraps around the post (affects effective stiffness)
switch klass.Diameter
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
L = klass.TendonL;      %tendon length

springrate = mult*Aeff*E/L;        
end

%% Subfunctions
function vhat = normalize_local(v)
N = size(v,1);
norms = vecnorm(v,2,2);
valid = norms > 1e-3 & all(~isnan(v), 2);
vhat = zeros(N, 3);
vhat(valid, :) = v(valid, :) ./ norms(valid);
end
