%% Optimize predicted torque for extensors.
function [f_all, bpa_all] = minimizeFlxPin(Xi0,Xi1,Xi2,idx_val)
% minimizeExt: calculates predicted torque and fit metrics for a given BPA index
%
% Inputs:
%   Xi0 - extra length correction (m)
%   Xi1 - bracket "axial" stiffness (N/m)
%   Xi2 - bracket "bending" stiffness (N/m)
%   whichIdx - index of the BPA to process (scalar)
%
% Outputs:
%   fitvec - [RMSE, FVU, MaxResidual] for the selected BPA
%   bpa    - updated BPA struct with prediction fields filled in

%% load
%kf = knee flexor, kf(1) = specific resting length, kf(2) = biomimetic;
%ke = knee extensor, same as above
%example: kf(1).Mz z-axis torque for pinned knee, flexor, 46cm length
%'exp' suffix means experimentally measured
%'_h' suffix means hybrid calculated
%'_p' suffix means prime, as in the new prediction values

load FlxPinBPASet.mat kf %This loads the following, which was ran and saved:

    % load KneeFlxPin_10mm_48cm.mat Bifemsh_Pam phiD
    % Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    % G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    % load Plot_KneeFlxPin10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(1) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % 
    % % 46cm length
    % load KneeFlxPin_10mm_46cm.mat Bifemsh_Pam phiD
    % Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    % G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    % load Plot_KneeFlxPin10mm_46cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % % 47cm length
    % load KneeFlxPin_10mm_47cm.mat Bifemsh_Pam phiD
    % Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    % G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    % load Plot_KneeFlxPin10mm_47cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(3) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
    % 
    % % 40cm length
    % load KneeFlxPin_10mm_40cm.mat Bifemsh_Pam phiD
    % Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    % G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    % load Plot_KneeFlxPin10mm_40cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    % A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    % kf(4) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
    %               'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
    %               'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
    %               'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
    %               'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
    %               'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
    %               'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
    %               'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
    % clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

%% Initialize output
nBPA = numel(ke);
% Default to all BPAs if none specified
if nargin < 4 || isempty(idx_val)
        idx_val = 1:nBPA;
end

% [ke.strain_f] = deal([]);
% [ke.Lm_f] = deal([]);

bpa_all = ke;  % initialize
f_all = NaN(nBPA, 3);


%% Evaluate each BPA
for i = idx_val
%     fprintf('Evaluating BPA #%d with [%.4f, %.2e, %.2e]\n', i, Xi0, Xi1, Xi2, Xi3);
    klass_i = kf(i);
    [bpa_all(i), f_all(i,:)] = evaluateBPA(klass_i, Xi0, Xi1, Xi2, Xi3);
    if any(isnan(bpa_all(i).strain_p))
        warning('NaNs in strain_p for BPA #%d', i);
    end
end

end


function [bpa_i, fitvec] = evaluateBPA(klass, Xi0, Xi1, Xi2)
%% Calculate locations and properties
bpa_i = klass;  % 
kspr = Spr(bpa_i); %Calculate spring rate (Infinite if no tendon is used)
strain_Xi0 = Contraction(bpa_i, [],Xi0); %Calculate contraction with constant length offset
[L_p, gemma] = Lok(bpa_i, Xi1, Xi2, kspr, strain_Xi0, Xi0);   %Bracket deformation changing geometry using stiffnesses and constant length offset
unitD_p = UD(bpa_i, L_p);   %New force direction
sL_p = seg(bpa_i, L_p);   %New segment lengths uses deformation but does not subtract length offset
Lmt_p = LMT(sL_p, Xi0);     %New musclulotendon length. Uses deformed geometry and constant length offset.
strain_p = Contraction(bpa_i, Lmt_p, []);  %*new contraction amount includes deformed geometry and constant length offset
F_p = Force(bpa_i, unitD_p, strain_p);  %new force vector
mA_p = Mom(bpa_i, L_p, unitD_p);   %new moment arm
% M_p = Tor(mA_p, F_p, bpa_i.Fm, strain_p);  %new torque
M_p = Tor(mA_p, F_p, strain_p);  %new torque

%% Package into output struct
bpa = bpa_i;
bpa.Lmt_p = Lmt_p;
bpa.mA_p = mA_p;
bpa.M_p = M_p;
bpa.F_p = F_p;
bpa.strain_p = strain_p;
bpa.L_p = L_p;
bpa.gama = gemma;

% GoF calculation
fitvec = SSE(bpa_i, M_p);

%% Nested functions, modified from MonoPamExplicit

%% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass, L_mt, X0)
rest   = klass.rest;
tendon = klass.ten;
fitting = klass.fitn;

if isempty(X0)
    X0 = 0;
end

if isempty(L_mt)
    L_mt = klass.Lmt;
end

contraction = (rest-(L_mt-tendon-2*fitting-X0))/rest;    %(minus Xi0 is also used in LMT function)
end

%% ------------- Location  ------------------------
function [LOC, gama] = Lok(klass,X1,X2,kSpr,strain,X0)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   Funit, force unit direction in the hip frame
%   strain – N×1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
L = klass.Loc;      %Location (wrapping, attachment points)
C = klass.CP;       %Cross point (moves from one frame to another)
%             T = klass.Tk;       %Transformation matrix
kmax = klass.Kmax;  %max contracted length
KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage

relstrain = strain/KMAX;            %relative strain
FF = festo4(klass.dBPA,relstrain,klass.P).*klass.Fm;        %Force magnitude
FF(FF<0) = 0;
%For muscle insertion
unitD = klass.unitD;            %unit direction of force vector, tibia frame
Fk = unitD.*FF;                  %Force vector, tibia frame
pB = L(C,:,(klass.Ak==0));                  %Distance from knee frame to muscle insertion
% Pbri = [-48.11, -107.81, 13.8]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever)
% Pbri = [-27.5, -107.81, -0.54]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever, but at tibial contact, no z offset)
Pbri = [-27.5, -125.91, -0.54]/1000;     %vector from knee ICR to upper bolt
pkbrB = pB-Pbri;                  %vector from bracket to point B, in the knee frame
thetabrB = atan2(pkbrB(2),pkbrB(1));   %angle between pbrB and x axis
RkbrZ = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
       sin(thetabrB) cos(thetabrB) 0; ...
       0    0   1];
pbrkB = RkbrZ'*pkbrB';       %Vector in the bracket frame
%             Now calculate angle from x-axis to this vector
thetaY = atan2(pbrkB(3), pbrkB(1));  % z vs x (in bracket frame)
% % Rotation matrix about y-axis (local frame adjustment)
Ry = [cos(thetaY) 0  sin(thetaY);
      0           1  0;
     -sin(thetaY) 0  cos(thetaY)];
Rkbr = RkbrZ*Ry';            %Rotate about y-axis in body frame
Tkbr = RpToTrans(RkbrZ, Pbri');    %Transformation matrix, flexor bracket frame in knee frame

LOC = L;            %new location matrix
N = size(L,3);
%             M = size(L,1);
Fbrk = zeros(N,3);       %Force vector represented in the tibial bracket frame  

parfor ii = 1:N                          %Repeat for each orientation
        Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),Fk(ii,:)); %Force vector in the tibia frame represented in the lower bracket frame
end

if isinf(X1) && isinf(X2)  && isinf(kSpr)
    [epsilon, delta, beta, gama] = deal(zeros(N,1));
else
    [epsilon, delta, beta, gama] = fortz(klass,Fbrk,X1,X2,kSpr,X0);  %strain from force divided by tensile stiffness
end

eB = [epsilon, delta, beta];
pbrBnew = [norm(pkbrB(1:2)), 0, pkbrB(3)] + eB; %new point B, in the bracket's frame
%             pbrBnew = [norm(pkbrB), 0, 0] + eB; %new point B, in the bracket's frame

pBnew = zeros(N,3);
for ii = 1:N                          %Repeat for each orientation
    pBnew(ii,:) = RowVecTrans(Tkbr, pbrBnew(ii,:));     %New point B, in the tibia frame
    LOC(2,:,ii) = pBnew(ii,:);
end
end

%% Force and length reduction due to deformation
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(klass,Fbr,X1,X2,kSpr,X0)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change
    N = size(Fbr,1);
    % Initialize outputs
    [e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N,1));
    
    if isinf(X1) && isinf(X2) && isinf(kSpr)
        return
    end
    
    if isempty(X0)
        X0=0;
    end
    D = klass.dBPA;         %BPA diameter
    rest = klass.rest;      %resting length
    tendon = klass.ten;     %tendon length
    fitn = klass.fitn;    %fitting length
    mL = klass.Lmt - X0 -tendon -2*fitn;       %Musculotendon length
    mif = klass.Fm;         %maximum force
    kmax = klass.Kmax;      %maximum contracted length
    KMAX = (rest-kmax)/rest; %turn it into a percentage
    P = klass.P;            %pressure
    
    % Normalize force vectors safely
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-4 & all(~isnan(Fbr), 2);
    u_hat_all = normalize(Fbr);
    
    % Vectorized k_b computation
    K_bracket = diag([X1, X2, X1]);       %bracket stiffness
    C_bracket = diag([1/X1, 1/X2, 1/X1]); %bracket compliance 
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    C_rep = repmat(C_bracket, [1, 1, N]);   % [3x3xN]
    c_b = pagemtimes(pagemtimes(u_hat, C_rep), permute(u_hat, [2, 1, 3]));
    c_b = reshape(c_b, [N, 1]);
    cSpr = 1/kSpr;
    c_eff = c_b+cSpr;        %effective compliance
    k_eff = 1 ./ c_eff;         % effective stiffness along u
    
    
    % Parallel root solve
   parfor i = 1:N
        if ~valid(i)
            continue;
        end
    
        % Per-instance constants
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
function SL = seg(klass, L)
C = klass.CP;
T = klass.Tk;
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
% Xi0 can be empty [] to skip correction (i.e., when already applied)
% N = size(sL_p, 1);
Lmt = sum(sL, 2);  % Nx1, sum across segments

if ~isempty(X0)
    Lmt = Lmt - X0;
end
end          

        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint (tibia frame).
function unitD = UD(klass, L)
            T = klass.Tk;
            C = klass.CP;
            direction = zeros(size(T, 3), 3);
%             unitD = zeros(size(direction));
            
            for i = 1:size(T, 3)
                pointA = L(C-1, :, i);
                pointB = L(C, :, i);
                direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
%                 unitD(i, :) = direction(i, :)/norm(direction(i, :));
            end
            unitD = normalize(direction);
end
        
        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
        %For every ViaPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
function mA = Mom(klass, L, unitD)
            T = klass.Tk;
            C = klass.CP;
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L(C, :, i);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
            end
end               


        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
function F = Force(klass, unitD, strain)
        %Inputs:
        %Lmt == muscle-tendon length, scalar
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %pres == measured pressure
        %kmax == maximum contraction length
        %Outputs:
        %F == Force, N           
           rest = klass.rest;
           kmax = klass.Kmax;  
           KMAX = (rest-kmax)/rest; %turn it into a percentage 
            
           rel = strain./KMAX;                    %relative strain        
           
           Fn = festo4(klass.dBPA,rel,klass.P);

           scalarForce = Fn.*klass.Fm;
           scalarForce(scalarForce <= 0) = 0;
%            scalarForce(scalarForce > maxF) = NaN;            
            
            F = scalarForce.*unitD;

end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
function Mz = Tor(mA, F, strain)
N = size(F, 1);
Mz = zeros(N, 3);

    for i = 1:N
        % if norm(F(i,:)) > maxF || strain(i,:) < -0.03
        if strain(i,:) < -0.03
            Mz(i,:) = NaN;
        else
            Mz(i,:) = cross(mA(i,:), F(i,:));
        end
    end
end   

%% tendon springrate
function springrate = Spr(klass)
    if klass.ten > 0
        mult = 2;           %Multiplier for number of cables used.
        Aeff = 1.51*10^-6;  %Effective area for 19-strand cable
        E = 193*10^9;       %Young's Modulus
        L = klass.ten;      %tendon length        
        springrate = mult*Aeff*E/L;
    else
        springrate = Inf;
    end
        
end

%% Subfunctions
function t = SSE(klass, M_p)
     [Ak_sorted, idx] = sort(klass.Ak);
     M_sorted = M_p(idx, 3);
     Mpredict2 = griddedInterpolant(Ak_sorted, M_sorted);
     M_opt = Mpredict2(klass.Aexp);
     [RMSE, fvu, maxResid] = Go_OfF(klass.Mexp,M_opt);
     t = [RMSE, fvu, maxResid];
end

function vhat = normalize(v)
    N = size(v,1);
    norms = vecnorm(v,2,2);
    valid = norms > 1e-4 & all(~isnan(v), 2);
    vhat = zeros(N, 3);
    vhat(valid, :) = v(valid, :) ./ norms(valid);
end


end   