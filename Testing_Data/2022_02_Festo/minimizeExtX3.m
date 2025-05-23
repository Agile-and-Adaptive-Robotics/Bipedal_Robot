%% Optimize predicted torque for extensors.
function [f_all, bpa_all] = minimizeExtX3(Xi0, Xi1, Xi2, Xi3, idx_val)
% minimizeExt: calculates predicted torque and fit metrics for a given BPA index
%
% Inputs:
%   Xi0 - extra length correction (m)
%   Xi1 - axial bracket stiffness (N/m)
%   Xi2 - bending bracket stiffness (N/m)
%   Xi3 - curvature factor
%   whichIdx - index of the BPA to process (scalar)
%
% Outputs:
%   fitvec - [RMSE, FVU, MaxResidual] for the selected BPA
%   bpa    - updated BPA struct with prediction fields filled in

%% load
%kf = knee flexor, kf(1) = pinned joint, kf(2) = biomimetic;
%ke = knee extensor, same as above
%ke.L := lengths = [42 42 46 48] cm
%example: ke(1).L(3).Mz z-axis torque for pinned knee, flexor, 46cm length
%'exp' suffix means experimentally measured
%'_h' suffix means hybrid calculation
%'_p' suffix means prime, as in the new prediction values

load ExtPinBPASet.mat ke %This loads the following, which was ran and saved:

% % 42cm length, no tendon
%     load KneeExtPin_10mm_all.mat Vas_Pam_42cm phiD
%     Ma = Vas_Pam_42cm.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeExtPin10mm_42cm.mat Angle0 Torque0 InflatedLength0 ICRtoMuscle0 TorqueHand0
%     A = sortrows([Angle0, Torque0, InflatedLength0, ICRtoMuscle0, TorqueHand0]);
%     ke(1) = struct('Ak',phiD,'Loc',Vas_Pam_42cm.Location,'CP',Vas_Pam_42cm.Cross,'dBPA',Vas_Pam_42cm.Diameter, ...
%                   'Tk',Vas_Pam_42cm.TransformationMat,'rest',Vas_Pam_42cm.RestingL,'Kmax',Vas_Pam_42cm.Kmax,...
%                   'fitn',Vas_Pam_42cm.FittingLength,'ten',Vas_Pam_42cm.TendonL,'P',Vas_Pam_42cm.Pressure, ...
%                   'Lmt',Vas_Pam_42cm.MuscleLength,'strain',Vas_Pam_42cm.Contraction, 'unitD',Vas_Pam_42cm.UnitDirection, ...
%                   'mA',G,'Fm',Vas_Pam_42cm.Fmax,'F',Vas_Pam_42cm.Force, 'seg',Vas_Pam_42cm.SegmentLengths, ...
%                   'M',Vas_Pam_42cm.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
%                   'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Vas_Pam_42cm phiD Ma G Angle0 Torque0 InflatedLength0 ICRtoMuscle0 TorqueHand0 A
% 
% % 42cm length, tendon
%     load KneeExtPin_10mm_all.mat Vas_Pam_42cm_tendon phiD
%     Ma = Vas_Pam_42cm_tendon.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeExtPin10mm_42cm.mat Angle1 Torque1 InflatedLength1 ICRtoMuscle1 TorqueHand1
%     A = sortrows([Angle1, Torque1, InflatedLength1, ICRtoMuscle1, TorqueHand1]);
%     ke(2) = struct('Ak',phiD,'Loc',Vas_Pam_42cm_tendon.Location,'CP',Vas_Pam_42cm_tendon.Cross,'dBPA',Vas_Pam_42cm_tendon.Diameter, ...
%                   'Tk',Vas_Pam_42cm_tendon.TransformationMat,'rest',Vas_Pam_42cm_tendon.RestingL,'Kmax',Vas_Pam_42cm_tendon.Kmax,...
%                   'fitn',Vas_Pam_42cm_tendon.FittingLength,'ten',Vas_Pam_42cm_tendon.TendonL,'P',Vas_Pam_42cm_tendon.Pressure, ...
%                   'Lmt',Vas_Pam_42cm_tendon.MuscleLength,'strain',Vas_Pam_42cm_tendon.Contraction, 'unitD',Vas_Pam_42cm_tendon.UnitDirection, ...
%                   'mA',G,'Fm',Vas_Pam_42cm_tendon.Fmax,'F',Vas_Pam_42cm_tendon.Force, 'seg',Vas_Pam_42cm_tendon.SegmentLengths, ...
%                   'M',Vas_Pam_42cm_tendon.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
%                   'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Vas_Pam_42cm_tendon phiD Ma G Angle1 Torque1 InflatedLength1 ICRtoMuscle1 TorqueHand1 A
%     
% 
% % 46cm length
%     load KneeExtPin_10mm_all.mat Vas_Pam_46cm phiD
%     Ma = Vas_Pam_46cm.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeExtPin10mm_46cm.mat AngleX Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([AngleX, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
%     ke(3) = struct('Ak',phiD,'Loc',Vas_Pam_46cm.Location,'CP',Vas_Pam_46cm.Cross,'dBPA',Vas_Pam_46cm.Diameter, ...
%                   'Tk',Vas_Pam_46cm.TransformationMat,'rest',Vas_Pam_46cm.RestingL,'Kmax',Vas_Pam_46cm.Kmax,...
%                   'fitn',Vas_Pam_46cm.FittingLength,'ten',Vas_Pam_46cm.TendonL,'P',Vas_Pam_46cm.Pressure, ...
%                   'Lmt',Vas_Pam_46cm.MuscleLength,'strain',Vas_Pam_46cm.Contraction, 'unitD',Vas_Pam_46cm.UnitDirection, ...
%                   'mA',G,'Fm',Vas_Pam_46cm.Fmax,'F',Vas_Pam_46cm.Force, 'seg',Vas_Pam_46cm.SegmentLengths, ...
%                   'M',Vas_Pam_46cm.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
%                   'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Vas_Pam_46cm phiD Ma G AngleX Torque InflatedLength ICRtoMuscle TorqueHand A
% 
%     load KneeExtPin_10mm_all.mat Vas_Pam_48cm phiD
%     Ma = Vas_Pam_48cm.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeExtPin10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
%     ke(4) = struct('Ak',phiD,'Loc',Vas_Pam_48cm.Location,'CP',Vas_Pam_48cm.Cross,'dBPA',Vas_Pam_48cm.Diameter, ...
%                   'Tk',Vas_Pam_48cm.TransformationMat,'rest',Vas_Pam_48cm.RestingL,'Kmax',Vas_Pam_48cm.Kmax,...
%                   'fitn',Vas_Pam_48cm.FittingLength,'ten',Vas_Pam_48cm.TendonL,'P',Vas_Pam_48cm.Pressure, ...
%                   'Lmt',Vas_Pam_48cm.MuscleLength,'strain',Vas_Pam_48cm.Contraction, 'unitD',Vas_Pam_48cm.UnitDirection, ...
%                   'mA',G,'Fm',Vas_Pam_48cm.Fmax,'F',Vas_Pam_48cm.Force, 'seg',Vas_Pam_48cm.SegmentLengths, ...
%                   'M',Vas_Pam_48cm.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
%                   'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Vas_Pam_48cm phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A


%% Initialize output
nBPA = numel(ke);
bpa_all = ke;  % initialize
f_all = NaN(nBPA, 3);
idx_opt = setdiff(1:nBPA, idx_val);

%% Evaluate each BPA
for i = 1:nBPA
%     fprintf('Evaluating BPA #%d with [%.4f, %.2e, %.2e]\n', i, Xi0, Xi1, Xi2, Xi3);
    klass_i = ke(i);
    [bpa_all(i), f_all(i,:)] = evaluateBPA(klass_i, Xi0, Xi1, Xi2, Xi3);
    if any(isnan(bpa_all(i).strain_p))
        warning('NaNs in strain_p for BPA #%d', i);
    end
end

end


function [bpa_i, fitvec] = evaluateBPA(klass, Xi0, Xi1, Xi2, Xi3)
%% Calculate locations and properties

bpa_i = klass;
kspr = Spr(bpa_i);          %tendon spring rate
Funit = computeForceVector(bpa_i);  %Force unit direction in the hip frame 
strain_Xi3 = Contraction(bpa_i, [], [], [], Xi3);
[L_p, gemma] = Lok(bpa_i, Xi1, Xi2,kspr, Funit, strain_Xi3, []); %Bracket deformation and new geometry
sL_p = seg(bpa_i, L_p); %segment lengths
Lmt_p = LMT(sL_p, Xi0); 
strain_p = Contraction(bpa_i, Lmt_p, [], gemma, []);
unitD_p = UD(bpa_i, L_p);           %unit direction, predicted with updated path
F_p = Force(bpa_i, unitD_p, strain_p); %Muscle force, new prediction
mA_p = Mom(bpa_i, L_p, unitD_p); %Moment arm vector
M_p = Tor(mA_p, F_p, bpa_i.Fm, strain_p); %New torque prediction


%% Package into output struct
bpa_i.L_p = L_p;
bpa_i.Lmt_p = Lmt_p;
bpa_i.mA_p = mA_p;
bpa_i.F_p = F_p;
bpa_i.M_p = M_p;
bpa_i.strain_p   = strain_p;
bpa_i.gama   = gemma;

% SSE calculation
fitvec = SSE(bpa_i, M_p);

%% Nested functions, modified from MonoPamExplicit
%% -------------Force unit direction ---------------
function F_unit = computeForceVector(klass)
%Calculate the force unit direction from muscle origin (hip frame) to the next
%real point. This takes into account if there are any additional via
%points between muscle origin and muscle insertion. It also takes into
%account if a homogenous transformation matrix needs to be used to
%convert the second point into the first points frame.

L = klass.Loc;      %Location (wrapping, attachment points)
C = klass.CP;       %Cross point (moves from one frame to another)
T = klass.Tk;       %Transformation matrix
    
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

%     % Debug output
%     fprintf("Frame %3d | pt1: [% .4f, % .4f, % .4f]  pt2: [% .4f, % .4f, % .4f]  norm: %.6f\n", ...
%         i, pt1(i,:), pt2(i,:), norm(pt2(i,:) - pt1(i,:)));
end
% Force Direction vector (hip frame)
F_vec = pt2 - pt1;
F_unit = normalize(F_vec);

end

        %% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass, Lmt, X0, gama, X3)
rest   = klass.rest;
tendon = klass.ten;
fitn = klass.fitn;
theta_k = klass.Ak(:);  % degrees
N = length(theta_k);

if isempty(X0)
    X0 = 0;
end

if isempty(Lmt)
    Lmt = LMT(klass.seg,X0);
end

% --- Delta_L from curvature ---
delta_L = zeros(N,1);
if ~isempty(X3)
    ang2 = -23;
    angleRad = deg2rad(theta_k - ang2);
    idx = angleRad < 0;
    KMAX = (rest - klass.Kmax)/rest;
    Lm = (Lmt - tendon -2*fitn); %length of the BPA
    strain = (rest - Lm)/rest;
    relstrain = strain/KMAX;
%     comp = 1-relstrain;  %additive complement to relative strain
    R = 0.03;           %minimum radius
%     delta_L(idx) = X3*(R)*angleRad(idx).*comp(idx);
    delta_L(idx) = X3*(R)*abs(angleRad(idx)).*Lm(idx)/rest;
end

% --- Gama (deformation) ---
if isempty(gama)
    gama = zeros(N,1);
end

adjusted_Lm = Lmt - (tendon + gama + delta_L) - 2 * fitn;
contraction = (rest - adjusted_Lm) / rest;
end

%% Curvature in bending BPA estimate
function kappa = Curvature(angleRad)
    % Parameters
    R_min = 0.03;               % [m] min radius of curvature (30 mm)
    kappa_max = 1 ./ R_min;      % Max curvature
    angleMid = deg2rad(-60);    % Midpoint of transition (steepest change)
    sharpness = 1;             % Larger value => steeper transition

    % Logistic-style transition: 0 → kappa_max
    kappa = kappa_max ./ (1 + exp(sharpness .* (angleRad - angleMid)));

    % Optional clamp (shouldn't be necessary, but safe)
    kappa = min(kappa, kappa_max);
end

%% ------------- Location  ------------------------
function [LOC, gama] = Lok(klass, X1, X2, kSpr, Funit, strain_predef, X0)
% Inputs:
%   strain_predef – N×1 strain vector (e.g., from Xi0 + Xi3 curvature-only effect)

L = klass.Loc;
C = klass.CP;
% T = klass.Tk;
rest = klass.rest;
Fm = klass.Fm;
P = klass.P;
D = klass.dBPA;
KMAX = (rest - klass.Kmax)/rest;
N = size(L,3);
M = size(L,1);

% Compute Force
relstrain = strain_predef ./ KMAX;
FF = festo4(D, relstrain, P) .* Fm;
FF(relstrain > 1) = 0;
Fh = Funit .* FF;  % N×3, already in hip frame

%Bracket transform
pA = L(1,:,1);
Pbr = [-6.26, -29.69, 75.06]/1000;                          %from hip origin to lower bolt hole on superior anterior bracket of the Bifemsh_Pam
phbrA = pA-Pbr;                                  %vector from bracket to point A (in the hip frame)
thetabrA = atan2(phbrA(2),phbrA(1));            %angle between pbrA and x axis
RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...     %Rotation matrix
       sin(thetabrA) cos(thetabrA) 0; ...
       0    0   1];
pbrhA = RhbrZ'*phbrA';       %Vector in the bracket frame
%             pbrhA = RhbrZ'*phbrA';       %Vector in the bracket frame
%             % Now calculate angle from x-axis to this vector
%             thetaY = atan2(pbrhA(3), pbrhA(1));  % z vs x (in bracket frame)
% 
%             % Rotation matrix about y-axis (local frame adjustment)
%             Ry = [cos(thetaY)  0  sin(thetaY);
%                   0                1  0;
%                  -sin(thetaY) 0  cos(thetaY)];
%             Rhbr = RhbrZ*Ry;            %Rotate about y-axis in body frame
            Thbr = RpToTrans(RhbrZ, Pbr');    %Transformation matrix, represent bracket frame in hip frame  
            
% %more complicated way to calculate vector and rotation matrix so that your
% %new x axis points to muscle origin.
% p0 = [-6.26, -29.69, 75.06]/1000;
% y0 = [-47.48, -35.63, 75.06]/1000;
% z0 = [-8.39, -14.85, 75.06]/1000;
% 
% z_hat = normalize(p0 - z0);
% y_hat_prelim = normalize(p0 - y0);
% v = pA - p0;
% lambda = -dot(v, z_hat);
% Pbr_A = p0 + lambda * z_hat;
% 
% x_hat = normalize(cross(y_hat_prelim, z_hat));
% y_hat = cross(z_hat, x_hat);
% 
% Rhbr = [x_hat(:), y_hat(:), z_hat(:)];
% Thbr = RpToTrans(Rhbr, Pbr_A');
% pbrA = (pA - Pbr_A) * Rhbr;

% Transform force into bracket frame
Fbrh = zeros(N,3);
for ii = 1:N
    Fbrh(ii,:) = RowVecTrans(Thbr\eye(4), Fh(ii,:));
end

% Bracket deformation
[epsilon, delta, beta, gama] = fortz(klass, Fbrh, X1, X2, kSpr, X0);
pbrAnew = [norm(pbrhA(1:2))+epsilon, delta, pbrhA(3)+ beta];

% Replace points
LOC = L;
pB = L(end,:,1);
pBnew = repmat(pB, N, 1);

for ii = 1:N
    for i = 1:M
        if i == 1
            pAnew = RowVecTrans(Thbr, pbrAnew(ii,:));
            if ~isequal(LOC(i,:,ii), LOC(i+1,:,ii))
                LOC(i,:,ii) = pAnew;
            else
                for k = 1:C-1
                    LOC(k,:,ii) = pAnew;
                end
            end
        elseif i == M
            LOC(i,:,ii) = pBnew(ii,:);
        end
    end
end

end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(klass,Fbr,X1,X2,kSpr,X0)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change
    if nargin < 6 || isempty(X0)
        X0 = 0;
    end
    
    D      = klass.dBPA;    %uninflated diameter
    mL     = klass.Lmt-X0;     %Length of musculo-tendon
    rest   = klass.rest;    %BPA resting length
    tendon = klass.ten;     %tendon length
    fitting = klass.fitn;   %end cap length
    mif    = klass.Fm;   %maximum bpa force
    kmax   = klass.Kmax; %maximum contraction length
    KMAX   = (rest-kmax)/rest; %maximum contraction percent
    pres   = klass.P;           %pressure
    P      = pres/620;          %normalized pressure
    
    N = size(Fbr,1);
    % Normalize force vectors safely
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-8 & all(~isnan(Fbr), 2);
    
    u_hat_all = zeros(N, 3);
    u_hat_all(valid, :) = Fbr(valid, :) ./ norms(valid);
    
    % Vectorized k_b computation
    K_bracket = diag([X1, X2, X1]);       %project bracket stiffness onto force direction
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    K_rep = repmat(K_bracket, [1, 1, N]);   % [3x3xN]
    k_b = pagemtimes(pagemtimes(u_hat, K_rep), permute(u_hat, [2, 1, 3]));
    k_b = reshape(k_b, [N, 1]);
    if isinf(kSpr) || isnan(kSpr)
        k_eff = k_b;
    else
        k_eff = 1 ./ (1 ./ k_b + 1 / kSpr);  % Nx1
    end
    
    % Initialize outputs
    [e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N,1));    
    
    % Parallel root solve
    for i = 1:N
        if ~valid(i)
            continue;
        end
        
        % Per-instance constants
%         keff = k_eff(i);
        unit_vec = u_hat_all(i, :);

        if isinf(X1) && isinf(X2)
            % Rigid bracket: no deformation, optional cable stretch
            e_axial(i) = 0;
            e_bendY(i) = 0;
            e_bendZ(i) = 0;
            e_cable(i) = (klass.ten > 0) * 0.0000;  % Zero or small value if needed
            continue;
        end
        
        % Flexible case: run Newton-Raphson
        r = 0.0001;  % Initial guess
        keff = k_eff(i);

        for iter = 1:50
            contraction = (rest - (mL(i) - (tendon + r) - 2 * fitting)) / rest;
            rel = contraction / KMAX;

            fM = f_festo(rel, P, D) * mif;
            fT = keff * r;
            Fbal = fM - fT;

            % Numerical derivative
            dr = 1e-6;
            contraction_d = (rest - (mL(i) - (tendon + r + dr) - 2 * fitting)) / rest;
            rel_d = contraction_d / KMAX;
            fM_d = f_festo(rel_d, P, D) * mif;
            Fbal_d = fM_d - keff * (r + dr);
            dF = (Fbal_d - Fbal) / dr;

            %Avoid zero slope
            if abs(dF) < 1e-12 || isnan(dF)
                r = NaN;
            break;
            end
            
            % Newton-Raphson update
            r = r - Fbal / dF;

            if abs(Fbal) < 1e-6
                break;
            end
        end

        % Final force magnitude
        contraction = (rest - (mL(i) - (tendon + r) - 2 * fitting)) / rest;
        relstrain = contraction / KMAX;
        F_mag = f_festo(relstrain, P, D) * mif;

        % Bracket displacement
        e_bkt = K_bracket \ (F_mag * unit_vec');

        e_axial(i) = e_bkt(1);
        e_bendY(i) = e_bkt(2);
        e_bendZ(i) = e_bkt(3);

        % Cable elongation
        r_bracket = unit_vec * e_bkt;
        e_cable(i) = r - r_bracket;
        if e_cable(i) < 0  || klass.ten == 0
            e_cable(i) = 0;
        end
    end
end

%% ------------- Segment Lengths ------------------------
function SL = seg(klass, L_p)
C = klass.CP;
T = klass.Tk;
N = size(T, 3);
M = size(L_p, 1);
SL = zeros(N, M-1);

for ii = 1:N                    %Repeat for each orientation
    for i = 1:M-1               %Calculate all segments
        pointA = L_p(i,:,ii);
        pointB = L_p(i+1,:,ii);
        if i+1 == C
            pointB = RowVecTrans(T(:,:,ii), pointB);
        end
        SL(ii,i) = norm(pointA - pointB);
    end
end
end
               
        %% ------------- Muscle Length ------------------------
function Lmt = LMT(sL, Xi0)
% Compute muscle-tendon length from segment lengths and offset Xi0
% Xi0 can be empty [] to skip correction (i.e., when already applied)
% N = size(sL_p, 1);
Lmt = sum(sL, 2);  % Nx1, sum across segments

if ~isempty(Xi0)
    Lmt = Lmt - Xi0;
end
end           

        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
function unitD = UD(klass, L_p)
            T = klass.Tk;
            C = klass.CP;
            direction = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointA = L_p(C-1, :, i);
                pointB = L_p(C, :, i);
                direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
            end
            unitD = normalize(direction);
end
        
        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
        %For every ViaPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
function mA = Mom(klass, L_p, unitD_p)
% Moment arm = perpendicular offset from joint to line of action
C = klass.CP;
N = size(klass.Tk, 3);
mA = zeros(N, 3);

    for i = 1:N
        pointB = L_p(C,:,i);
        mA(i,:) = pointB - unitD_p(i,:) * dot(unitD_p(i,:), pointB);
    end
end     


        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
function F = Force(klass, unitD_p, strain_p)
    rest = klass.rest;
    KMAX = (rest - klass.Kmax) / rest;
    rel = strain_p ./ KMAX;  % normalized strain
    Fn = festo4(klass.dBPA, rel, klass.P);  % Force, normalized
    scalarForce = Fn .* klass.Fm;  %Redimensionalize
    scalarForce(scalarForce < 0) = 0;  % Eliminate negatives

    F = scalarForce .* unitD_p;  % Nx3
end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
function Mz = Tor(mA_p, F_p, maxF, strain_p)
N = size(F_p, 1);
Mz = zeros(N, 3);

    for i = 1:N
        if norm(F_p(i,:)) > maxF || strain_p(i,:) < -0.03
            Mz(i,:) = NaN;
        else
            Mz(i,:) = cross(mA_p(i,:), F_p(i,:));
        end
    end
end    

%% tendon springrate
function springrate = Spr(klass)
    if klass.ten > 0
        mult = 2;
        Aeff = 1.51*10^-6;%Effective area for 19-strand cable
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
    valid = norms > 1e-6 & all(~isnan(v), 2);
    vhat = zeros(N, 3);
    vhat(valid, :) = v(valid, :) ./ norms(valid);
end


end   