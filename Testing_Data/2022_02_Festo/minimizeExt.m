%% Optimize predicted torque for extensors.
function [f_all, bpa_all] = minimizeExt(Xi0, Xi1, Xi2, idx_val)
% minimizeExt: calculates predicted torque and fit metrics for a given BPA index
%
% Inputs:
%   Xi0 - extra length correction (m)
%   Xi1 - axial bracket stiffness (N/m)
%   Xi2 - bending bracket stiffness (N/m)
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
%     fprintf('Evaluating BPA #%d with [%.4f, %.2e, %.2e]\n', i, Xi0, Xi1, Xi2);
    klass_i = ke(i);
    [bpa_all(i), f_all(i,:)] = evaluateBPA(klass_i, Xi0, Xi1, Xi2);
end

end


function [bpa_i, fitvec] = evaluateBPA(klass, Xi0, Xi1, Xi2)
%% Calculate locations and properties
% fprintf('>> evaluateBPA: Xi0=%.4f, Xi1=%.2e, Xi2=%.2e\n', Xi0, Xi1, Xi2);


bpa_i = klass;
kspr = Spr(bpa_i);
Funit = computeForceVector(klass);
[L_p, gemma] = Lok(bpa_i, Xi1, Xi2,kspr, Funit);
unitD_p = UD(bpa_i, L_p);
sL_p = seg(bpa_i, L_p);
Lmt_p = LMT(bpa_i,L_p,sL_p, Xi0);
strain_p = Contraction(bpa_i, Lmt_p, gemma);
F_p = Force(bpa_i, unitD_p, strain_p);
mA_p = Mom(bpa_i, L_p, unitD_p);
M_p = Tor(mA_p, F_p, bpa_i.Fm, strain_p);

% fprintf('   Contraction range: min=%.4f max=%.4f\n', min(strain_p), max(strain_p));
% fprintf('   Original calculated torque mean: %.4f\n', mean(bpa_i.M,'omitnan'));
% fprintf('   Predicted Torque mean: %.4f\n', mean(M_p(:,3),'omitnan'));
% fprintf('   BPA resting length =%.4f and tendon length =%.4f\n', klass.rest, klass.ten);
% fprintf('   Experimental torque mean: %.4f (n=%d)\n', mean(bpa_i.Mexp,'omitnan'), numel(bpa_i.Mexp));
% fprintf('   Aexp range: %.2f to %.2f\n', min(bpa_i.Aexp), max(bpa_i.Aexp));
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
%Calculate the force unit direction from muscle origin to the next
%real point. This takes into account if there are any additional via
%points between muscle origin and muscle insertion. It also takes into
%account if a homogenous transformation matrix needs to be used to
%convert the second point into the first points frame.

L = klass.Loc;      %Location (wrapping, attachment points)
C = klass.CP;       %Cross point (moves from one frame to another)
T = klass.Tk;       %Transformation matrix
    
% Step 1: Detect the first valid segment (non-repeated)
pt1 = zeros(size(L, [3 2]));
pt2 = zeros(size(pt1));
for i = 1:size(L,3)
    for k = 2:size(L,1)
        if norm(L(k,:,i) - L(k-1,:,i)) > 1e-6  % tolerance to avoid numerical noise
            pt1(i,:) = L(k-1,:,i);
            if k == C
                pt2(i,:) = RowVecTrans(T(:,:,i),L(k,:,i));
            else
                pt2(i,:) = L(k,:,i);
            end
          break;
        end
    end
end

% Force Direction vector (hip frame)
F_vec = pt2 - pt1;
F_unit = normalize(F_vec);

end

%% ------------- Location  ------------------------
function [LOC, gama] = Lok(klass,X1,X2, kSpr, Funit)
%Point A - muscle origin
%Point B - muscle insertion
%note that Points A and B have different meanings in the other
%subfunctions.

            L = klass.Loc;      %Location (wrapping, attachment points)
            C = klass.CP;       %Cross point (moves from one frame to another)
            T = klass.Tk;       %Transformation matrix
            kmax = klass.Kmax;
            KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage 
            relstrain = klass.strain/KMAX;
            FF = festo4(klass.dBPA,relstrain,klass.P).*klass.Fm;        %Force magnitude
            FF(relstrain>1) = 0;
%             unitD = klass.unitD;            %unit direction of force vector
%             F = unitD.*FF;                  %Force vector
            Fh = Funit.*FF.*klass.Fm;     %Force vector represented in the hip frame
            
            %BPA origin bracket
            pA = L(1,:,1);                                  %Distance from hip origin to muscle insertion
            p0 = [-6.26, -29.69, 75.06] / 1000;       %from hip origin to lower bolt hole on superior anterior bracket of the Bifemsh_Pam
            y0 = [-47.48, -35.63, 75.06] / 1000;             % y reference, lower bolt hole on superior posterior bracket
            z0 = [-8.39, -14.85, 75.06] / 1000;              % z reference, upper bolt hole on superior anterior bracket
            
            %Compute direction vectors
            z_hat = (p0 - z0) / norm(p0 - z0);               % z' direction
            y_hat_prelim = (p0 - y0) / norm(p0 - y0);        % y' approximate
            
            %Project pA-p0 onto z' to ensure its on x'-y' plane
            v = pA - p0;
            lambda = -dot(v, z_hat);                        %project v onto z'
            Pbr_A = p0 + lambda * z_hat;                      % bracket frame origin            
            
            x_hat = cross(y_hat_prelim, z_hat);     %Orthogonal to y_hat_prelim and z_hat
            x_hat = x_hat / norm(x_hat);            %unit vector
            y_hat = cross(z_hat, x_hat);            %Recomputed to ensure right-handedness

            Rhbr = [x_hat(:), y_hat(:), z_hat(:)];  %Rotation matrix for bracket frame described in hip frame
            Thbr = RpToTrans(Rhbr, Pbr_A');         %Transformation matrix, represent bracket frame in hip frame 
            pbrA = pA-Pbr_A;                        %vector from bracket to point A (i.e. Point A in the bracket frame)
            %BPA insertion bracket
            pB = L(end,:,1);                  %Distance from knee frame to muscle insertion
%             Pbr_B = [27.61, -125.91, -0.54]/1000;    %from knee ICR to extensor insertion bracket (where it starts to cantilever)
%             pkbrB = pB-Pbr_B;                  %vector from bracket to point B, in the knee frame
%             thetabrB = atan2(pkbrB(2),pkbrB(1));   %angle between pbrB and x axis
%             Rkbr = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
%                    sin(thetabrB) cos(thetabrB) 0; ...
%                    0    0   1];
%             Tkbr = RpToTrans(Rkbr, Pbr_B');    %Transformation matrix, flexor bracket frame in knee frame            
            
            %setup arrays for loop
            LOC = L;            %new location matrix            
            N = size(L,3);
            M = size(L,1);
%             Fh = zeros(N,3);
            Fbrh = zeros(N,3);
            pAnew = zeros(N,3);     %New point A, in the hip frame
%             Fbrk = zeros(N,3);       %Force vector represented in the tibia bracket frame  
%             pBnew = zeros(N,3);
            pBnew = repmat(pB, N, 1); %Force pBnew = pB across all orientations
%             gama = zeros(N,1);  % Ensure gamma is returned even if bracket B is rigid

            %Transform force vector into each bracket's frame
            for ii = 1:N                          %Repeat for each orientation
%                 Fh(ii,:) = -RowVecTrans(T(:,:,ii),F(ii,:));       %Force vector represented in the hip frame
                Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:));   %Force vector in the hip frame represented in the bracket frame    
%                 Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),F(ii,:));    %Force vector in the tibia frame represented in the bracket frame                    
            end
            [epsilon1, delta1, beta1, gama] = fortz(klass,Fbrh,X1,X2,kSpr);  %strain from force divided by tensile stiffness
            pbrAnew = [pbrA(1)+epsilon1, pbrA(2)+delta1, pbrA(3)+beta1]; %New point A, represented in the bracket frame
            
% --- (Optional) Compute deformation for insertion bracket
% [epsilon2, delta2, beta2, gama] = fortz(klass, Fbrk, X1, X2, kSpr);
% pbrBnew = [norm(pkbrB(1:2)) + epsilon2, delta2, pkbrB(3) + beta2];
            
            for ii = 1:N                          %Repeat for each orientation
                for i = 1:M                      %Repeat for all muscle segments
                    if i == 1
                        pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
                        if ~isequal(LOC(i,:,ii), LOC(i+1,:,ii))
                            LOC(i,:,ii) = pAnew(ii,:);      %Update location matrix
                        else
                            for k = 1:C-1
                                LOC(k,:,ii) = pAnew(ii,:);
                            end
                        end
                    elseif i == M
%                         pBnew(ii,:) = RowVecTrans(Tkbr, pbrBnew(ii,:));     %New point B, in the tibia frame
                        LOC(i,:,ii) = pBnew(ii,:);
                    else
                    end

                end
            end
            
end

%% ------------- Segment Lengths ------------------------
function SL = seg(klass, L_p)
            C = klass.CP;
            T = klass.Tk;
            SL = zeros(size(T, 3), size(L_p, 1) - 1);
            
            for ii = 1:size(T, 3)                          %Repeat for each orientation
                for i = 1:size(L_p, 1)-1                      %Repeat for all muscle segments
                    pointA = L_p(i, :,ii);
                    pointB = L_p(i+1, :,ii);
                    if i+1 == C
                        pointB = RowVecTrans(T(:, :, ii), pointB);
                    end
                    SL(ii, i) = norm(pointA - pointB);
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
            unitD = zeros(size(direction));
            
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
            T = klass.Tk;
            C = klass.CP;
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L_p(C, :, i);
                mA(i, :) = pointB - unitD_p(i, :)*dot(unitD_p(i, :), pointB);
            end
end       
        
        %% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass,Lmt_p,gama)
            rest = klass.rest;
            tendon = klass.ten;
            fitting = klass.fitn;
%             theta_k = klass.Ak;
            
            adjusted_Lmt = Lmt_p - (tendon + gama + delta_L) - 2*fitting;
            contraction = (rest - adjusted_Lmt) / rest;
end

        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
function F = Force(klass, unitD_p, strain_p)
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
            
           rel = strain_p./KMAX;                    %relative strain        
           
           Fn = festo4(klass.dBPA,rel,klass.P);

           scalarForce = Fn.*klass.Fm;
           scalarForce(scalarForce < 0) = 0;
%            scalarForce(scalarForce > maxF) = NaN;            
            
            F = scalarForce.*unitD_p;

end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
function Mz = Tor(mA_p, F_p, maxF, strain_p)  
            Mz = zeros(size(F_p));
           
            for i = 1:size(F_p, 1)
                if norm(F_p(i,:)) > maxF
                    Mz(i,:) = NaN;
                elseif strain_p(i,:) < -0.03
                    Mz(i,:) = NaN;
                else
                    Mz(i, :) = cross(mA_p(i, :), F_p(i, :));
                end
            end


end      

%% tendon springrate
function springrate = Spr(klass)
    if klass.ten > 0
        mult = 1;
        Aeff = 1.51*10^-6;%Effective area for 19-strand cable
        E = 193*10^9;       %Young's Modulus
        L = klass.ten;      %tendon length        
        springrate = mult*Aeff*E/L;
    else
        springrate = Inf;
    end
        
end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(klass,Fbr,X1,X2,kSpr)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change

    D      = klass.dBPA;
    mL     = klass.Lmt;
    rest   = klass.rest;
    tendon = klass.ten;
    fitting = klass.fitn;
    mif    = klass.Fm;
    kmax   = klass.Kmax;
    KMAX   = (rest-kmax)/rest;
    pres   = klass.P;
    P      = pres/620;
    
    N = size(Fbr,1);
    % Normalize force vectors safely
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-8 & all(~isnan(Fbr), 2);
    
    u_hat_all = zeros(N, 3);
    u_hat_all(valid, :) = Fbr(valid, :) ./ norms(valid);
    
    % Vectorized k_b computation
    K_bracket = diag([X2, X2, X1]);       %project bracket stiffness onto force direction
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
    end
end

%% Subfunctions
function t = SSE(klass, M_p)
     Mpredict1 = M_p;
     Mpredict2 = griddedInterpolant(klass.Ak,Mpredict1(:,3));
     M_opt = Mpredict2(klass.Aexp);
     [RMSE, fvu, maxResid] = Go_OfF(klass.Mexp,M_opt);
     t = [RMSE, fvu, maxResid];
end

function vhat = normalize(v)
    N = size(v,1);
    norms = vecnorm(v,2,2);
    valid = norms > 1 & all(~isnan(v), 2);
    vhat = zeros(N, 3);
    vhat(valid, :) = v(valid, :) ./ norms(valid);
end


end   