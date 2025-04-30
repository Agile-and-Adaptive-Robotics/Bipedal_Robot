%% Optimize predicted torque for extensors.
function [f, varargout] = minimizeExt(Xi0,Xi1,Xi2)

%% load
%kf = knee flexor, kf(1) = pinned joint, kf(2) = biomimetic;
%ke = knee extensor, same as above
%ke.L := lengths = [42 42 46 48] cm
%example: ke(1).L(3).Mz z-axis torque for pinned knee, flexor, 46cm length
%'exp' suffix means experimentally measured
%'_h' suffix means hybrid calculation
%'_p' suffix means prime, as in the new prediction values

% 42cm length, no tendon
    load KneeExtPin_10mm_all.mat Vas_Pam_42cm phiD
    Ma = Vas_Pam_42cm.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeExtPin10mm_42cm.mat Angle0 Torque0 InflatedLength0 ICRtoMuscle0 TorqueHand0
    A = sortrows([Angle0, Torque0, InflatedLength0, ICRtoMuscle0, TorqueHand0]);
    ke(1).L(1) = struct('Ak',phiD,'Loc',Vas_Pam_42cm.Location,'CP',Vas_Pam_42cm.Cross,'dBPA',Vas_Pam_42cm.Diameter, ...
                  'Tk',Vas_Pam_42cm.TransformationMat,'rest',Vas_Pam_42cm.RestingL,'Kmax',Vas_Pam_42cm.Kmax,...
                  'fitn',Vas_Pam_42cm.FittingLength,'ten',Vas_Pam_42cm.TendonL,'P',Vas_Pam_42cm.Pressure, ...
                  'Lmt',Vas_Pam_42cm.MuscleLength,'strain',Vas_Pam_42cm.Contraction, 'unitD',Vas_Pam_42cm.UnitDirection, ...
                  'mA',G,'Fm',Vas_Pam_42cm.Fmax,'F',Vas_Pam_42cm.Force, 'seg',Vas_Pam_42cm.SegmentLengths, ...
                  'M',Vas_Pam_42cm.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
    clear Vas_Pam_42cm phiD Ma G Angle0 Torque0 InflatedLength0 ICRtoMuscle0 TorqueHand0 A

% 42cm length, tendon
    load KneeExtPin_10mm_all.mat Vas_Pam_42cm_tendon phiD
    Ma = Vas_Pam_42cm_tendon.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeExtPin10mm_42cm.mat Angle1 Torque1 InflatedLength1 ICRtoMuscle1 TorqueHand1
    A = sortrows([Angle1, Torque1, InflatedLength1, ICRtoMuscle1, TorqueHand1]);
    ke(1).L(2) = struct('Ak',phiD,'Loc',Vas_Pam_42cm_tendon.Location,'CP',Vas_Pam_42cm_tendon.Cross,'dBPA',Vas_Pam_42cm_tendon.Diameter, ...
                  'Tk',Vas_Pam_42cm_tendon.TransformationMat,'rest',Vas_Pam_42cm_tendon.RestingL,'Kmax',Vas_Pam_42cm_tendon.Kmax,...
                  'fitn',Vas_Pam_42cm_tendon.FittingLength,'ten',Vas_Pam_42cm_tendon.TendonL,'P',Vas_Pam_42cm_tendon.Pressure, ...
                  'Lmt',Vas_Pam_42cm_tendon.MuscleLength,'strain',Vas_Pam_42cm_tendon.Contraction, 'unitD',Vas_Pam_42cm_tendon.UnitDirection, ...
                  'mA',G,'Fm',Vas_Pam_42cm_tendon.Fmax,'F',Vas_Pam_42cm_tendon.Force, 'seg',Vas_Pam_42cm_tendon.SegmentLengths, ...
                  'M',Vas_Pam_42cm_tendon.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
    clear Vas_Pam_42cm_tendon phiD Ma G Angle1 Torque1 InflatedLength1 ICRtoMuscle1 TorqueHand1 A
    

% 46cm length
    load KneeExtPin_10mm_all.mat Vas_Pam_46cm phiD
    Ma = Vas_Pam_46cm.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeExtPin10mm_46cm.mat AngleX Torque InflatedLength ICRtoMuscle TorqueHand
    A = sortrows([AngleX, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    ke(1).L(3) = struct('Ak',phiD,'Loc',Vas_Pam_46cm.Location,'CP',Vas_Pam_46cm.Cross,'dBPA',Vas_Pam_46cm.Diameter, ...
                  'Tk',Vas_Pam_46cm.TransformationMat,'rest',Vas_Pam_46cm.RestingL,'Kmax',Vas_Pam_46cm.Kmax,...
                  'fitn',Vas_Pam_46cm.FittingLength,'ten',Vas_Pam_46cm.TendonL,'P',Vas_Pam_46cm.Pressure, ...
                  'Lmt',Vas_Pam_46cm.MuscleLength,'strain',Vas_Pam_46cm.Contraction, 'unitD',Vas_Pam_46cm.UnitDirection, ...
                  'mA',G,'Fm',Vas_Pam_46cm.Fmax,'F',Vas_Pam_46cm.Force, 'seg',Vas_Pam_46cm.SegmentLengths, ...
                  'M',Vas_Pam_46cm.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
    clear Vas_Pam_46cm phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

    load KneeExtPin_10mm_all.mat Vas_Pam_48cm phiD
    Ma = Vas_Pam_48cm.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeExtPin10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    ke(1).L(4) = struct('Ak',phiD,'Loc',Vas_Pam_48cm.Location,'CP',Vas_Pam_48cm.Cross,'dBPA',Vas_Pam_48cm.Diameter, ...
                  'Tk',Vas_Pam_48cm.TransformationMat,'rest',Vas_Pam_48cm.RestingL,'Kmax',Vas_Pam_48cm.Kmax,...
                  'fitn',Vas_Pam_48cm.FittingLength,'ten',Vas_Pam_48cm.TendonL,'P',Vas_Pam_48cm.Pressure, ...
                  'Lmt',Vas_Pam_48cm.MuscleLength,'strain',Vas_Pam_48cm.Contraction, 'unitD',Vas_Pam_48cm.UnitDirection, ...
                  'mA',G,'Fm',Vas_Pam_48cm.Fmax,'F',Vas_Pam_48cm.Force, 'seg',Vas_Pam_48cm.SegmentLengths, ...
                  'M',Vas_Pam_48cm.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
    clear Vas_Pam_48cm phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A


%% What to optimize based on selection
% M_opt = zeros(size(kf(1).L(3).Mexp)); %Optimized torque prediction

%% RMSE, fvu, and Max Residual
f = NaN(1,3);
g = NaN(1,3);
a = 4;
h = cell(1,2);

klaus_temp = ke(1).L(3);   % Grab a template struct
% Initialize prediction fields to empty
klaus_temp.L_p   = [];
klaus_temp.Lmt_p = [];
klaus_temp.strain_p   = [];
klaus_temp.gamma   = [];
klaus_temp.mA_p  = [];
klaus_temp.F_p   = [];
klaus_temp.M_p   = [];

% Initialize struct arrays with proper shape and fields
klaus = repmat(klaus_temp, 1, a);
bpa   = repmat(klaus_temp, 1, a);

L_p = cell(1,a);
gemma = cell(1,a);
kspr = cell(1,a);
unitD_p = cell(1,a);
sL_p = cell(1,a);
Lmt_p = cell(1,a);
strain_p = cell(1,a);
F_p = cell(1,a);
mA_p = cell(1,a);
M_p = cell(1,a);

for j = 1:a
    klaus(j) = ke(1).L(j);
    % Calculate locations and properties
    kspr{j} = Spr(klaus(j));
    [L_p{j}, gemma{j}] = Lok(klaus(j), Xi1, Xi2,kspr{j});
    unitD_p{j} = UD(klaus(j), L_p{j});
    sL_p{j} = seg(klaus(j), L_p{j});
    Lmt_p{j} = LMT(klass(j),L_p,sL_p, Xi0);
    strain_p{j} = Contraction(klaus(j), Lmt_p{j}, gemma{j});
    F_p{j} = Force(klaus(j), unitD_p{j}, strain_p{j});
    mA_p{j} = Mom(klaus(j), L_p{j}, unitD_p{j});
    M_p{j} = Tor(mA_p{j}, F_p{j}, klaus(j).Fm, strain_p{j});
    
    % Package into output struct
    bpa(j) = klaus(j);
    bpa(j).L_p = L_p{j};
    bpa(j).Lmt_p = Lmt_p{j};
    bpa(j).mA_p = mA_p{j};
    bpa(j).F_p = F_p{j};
    bpa(j).M_p = M_p{j};
    bpa(j).strain_p   = strain_p{j};
    bpa(j).gamma   = gemma{j};

    % SSE calculation
    h{j} = SSE(klaus(j), M_p{j});
end


    f = h{1}; % Optimization fit values
    if nargout > 1
        varargout{1} = h{2}; % Additional BPA fit values
        varargout{2} = h{3}; % Additional BPA fit values
        varargout{3} = h{4}; % Additional BPA fit values
        varargout{4} = bpa;  % Full structure with prediction info
    end


%% Nested functions, modified from MonoPamExplicit
        %% ------------- Location  ------------------------
function [LOC, gamma] = Lok(klass,X1,X2, kSpr)
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
            unitD = klass.unitD;            %unit direction of force vector
            F = unitD.*FF;                  %Force vector
            
            %BPA origin bracket
            pA = L(1,:,1);                                  %Distance from hip origin to muscle insertion
            Pbr_A = [-0.8100  -20.2650   32.2100]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
            pbrA = pA-Pbr_A;                                  %vector from bracket to point A (in the hip frame)
            thetabrA = norm(wrapToPi(acos(dot([pbrA(1) pbrA(2) 0],[1,0,0])/(norm([pbrA(1) pbrA(2) 0])))));   %angle between pbrA and x axis
            Rhbr = [cos(thetabrA) -sin(thetabrA) 0; ...     %Rotation matrix
                   sin(thetabrA) cos(thetabrA) 0; ...
                   0    0   1];
            Thbr = RpToTrans(Rhbr, Pbr_A');    %Transformation matrix, represent bracket frame in hip frame 
            
            %BPA insertion bracket
            pB = L(end,:,1);                  %Distance from knee frame to muscle insertion
            Pbr_B = [0.2575, -104.25, 0];    %from knee ICR to extensor insertion bracket (where it starts to cantilever)
            pkbrB = pB-Pbr_B;                  %vector from bracket to point B, in the knee frame
            thetabrB = wrapToPi(acos(dot([pkbrB(1) pkbrB(2) 0],[1,0,0])/(norm([pkbrB(1) pkbrB(2) 0]))));   %angle between pbrB and x axis
            Rkbr = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
                   sin(thetabrB) cos(thetabrB) 0; ...
                   0    0   1];
            Tkbr = RpToTrans(Rkbr, Pbr_B');    %Transformation matrix, flexor bracket frame in knee frame            
            
            %setup arrays for loop
            LOC = L;            %new location matrix            
            N = size(L,3);
            M = size(L,1);
            Fh = zeros(N,3);
            Fbrh = zeros(N,3);
            pAnew = zeros(N,3);     %New point A, in the hip frame
            Fbrk = zeros(N,3);       %Force vector represented in the tibia bracket frame  
            pBnew = zeros(N,3);
            
            %Transform force vector into each bracket's frame
            for ii = 1:N                          %Repeat for each orientation
                Fh(ii,:) = -RowVecTrans(T(:,:,ii),F(ii,:));       %Force vector represented in the hip frame
                Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:));   %Force vector in the hip frame represented in the bracket frame    
                Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),F(ii,:));    %Force vector in the tibia frame represented in the bracket frame                    
            end
            [epsilon1, delta1, beta1, ~] = fortz(klass,Fbrh,X1,X2,kSpr);  %strain from force divided by tensile stiffness
            pbrAnew = [norm([pbrA(1) pbrA(2)])+epsilon1, delta1, pbrA(3)+beta1]; %New point A, represented in the bracket frame
            
            [epsilon2, delta2, beta2, gamma] = fortz(klass,Fbrk,X1,X2);  %strain from force divided by tensile stiffness
            pbrBnew = [norm([pkbrB(1) pkbrB(2)])+epsilon2, delta2, pkbrB(3)+beta2]; %new point B, in the bracket frame
            
            for ii = 1:N                          %Repeat for each orientation
                for i = 1:M                      %Repeat for all muscle segments
                    if i == 1
                        pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
                        LOC(i,:,ii) = pAnew(ii,:);      %Update location matrix                        
                    elseif i == M
                        pBnew(ii,:) = RowVecTrans(Tkbr, pbrBnew(ii,:));     %New point B, in the tibia frame
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
        %Function that calculates the musclutendon length
function Lmt = LMT(klass,L_p,sL_p, Xi0)
            T = klass.Tk;
            Lmt = zeros(size(T, 3), 1);
            
            for ii = 1:size(Lmt, 1)                          %Repeat for each orientation
                for i = 1:size(L_p, 1, 1)-1                      %Repeat for all muscle segments
                    Lmt(ii, 1) = Lmt(ii, 1) + sL_p(ii, i);
                end
            end
            Lmt = Lmt - Xi0;
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
                unitD(i, :) = direction(i, :)/norm(direction(i, :));
            end
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
function contraction = Contraction(klass,Lmt_p,gamma)
            rest = klass.rest;
            tendon = klass.ten;
            fitting = klass.fitn;
            
            contraction = (rest-(Lmt_p-(tendon+gamma)-2*fitting))/rest;    %(minus Xi0 is used in LMT function, above)
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
        switch klass.dBPA
            case 10
                mult = 2;
            case 20
                mult = 6;
        end
        Aeff = 1.51*10^-6;%Effective area for 19-strand cable
        E = 193*10^9;       %Young's Modulus
        L = klass.ten+.015;      %tendon length
        
        springrate = mult*Aeff*E/L;
%         springrate = Inf;
        
end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(klass,Fbr,X1,X2,kSpr)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change

    D = klass.dBPA;         %BPA diameter
    mL = klass.Lmt;         %musculotendon length
    rest = klass.rest;      %resting length
    tendon = klass.ten;     %tendon length
    fitting = klass.fitn;    %fitting length
    mif = klass.Fm;         %maximum force
    kmax = klass.Kmax;      %maximum contracted length
    KMAX = (rest-kmax)/rest; %turn it into a percentage
    pres = klass.P;         %pressure
    P = pres/620;            %normalized pressure
    
    N = size(Fbr,1);
    % Normalize force vectors safely
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-8 & all(~isnan(Fbr), 2);
    u_hat_all = zeros(N, 3);
    u_hat_all(valid, :) = Fbr(valid, :) ./ norms(valid);
    
    % Vectorized k_b computation
    K_bracket = diag([X1, X2, X2]);       %project bracket stiffness onto force direction
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    K_rep = repmat(K_bracket, [1, 1, N]);   % [3x3xN]
    k_b = pagemtimes(pagemtimes(u_hat, K_rep), permute(u_hat, [2, 1, 3]));
    k_b = reshape(k_b, [N, 1]);
    if isinf(kSpr) || isnan(kSpr)
    k_eff = k_b;
    else
    k_eff = 1 ./ (1 ./ k_b + 1 / kSpr);  % Nx1
    end
    
    % Allocate outputs
    e_axial = zeros(N, 1);
    e_bendY = zeros(N, 1);
    e_bendZ = zeros(N, 1);
    e_cable = zeros(N, 1);    
    
    % Parallel root solve
    for i = 1:N
        if ~valid(i)
            continue;
        end
        
        % Per-instance constants
        keff = k_eff(i);
        unit_vec = u_hat_all(i, :);

        
        r = 0.0001;     %initial guess
        
        
        if isinf(X1) || isinf(X2)
            % Rigid body: no bracket deformation
            e_axial(i) = 0;
            e_bendY(i) = 0;
            e_bendZ(i) = 0;
            if klass.ten > 0
            e_cable(i) = r;  % All elongation goes to cable
            else
            e_cable(i) = 0;
            end
        continue;
        end

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




end   