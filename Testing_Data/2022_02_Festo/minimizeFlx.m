%% Optimize predicted torque for extensors.
function [f, varargout] = minimizeFlx(Xi0,Xi1,Xi2)

%% load
%kf = knee flexor, kf(1) = pinned joint, kf(2) = biomimetic;
%ke = knee extensor, same as above
%kf.L := lengths = [42 42] cm   where 1st length is for 10mm diameter,
%2nd length is for 20mm diameter.
%example: kf(2).L(1).Mz z-axis torque for biomimetic knee, flexor, 42cm length
%'exp' suffix means experimental
%'_h' suffix means hybrid
%'_p' suffix means prime, as in the new prediction values

load FlxBioBPASet.mat kf %This loads the following, which was ran and saved:

%     load KneeFlx_10mm_42cm.mat Bifemsh_Pam phiD
%     Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeFlx_10mm_42cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([Angle', Torque', InflatedLength', ICRtoMuscle', TorqueHand']);
%     kf(1) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
%                   'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
%                   'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
%                   'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
%                   'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
%                   'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
%                   'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
% 
%     % 42cm length, 20mm diameter, 620 pressure.
%     load KneeFlx_20mm_42cm.mat Bifemsh_Pam3 phiD
%     Bifemsh_Pam = Bifemsh_Pam3;
%     Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeFlx_20mm_42cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([Angle', Torque, InflatedLength', ICRtoMuscle', TorqueHand']);
%     kf(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
%                   'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
%                   'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
%                   'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
%                   'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
%                   'M',Bifemsh_Pam.Torque(:,3),'Aexp',A([1:10, 12],1),'Mexp',A([1:10, 12],2),...
%                   'A_h',A([1:10, 12],1),'Lm_h',A([1:10, 12],3),'mA_h',A([1:10, 12],4),'M_h',A([1:10, 12],5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Bifemsh_Pam3 phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
%     
%     % 42cm length, 20mm diameter, 325 pressure.
%     load KneeFlx_20mm_42cm.mat Bifemsh_Pam2 phiD
%     Bifemsh_Pam = Bifemsh_Pam2;
%     Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeFlx_20mm_42cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([Angle', Torque, InflatedLength', ICRtoMuscle', TorqueHand']);
%     kf(3) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
%                   'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
%                   'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
%                   'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
%                   'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
%                   'M',Bifemsh_Pam.Torque(:,3),'Aexp',A([17, 20:22],1),'Mexp',A([17, 20:22],2),...
%                   'A_h',A([17, 20:22],1),'Lm_h',A([17, 20:22],3),'mA_h',A([17, 20:22],4),'M_h',A([17, 20:22],5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Bifemsh_Pam2 phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A


%% RMSE, fvu, and Max Residual
f = NaN(1,3);
g = NaN(1,3);
a = 3;
h = cell(1,3);

% Use first loaded struct as template
klaus_template = kf(1);  

% Initialize struct arrays with proper shape and fields
klaus = repmat(klaus_template, 1, a);
bpa   = repmat(klaus_template, 1, a);

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
    klaus(j) = kf(j);
    % Calculate locations and properties
    kspr{j} = Spr(klaus(j));
    [L_p{j}, gemma{j}] = Lok(klaus(j), Xi1, Xi2,kspr{j});
    unitD_p{j} = UD(klaus(j), L_p{j});
    sL_p{j} = seg(klaus(j), L_p{j});
    Lmt_p{j} = LMT(sL_p{j}, Xi0);
    strain_p{j} = Contraction(klaus(j), Lmt_p{j},gemma{j});
    F_p{j} = Force(klaus(j), unitD_p{j}, strain_p{j});
    mA_p{j} = Mom(klaus(j), L_p{j}, unitD_p{j});
    M_p{j} = Tor(mA_p{j}, F_p{j}, klaus(j).Fm, strain_p{j});
    
    % Package into output struct
    bpa(j) = klaus(j);
    bpa(j).Lmt_p = Lmt_p{j};
    bpa(j).mA_p = mA_p{j};
    bpa(j).M_p = M_p{j};
    bpa(j).F_p = F_p{j};
    bpa(j).strain_p = strain_p{j};
    bpa(j).L_p = L_p{j};
    bpa(j).gama = gemma{j};

    % SSE calculation
    h{j} = SSE(klaus(j), M_p{j});
end


    f = h{1}; % 10mm biomimetic knee fit values
    if nargout > 1
        varargout{1} = h{2}; % 20mm biomimetic knee fit values, full pressure
        varargout{2} = h{3}; % 20mm biomimetic knee fit values, full pressure
        varargout{3} = bpa;  % Full structure with prediction info
    end

%% Nested functions, mostly modified from MonoPamExplicit

%% ------------- Location  ------------------------
function [LOC, gamma] = Lok(klass,X1,X2,kSpr)
            L = klass.Loc;      %Location (wrapping, attachment points)
            C = klass.CP;       %Cross point (moves from one frame to another)
            T = klass.Tk;       %Transformation matrix, tibia frame represented in the hip frame
            kmax = klass.Kmax;  
            KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage 
            FF = festo4(klass.dBPA,klass.strain/KMAX,klass.P).*klass.Fm;       %Force magnitude
%             unitD = klass.unitD;                                                %Unit direction of force, tibia frame
%             F = unitD.*FF;                                            %Force vector, tibia frame
            pA = L(1,:,1);                                  %Distance from hip origin to muscle insertion
%             Pbr = [-0.8100  -20.222   31.66]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
            Pbr = [9.48  -36.15   30.27]/1000;       %from hip origin to bracket bolt pattern centroid
            phbrA = pA-Pbr;                                  %vector from bracket to point A (in the hip frame)
            thetabrA = atan2(phbrA(2),phbrA(1));   %angle between pbrA and x axis
            RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...     %Rotation matrix
                   sin(thetabrA) cos(thetabrA) 0; ...
                   0    0   1];
            pbrhA = RhbrZ'*phbrA';       %Vector in the bracket frame
            % Now calculate angle from x-axis to this vector
            thetaY = atan2(pbrhA(3), pbrhA(1));  % z vs x (in bracket frame)

            % Rotation matrix about y-axis (local frame adjustment)
            Ry = [cos(thetaY)  0  sin(thetaY);
                  0                1  0;
                 -sin(thetaY) 0  cos(thetaY)];
            Rhbr = RhbrZ*Ry;            %Rotate about y-axis in body frame
            Thbr = RpToTrans(RhbrZ, Pbr');    %Transformation matrix, represent bracket frame in hip frame              
            LOC = L;
            N = size(L,3);
            M = size(L,1);
            Fh = zeros(N,3);
            Fbrh = zeros(N,3);
            pAnew = zeros(N,3);     %New point A, in the hip frame
            for ii = 1:N                          %Repeat for each orientation
                        Fh(ii,:) = norm(RowVecTrans(T(:,:,ii),L(C,:,ii))-L(C-1,:,ii)).*FF(ii,1);               %Force vector represented in the hip frame
                        Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:));            %Force vector in the hip frame represented in the bracket frame
            end
            [epsilon, delta, beta, gamma] = fortz(klass,Fbrh,X1,X2,kSpr);  %strain from force divided by tensile stiffness
            pbrAnew = [norm(pbrhA(1:2))+epsilon, delta, pbrhA(3)+ beta]; %New point A, represented in the bracket frame
                        

            for ii = 1:N                          %Repeat for each orientation
                for i = 1:M                       %Repeat for all muscle segments
                    if i == C-1                  
                        pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
                        LOC(i,:,ii) = pAnew(ii,:);      %Update location matrix
                    else
                    end

                end
            end
end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(klass,Fbr,X1,X2,kSpr)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change

    D = klass.dBPA;         %BPA diameter
    sL = seg(klass, klass.Loc);
    Lmt = LMT(sL,[]);         %musculotendon length
    rest = klass.rest;      %resting length
    tendon = klass.ten;     %tendon length
    fitn = klass.fitn;    %fitting length
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
    K_bracket = diag([X1, X2, X1]);       %project bracket stiffness onto force direction
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    K_rep = repmat(K_bracket, [1, 1, N]);   % [3x3xN]
    k_b = pagemtimes(pagemtimes(u_hat, K_rep), permute(u_hat, [2, 1, 3]));
    k_b = reshape(k_b, [N, 1]);
    k_eff = 1 ./ (1 ./ k_b + 1 / kSpr);  % Nx1
    
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
        
        
        if isinf(X1) && isinf(X2)
            % Rigid body: no bracket deformation
            e_axial(i) = 0;
            e_bendY(i) = 0;
            e_bendZ(i) = 0;
            e_cable(i) = r;  % All elongation goes to cable
            continue;
        end

        % Solve for r (deflection) using fzero
    % First, compute contraction from Xi0 if deflection is zero 
    contraction0    = ( rest ...
                   - ( Lmt(i) - (tendon + 0) - 2*fitn ) ) ...
                   / rest;
    relstrain0      = contraction0 / KMAX;  %relative strain
    if relstrain0 >= 1
        r = 0;
%     elseif contraction0 <= -0.03
%         r = NaN;
    else
    
        relfun = @(r) ...
            f_festo( ...
                (rest - (Lmt(i) - (tendon + r) - 2*fitn)) / rest / KMAX, ...
                P, D ...
            ) * mif - keff * r;

        try
            r = fzero(relfun, [0, 0.2]);
        catch
            r = 0;
        end
        r = max(r,0); %guard against r being slightly negative.
    end

    % Final force magnitude
        contraction = (rest - (Lmt(i) - (tendon + r) - 2 * fitn)) / rest;
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
        L = klass.ten;      %tendon length
%         L = klass.ten+.015;      %tendon length
        
        springrate = mult*Aeff*E/L;
%         springrate = Inf;
        
end



%% Subfunctions
function t = SSE(klass, M_p)
     [Ak_sorted, idx] = sort(klass.Ak);
     M_sorted = M_p(idx, 3);
     Mpredict2 = griddedInterpolant(Ak_sorted, M_sorted);
     M_opt = Mpredict2(klass.Aexp);
     [RMSE, fvu, maxResid] = Go_OfF(klass.Mexp,M_opt);
     t = [RMSE, fvu, maxResid];
% --- New: Focused SSE where |Mexp| < 3 Nm ---
%     lowTorqueMask = abs(klass.Mexp) < 3;
%     if any(lowTorqueMask)
%         SSE_low = mean((klass.Mexp(lowTorqueMask) - M_opt(lowTorqueMask)).^2,'omitnan');
%     else
%         SSE_low = 0; % Or NaN if you'd rather it not contribute
%     end
% 
%     t = [RMSE, fvu, SSE_low];  % Replace maxResid with low-torque SSE
end

end