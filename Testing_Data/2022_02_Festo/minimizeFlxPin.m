%% Optimize predicted torque for extensors.
function [f, varargout] = minimizeFlxPin(Xi0,Xi1,Xi2)

%% load
%kf = knee flexor, kf(1) = pinned joint, kf(2) = biomimetic;
%ke = knee extensor, same as above
%kf.L := lengths = [42 46 48] cm
%example: kf(1).L(2).Mz z-axis torque for pinned knee, flexor, 46cm length
%'exp' suffix means experimental
%'_h' suffix means hybrid
%'_p' suffix means prime, as in the new prediction values

load FlxPinBPASet.mat kf %This loads the following, which was ran and saved:

%     load KneeFlxPin_10mm_48cm.mat Bifemsh_Pam phiD
%     Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
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
% % if nargout > 1
%     % 46cm length
%     load KneeFlxPin_10mm_46cm.mat Bifemsh_Pam phiD
%     Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
%     G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
%     load Plot_KneeFlxPin_10mm_46cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
%     A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
%     kf(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
%                   'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
%                   'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
%                   'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
%                   'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
%                   'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
%                   'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
%                   'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'strain_p',[],'L_p',[],'gama',[]);
%     clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
% % end


%% What to optimize based on selection
% M_opt = zeros(size(kf(1).L(3).Mexp)); %Optimized torque prediction

%% RMSE, fvu, and Max Residual
f = NaN(1,3);
g = NaN(1,3);
a = 1 + ( nargout>1 );
h = cell(1,2);

klaus_temp = kf(1);   % Grab a template struct
klaus = repmat(klaus_temp, 1, a);  % ✅ GOOD INIT
bpa = repmat(klaus_temp, 1, a);  % ✅ GOOD INIT

L_p = cell(1,a);
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
    L_p{j} = Lok(klaus(j), Xi1, Xi2, []);   %Bracket deformation changing geometry
    unitD_p{j} = UD(klaus(j), L_p{j});   %New force direction
    sL_p{j} = seg(klaus(j), L_p{j});   %New segment lengths
    Lmt_p{j} = LMT(sL_p{j}, Xi0);     %New musclulotendon length
    strain_p{j} = Contraction(klaus(j), Lmt_p{j});  %new contraction amount
    F_p{j} = Force(klaus(j), unitD_p{j}, strain_p{j});  %new force amount
    mA_p{j} = Mom(klaus(j), L_p{j}, unitD_p{j});   %new moment arm
    M_p{j} = Tor(mA_p{j}, F_p{j}, klaus(j).Fm, strain_p{j});  %new torque
    
    % Package into output struct
    bpa(j) = klaus(j);
    bpa(j).Lmt_p = Lmt_p{j};
    bpa(j).mA_p = mA_p{j};
    bpa(j).M_p = M_p{j};
    bpa(j).F_p = F_p{j};
    bpa(j).strain_p = strain_p{j};
    bpa(j).L_p = L_p{j};

    % SSE calculation
    h{j} = SSE(klaus(j), M_p{j});
end


    f = h{1}; % Optimization fit values
    if nargout > 1
        varargout{1} = h{2}; % Validation fit values
        varargout{2} = bpa;  % Full structure with prediction info
    end


%% Nested functions, modified from MonoPamExplicit
        %% ------------- Location  ------------------------
function LOC = Lok(klass,X1,X2,X0)
            L = klass.Loc;      %Location (wrapping, attachment points)
            C = klass.CP;       %Cross point (moves from one frame to another)
            T = klass.Tk;       %Transformation matrix
            kmax = klass.Kmax;  %max contracted length
            KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage
            sL = seg(klass, L);
            if isempty(X0)
                X0=0;
            end
            Lmt = LMT(sL, X0);   %musculotendon length
            strain = Contraction(klass, Lmt);
            relstrain = strain/KMAX;
            FF = festo4(klass.dBPA,relstrain,klass.P).*klass.Fm;        %Force magnitude
            FF(relstrain>=1) = 0;
            %For muscle insertion
            unitD = klass.unitD;            %unit direction of force vector, tibia frame
            Fk = unitD.*FF;                  %Force vector, tibia frame
            pB = L(C,:,1);                  %Distance from knee frame to muscle insertion
            Pbri = [-27.5, -125.91, -0.54]/1000;     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever)
            pkbrB = pB-Pbri;                  %vector from bracket to point B, in the knee frame
            thetabrB = atan2(pkbrB(2),pkbrB(1));   %angle between pbrB and x axis
            RkbrZ = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
                   sin(thetabrB) cos(thetabrB) 0; ...
                   0    0   1];
            pbrkB = RkbrZ'*pkbrB';       %Vector in the bracket frame
            % Now calculate angle from x-axis to this vector
%             thetaY = atan2(pbrkB(3), pbrkB(1));  % z vs x (in bracket frame)
%             % Rotation matrix about y-axis (local frame adjustment)
%             Ry = [cos(thetaY)  0  sin(thetaY);
%                   0                1  0;
%                  -sin(thetaY) 0  cos(thetaY)];
%             Rkbr = RkbrZ*Ry;            %Rotate about y-axis in body frame
            Tkbr = RpToTrans(RkbrZ, Pbri');    %Transformation matrix, flexor bracket frame in knee frame
%             %For muscle origin
%             pA = L(C-1,:,1); %Distance from hip frame origin to muscle origin
%             Pbro = [-49.62 -20.79 75.06]/1000; %Distance from hip origin to upper bracket (upper bolt hole)
%             phbrA = pA - Pbro;  %vector from bracket to point A, in the hip frame
%             thetabrA = atan2(phbrA(2),phbrA(1));   %angle between pbrB and x axis
%             RhbrZ = [cos(thetabrA) -sin(thetabrA) 0; ...     %Rotation matrix
%                    sin(thetabrA) cos(thetabrA) 0; ...
%                    0    0   1];
%             pbrhA = RhbrZ'*phbrA';       %Vector in the bracket frame
%             % Now calculate angle from x-axis to this vector
%             thetaYh = atan2(pbrhA(3), pbrhA(1));  % z vs x (in bracket frame)
%             % Rotation matrix about y-axis (local frame adjustment)
%             Ryh = [cos(thetaYh)  0  sin(thetaYh);
%                   0                1  0;
%                  -sin(thetaYh) 0  cos(thetaYh)];
%             Rhbr = RhbrZ*Ryh;            %Rotate about y-axis in body frame
%             Thbr = RpToTrans(RhbrZ, Pbro');    %Transformation matrix, flexor bracket frame in knee frame
            LOC = L;            %new location matrix
            N = size(L,3);
            M = size(L,1);
%             Fh = zeros(N,3);        %Force vector in the hip frame
%             Fbrh = zeros(N,3);       %Force vector represented in the tibial bracket frame  
            pAnew = zeros(N,3);
            Fbrk = zeros(N,3);       %Force vector represented in the tibial bracket frame  
            pBnew = zeros(N,3);
            for ii = 1:N                          %Repeat for each orientation
                    Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),Fk(ii,:)); %Force vector in the tibia frame represented in the lower bracket frame
%                     Fh(ii,:) = norm(RowVecTrans(T(:,:,ii),L(C,:,ii))-L(C-1,:,ii)).*FF(ii,1);    %Force vector in the hip frame
%                     Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:)); %Force vector in the hip frame represented in the lower bracket frame  
            end
%             idx = klass.strain < -0.01;     %when force is too high
%             Fbrk(idx,:) = 1e3*Fbrk(idx,:);                %penalty
            [epsilon, delta, beta] = fortz(klass,Fbrk,X1,X2,X0);  %strain from force divided by tensile stiffness
%             [eA, eB] = fortz_dual(klass, Fbrk, Fbrh, X1, X2, pkbrB, phbrA,X0);
            eB = [epsilon, delta, beta];
            pbrBnew = [norm(pkbrB(1:2)), 0, pkbrB(3)]+eB; %new point B, in the bracket's frame
%             pbrAnew = [norm(phbrA(1:2)), 0, phbrA(3)]+eA; %new point A, in the bracket's frame            
            for ii = 1:N                          %Repeat for each orientation
                for i = 1:M                      %Repeat for all muscle segments
                    if i == C
                        pBnew(ii,:) = RowVecTrans(Tkbr, pbrBnew(ii,:));     %New point B, in the tibia frame
                        LOC(i,:,ii) = pBnew(ii,:);
                    else
                    end

                end
            end
end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ] = fortz(klass,Fbr,X1,X2,X0)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change
    
    % Segment-based muscle length with Xi0 applied via LMT()
    sL = seg(klass, klass.Loc); % current, undeformed path
    if isempty(X0)
        X0=0;
    end
    Lmt = LMT(sL, X0);          % muscle length from geometry
    
    D = klass.dBPA;         %BPA diameter
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
    valid = norms > 1 & all(~isnan(Fbr), 2);
    u_hat_all = zeros(N, 3);
    u_hat_all(valid, :) = Fbr(valid, :) ./ norms(valid);
    
    % Vectorized k_b computation
    K_bracket = diag([X1, X2, X1]);       %project bracket stiffness onto force direction
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    K_rep = repmat(K_bracket, [1, 1, N]);   % [3x3xN]
    k_b = pagemtimes(pagemtimes(u_hat, K_rep), permute(u_hat, [2, 1, 3]));
    k_b = reshape(k_b, [N, 1]);
%     k_eff = 1 ./ (1 ./ k_b );  % Nx1 (If you end up doing one for lower and upper brackets, it would be 1./(1./k_b1+1./k_b2)
    k_eff = k_b;
    
    % Allocate outputs
    e_axial = zeros(N, 1);
    e_bendY = zeros(N, 1);
    e_bendZ = zeros(N, 1); 
    
    % Parallel root solve
   for i = 1:N
    if ~valid(i)
        continue;
    end

    % Per-instance constants
    keff = k_eff(i);
    unit_vec = u_hat_all(i, :);

    % Bracket stiffness disabled (rigid bracket)
    if isinf(X1) && isinf(X2)
        e_axial(i) = 0;
        e_bendY(i) = 0;
        e_bendZ(i) = 0;
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
   end
end

%% Force and BPA length reduction due to bracket deflection (two brackets)
function [eA, eB] = fortz_dual(klass, FbrB, FbrA, X1, X2, pbrB, pbrA,X0)

    % Constants
    rest   = klass.rest; %resting length
    tendon = klass.ten; %tendon length
    fitn   = klass.fitn; %BPA end cap length
    D      = klass.dBPA; %BPA diameter, uninflated
    P      = klass.P / 620;  % Pressure, normalized
    mif    = klass.Fm;  %max. force
    KMAX   = (rest - klass.Kmax)/rest; %max. contraction percent
    sL     = seg(klass, klass.Loc); %length of segments of muscle
    if isempty(X0)
        X0=0;
    end
    Lmt    = LMT(sL, X0);  % Xi0 not yet applied

    % Geometry
    lenB = norm(pbrB);  %length of bracket
    lenA = norm(pbrA);  %length of bracket
    KB = diag([X1, X2, X2]); %Stiffness of lower bracket
    KA = diag([X1, X2, X2]); %Stiffness of upper bracket

    % Normalize forces
    uB = normalize(FbrB);
    uA = normalize(FbrA);
    N  = size(FbrB,1);
    % ✅ Preallocate as Nx3
    eA = zeros(N,3);
    eB = zeros(N,3);

    for i = 1:N
        if any(isnan(FbrA(i,:))) || any(isnan(FbrB(i,:))) || norm(FbrA(i,:)) < 1
            continue;
        end

        % Projections
        u_hatA = uA(i,:)'; 
        u_hatB = uB(i,:)';
        ka = u_hatA' * KA * u_hatA;
        kb = u_hatB' * KB * u_hatB;
        keff = 1 / (1/ka + 1/kb);  % effective spring constant
        
            % Bracket stiffness disabled (rigid bracket)
        if isinf(X1) && isinf(X2)
            eA(i,:) = zeros(1,3);
            eB(i,:) = zeros(1,3);
            continue;
        end

        % Define relfun for fzero
        relfun = @(r_total) ...
            f_festo( ...
                (rest - (Lmt(i) - (tendon + r_total) - 2*fitn)) / rest / KMAX, ...
                P, D ...
            ) * mif - keff * r_total;

        % Initial check
        contraction0 = (rest - (Lmt(i) - tendon - 2*fitn)) / rest;
        relstrain0 = contraction0 / KMAX;

        if relstrain0 >= 1
            r_total = 0;
        else
            try
                r_total = fzero(relfun, [0, 0.2]);
            catch
                r_total = 0;
            end
        end

        r_total = max(r_total, 0);  % safety

        % Split displacements
        rA = r_total * (1/ka) / (1/ka + 1/kb);
        rB = r_total * (1/kb) / (1/ka + 1/kb);

        % Final force value
        contraction = (rest - (Lmt(i) - tendon - 2*fitn - r_total)) / rest;
        relstrain = contraction / KMAX;
        Fmag = f_festo(relstrain, P, D) * mif;

        % Deformation vectors
        eA(i,:) = (KA \ (Fmag * u_hatA))';
        eB(i,:) = (KB \ (Fmag * u_hatB))';
%         e_total = eA + eB;
% 
%         % Store output
%         e_axial(i) = e_total(1);
%         e_bendY(i) = e_total(2);
%         e_bendZ(i) = e_total(3);
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
function mA = Mom(klass, L, unitD)
            T = klass.Tk;
            C = klass.CP;
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L(C, :, i);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
            end
end        
        
        %% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass, L_mt)
rest   = klass.rest;
tendon = klass.ten;
fitting = klass.fitn;

contraction = (rest-(L_mt-tendon-2*fitting))/rest;    %(minus Xi0 is used in LMT function, above)
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
function Mz = Tor(mA, F, maxF, strain)
N = size(F, 1);
Mz = zeros(N, 3);

    for i = 1:N
        if norm(F(i,:)) > maxF || strain(i,:) < -0.03
            Mz(i,:) = NaN;
        else
            Mz(i,:) = cross(mA(i,:), F(i,:));
        end
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

function vhat = normalize(v)
    N = size(v,1);
    norms = vecnorm(v,2,2);
    valid = norms > 1 & all(~isnan(v), 2);
    vhat = zeros(N, 3);
    vhat(valid, :) = v(valid, :) ./ norms(valid);
end


end   