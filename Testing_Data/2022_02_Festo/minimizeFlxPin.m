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
    load KneeFlxPin_10mm_48cm.mat Bifemsh_Pam phiD
    Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeFlxPin_10mm_48cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    kf(1).L(3) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
                  'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
                  'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
                  'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
                  'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
                  'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'L_p',[]);
    clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

if nargout > 1
    % 46cm length
    load KneeFlxPin_10mm_46cm.mat Bifemsh_Pam phiD
    Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeFlxPin_10mm_46cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    A = sortrows([Angle, Torque, InflatedLength, ICRtoMuscle, TorqueHand]);
    kf(1).L(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
                  'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
                  'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
                  'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
                  'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
                  'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[],'F_p',[],'L_p',[]);
    clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
end


%% What to optimize based on selection
% M_opt = zeros(size(kf(1).L(3).Mexp)); %Optimized torque prediction

%% RMSE, fvu, and Max Residual
f = NaN(1,3);
g = NaN(1,3);
a = 1 + ( nargout>1 );
h = cell(1,2);

klaus_temp = kf(1).L(3);   % Grab a template struct
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
    klaus(j) = kf(1).L(4-j);
    % Calculate locations and properties
    L_p{j} = Lok(klaus(j), Xi1, Xi2);
    unitD_p{j} = UD(klaus(j), L_p{j});
    sL_p{j} = seg(klaus(j), L_p{j});
    Lmt_p{j} = LMT(sL_p{j}, Xi0);
    strain_p{j} = Contraction(klaus(j), Lmt_p{j});
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
function LOC = Lok(klass,X1,X2)
            L = klass.Loc;      %Location (wrapping, attachment points)
            C = klass.CP;       %Cross point (moves from one frame to another)
%             T = klass.Tk;       %Transformation matrix
            kmax = klass.Kmax;
            KMAX = (klass.rest-kmax)/klass.rest; %turn it into a percentage 
            relstrain = klass.strain/KMAX;
            FF = festo4(klass.dBPA,relstrain,klass.P).*klass.Fm;        %Force magnitude
            FF(relstrain>1) = 0;
            unitD = klass.unitD;            %unit direction of force vector
            F = unitD.*FF;                  %Force vector
            pB = L(end,:,1);                  %Distance from knee frame to muscle insertion
            Pbr = [-0.043, -.1115, 0];     %vector from knee ICR to flexor insertion bracket (where it starts to cantilever)
%             pkI = [0.2575, -104.25, 0];    %from knee ICR to extensor insertion bracket (where it starts to cantilever)
            pkbrB = pB-Pbr;                  %vector from bracket to point B, in the knee frame
            thetabrB = wrapToPi(acos(dot([pkbrB(1) pkbrB(2) 0],[1,0,0])/(norm([pkbrB(1) pkbrB(2) 0]))));   %angle between pbrB and x axis
            Rkbr = [cos(thetabrB) -sin(thetabrB) 0; ...     %Rotation matrix
                   sin(thetabrB) cos(thetabrB) 0; ...
                   0    0   1];
            Tkbr = RpToTrans(Rkbr, Pbr');    %Transformation matrix, flexor bracket frame in knee frame                
            LOC = L;            %new location matrix
            N = size(L,3);
            M = size(L,1);
            Fbrk = zeros(N,3);       %Force vector represented in the bracket frame  
            pBnew = zeros(N,3);
            for ii = 1:N                          %Repeat for each orientation
                    Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),F(ii,:)); %Force vector in the tibia frame represented in the bracket frame
                    
            end
%             idx = klass.strain < -0.01;     %when force is too high
%             Fbrk(idx,:) = 1e3*Fbrk(idx,:);                %penalty
            [epsilon, delta, beta] = fortz(klass,Fbrk,X1,X2);  %strain from force divided by tensile stiffness
            pbrBnew = [norm([pkbrB(1) pkbrB(2)])+epsilon, delta, pkbrB(3)+beta]; %new point B, in the bracket frame
                        
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
function Lmt = LMT(sL_p, Xi0)
            Lmt1 = sL_p;
            Lmt = Lmt1 - Xi0;
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
                %mA(i, :) = cross(pointB, unitD_p(i, :));
            end
end        
        
        %% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass, Lmt_p)
            rest = klass.rest;
            tendon = klass.ten;
            fitting = klass.fitn;
            
            contraction = (rest-(Lmt_p-tendon-2*fitting))/rest;    %(minus Xi0 is used in LMT function, above)
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

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ] = fortz(klass,Fbr,X1,X2)
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
    valid = norms > 1 & all(~isnan(Fbr), 2);
    u_hat_all = zeros(N, 3);
    u_hat_all(valid, :) = Fbr(valid, :) ./ norms(valid);
    
    % Vectorized k_b computation
    K_bracket = diag([X1, X2, X2]);       %project bracket stiffness onto force direction
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
    if isinf(X1) || isinf(X2)
        e_axial(i) = 0;
        e_bendY(i) = 0;
        e_bendZ(i) = 0;
        continue;
    end

    % Solve for r using fzero
    relfun = @(r) ...
        f_festo( ...
            (rest - (mL(i) - (tendon + r) - 2*fitting)) / rest / KMAX, ...
            P, D ...
        ) * mif - keff * r;

    try
        r = fzero(relfun, [0, 0.2]);
    catch
        r = 0;
    end

    if isnan(r)
        e_axial(i) = 0;
        e_bendY(i) = 0;
        e_bendZ(i) = 0;
        continue;
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