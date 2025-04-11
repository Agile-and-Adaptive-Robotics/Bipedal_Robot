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
    load KneeFlx_10mm_42cm.mat Bifemsh_Pam phiD
    Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeFlx_10mm_42cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    A = sortrows([Angle', Torque', InflatedLength', ICRtoMuscle', TorqueHand']);
    kf(2).L(1) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
                  'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
                  'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
                  'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
                  'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
                  'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
    clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand Angle

if nargout > 1
    % 42cm length, 20mm diameter.
    load KneeFlx_20mm_42cm.mat Bifemsh_Pam3 phiD
    Bifemsh_Pam = Bifemsh_Pam3;
    Ma = Bifemsh_Pam.MomentArm;                 %Calculated moment arm
    G = (Ma(:,1).^2+Ma(:,2).^2).^(1/2);         %Moment arm for z-axis torque
    load Plot_KneeFlx_20mm_42cm.mat Angle Torque InflatedLength ICRtoMuscle TorqueHand
    A = sortrows([Angle', Torque, InflatedLength', ICRtoMuscle', TorqueHand']);
    kf(2).L(2) = struct('Ak',phiD,'Loc',Bifemsh_Pam.Location,'CP',Bifemsh_Pam.Cross,'dBPA',Bifemsh_Pam.Diameter, ...
                  'Tk',Bifemsh_Pam.TransformationMat,'rest',Bifemsh_Pam.RestingL,'Kmax',Bifemsh_Pam.Kmax,...
                  'fitn',Bifemsh_Pam.FittingLength,'ten',Bifemsh_Pam.TendonL,'P',Bifemsh_Pam.Pressure, ...
                  'Lmt',Bifemsh_Pam.MuscleLength,'strain',Bifemsh_Pam.Contraction, 'unitD',Bifemsh_Pam.UnitDirection, ...
                  'mA',G,'Fm',Bifemsh_Pam.Fmax,'F',Bifemsh_Pam.Force, 'seg',Bifemsh_Pam.SegmentLengths, ...
                  'M',Bifemsh_Pam.Torque(:,3),'Aexp',A(:,1),'Mexp',A(:,2),...
                  'A_h',A(:,1),'Lm_h',A(:,3),'mA_h',A(:,4),'M_h',A(:,5),...
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
    clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand Angle
end


%% What to optimize based on selection
M_opt = zeros(size(kf(2).L(1).Mexp)); %Optimized torque prediction

%% RMSE, fvu, and Max Residual
if nargout>1
    a = 2;
else
    a = 1;
end
h = cell(1,2);

for j = 1:a
    klaus(j) = kf(2).L(j);
    bpa(j) = klaus(j);
    kspr = Spr(klaus(j));
    [L_p, gemma] = Lok(klaus(j));
    unitD_p = UD(klaus(j));
    sL_p = seg(klaus(j));
    Lmt_p = LMT(klaus(j));
    bpa(j).Lmt_p = Lmt_p;
    mA_p = Mom(klaus(j));
    strain_p = Contraction(klaus(j));
    F_p = Force(klaus(j));
    bpa(j).mA_p = mA_p;
    M_p = Tor;
    bpa(j).M_p = M_p;
    h{j} = SSE(klaus(j));
end

f = h{1};
g = h{2};
varargout{1} = g;
varargout{2} = bpa;

%% Nested functions, mostly modified from MonoPamExplicit

%% ------------- Location  ------------------------
function [LOC, gamma] = Lok(klass)
            L = klass.Loc;      %Location (wrapping, attachment points)
            C = klass.CP;       %Cross point (moves from one frame to another)
            T = klass.Tk;       %Transformation matrix, tibia frame represented in the hip frame
            switch klass.dBPA
                case 10
                    D = 10;
                case 20
                    D = 20;
            end
            FF = festo4(D,klass.strain/(klass.rest-klass.Kmax),klass.P).*klass.Fm;       %Force magnitude
            unitD = klass.unitD;                                                %Unit direction of force, tibia frame
            F = unitD.*FF;                                            %Force vector, tibia frame
            pA = L(1,:,1);                                  %Distance from hip origin to muscle insertion
            Pbr = [-0.8100  -20.2650   32.2100]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
            pbrA = pA-Pbr;                                  %vector from bracket to point A (in the hip frame)
            thetabrA = norm(wrapToPi(acos(dot(pbrA,[1,0,0])/(norm(pbrA)))));   %angle between pbrA and x axis
            Rhbr = [cos(thetabrA) -sin(thetabrA) 0; ...     %Rotation matrix
                   sin(thetabrA) cos(thetabrA) 0; ...
                   0    0   1];
            Thbr = RpToTrans(Rhbr, Pbr');    %Transformation matrix, represent bracket frame in hip frame              
            LOC = zeros(size(L));
            Fh = zeros(size(L,3),3);
            Fbrh = zeros(size(L,3),3);
            epsilon = zeros(size(L,3),1);
            delta = zeros(size(L,3),1);
            gamma = zeros(size(L,3),1);
            pbrAnew = zeros(size(L,3),3);   %New point A, represented in bracket frame
            pAnew = zeros(size(L,3),3);     %New point A, in the hip frame
            for ii = 1:size(L, 3)                          %Repeat for each orientation
                for i = 1:size(L, 1)                      %Repeat for all muscle segments
                    if i == C-1
                        Fh(ii,:) = -RowVecTrans(T(:,:,ii),F(ii,:));               %Force vector represented in the hip frame
                        Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:));   %Force vector in the hip frame represented in the bracket frame
                        [epsilon(ii), delta(ii), gamma(ii)] = fortz(klass,Fbrh(ii,1:2));  %strain from force divided by tensile stiffness
                        pbrAnew(ii,:) = [norm([pbrA(1) pbrA(2)])+epsilon(ii,:), delta(ii,:), pbrA(3)]; %New point A, represented in the bracket frame
                        pAnew(ii,:) = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
                        LOC(i,:,ii) = pAnew(ii,:);      %Update location matrix
                    else
                        LOC(i,:,ii) = L(i, :,ii);
                    end

                end
            end
end

%% ------------- Segment Lengths ------------------------
function SL = seg(klass)
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
function Lmt = LMT(klass)
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
function unitD = UD(klass)
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
function mA = Mom(klass)
            T = klass.Tk;
            C = klass.CP;
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L_p(C, :, i);
                mA(i, :) = pointB - unitD_p(i, :)*dot(unitD_p(i, :), pointB);
            end
end        
        
        %% -------------- Contraction of the PAM --------------------------
function contraction = Contraction(klass)
            rest = klass.rest;
            tendon = klass.ten;
            fitting = klass.fitn;
            
            contraction = (rest-(Lmt_p-tendon-2*fitting)-gemma)/rest;    %(minus Xi0 is used in LMT function, above)
end


        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
function F = Force(klass)
        %Inputs:
        %Lmt == muscle-tendon length, scalar
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %pres == measured pressure
        %kmax == maximum contraction length
        %Outputs:
        %F == Force, N           
            dia = klass.dBPA;
            rest = klass.rest;
            pres = klass.P;
            kmax = klass.Kmax;  
            KMAX = (rest-kmax)/rest; %turn it into a percentage 
            maxF = klass.Fm;
            
           rel = strain_p./KMAX;                    %relative strain        
           relPres = pres/620;                      %relative pressure
           
           Fn = festo4(dia,rel,relPres);

           scalarForce = Fn.*maxF;
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
function Mz = Tor(~)
            Mz = zeros(size(F_p));
            
            for i = 1:size(F_p, 1)
                Mz(i, :) = cross(mA_p(i, :), F_p(i, :));
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
        Aeff = 1.555*10^-6;%Effective area for 19-strand cable
        E = 193*10^9;       %Young's Modulus
        L = klass.ten;      %tendon length
        
        springrate = mult*Aeff*E/L;
        
end

%% Force and length reduction due to tendon
function [e_axial, e_bending, e_cable] = fortz(klass,Fbrh)
% e_axial, bracket axial elongation
% e_bending, bracket bending displacement
% e_cable, tendon cable stretch

    D = klass.dBPA;         %BPA diameter
    mL = klass.Lmt;         %musculotendon length
    rest = klass.rest;      %resting length
    tendon = klass.ten;     %tendon length
    fitting = klass.fitn;    %fitting length
    mif = klass.Fm;         %maximum force
    kmax = klass.Kmax;      %maximum contracted length
    KMAX = (rest-kmax)/rest; %turn it into a percentage 
    
    % Ensure Fbrh_vec is a column vector
    if isrow(Fbrh)
        Fbrh = Fbrh';
    end
    
    u_hat = Fbrh./norm(Fbrh);      %Normalized force direction vector
    k_b = u_hat' * diag([Xi1, Xi2]) *u_hat;       %project bracket stiffness onto cable direction
    k_eff = 1 / (1/k_b + 1/kspr);   %effective stiffness
    function Fbal = myfunc(r)
        contraction = (rest-(mL-tendon-2*fitting)-r)/rest;
        rel = contraction/KMAX;
        fM = festo4(D,rel,klass.P)*mif;     % Muscle force
        fT = k_eff*r;                       % tendon force
        Fbal = fM-fT;

    end 

    x0 = 1e-4; %Initial guess    
    options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-6);    
    r2 = fsolve(@myfunc,x0,options);
    
    contract = (rest-(mL-tendon-2*fitting)-r2)/rest;
    relstrain = contract/KMAX;
    F_mag = festo4(D,relstrain,klass.P)*mif;     % Muscle force at equilibrium
%     Fbrh_final = F_mag*u_hat;
    
    %bracket displacement vector
    e_bkt = k_b \ (F_mag * u_hat);
    e_axial = e_bkt(1);
    e_bending = e_bkt(2);
    
    %Project bracket displacement into cable direction
    r_bracket = u_hat' * e_bkt;
    
    % Cable elongation
    e_cable = r2 - r_bracket;

end

%% Subfunctions
function t = SSE(klass)
     Mpredict1 = Tor(klass);
%      Mpredict1 = Mom(klass);
%      G_p = (Mpredict1(:,1).^2+Mpredict1(:,2).^2).^(1/2);
     Mpredict2 = griddedInterpolant(klass.Ak,Mpredict1(:,3));
%      Mpredict2 = griddedInterpolant(klass.Ak,G_p);
     M_opt = Mpredict2(klass.Aexp);
     [RMSE, fvu, maxResid] = Go_OfF(klass.Mexp,M_opt);
     t = [RMSE, fvu, maxResid];
end

end