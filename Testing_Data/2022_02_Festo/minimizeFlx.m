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
%     clear Bifemsh_Pam3 Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A
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

Funit = cell(1,a);
strain_Xi0 = cell(1,a);
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
    %Borrowed from minimizeExt3, for reference
%     bpa_i = klass;
%     kspr = Spr(bpa_i);          %tendon spring rate
%     Funit = computeForceVector(bpa_i);  %Force unit direction, calculate in the hip frame 
%     [strain_Xi3, delta_L] = Contraction(bpa_i, [], Xi0, [], Xi3); %Calculate contraction and loss of length due to constant length offset Xi0 and wrapping factor Xi3
%     [L_p, gemma] = Lok(bpa_i, Xi1, Xi2,kspr, Funit, strain_Xi3, delta_L+Xi0); %Bracket deformation and new geometry
%     sL_p = seg(bpa_i, L_p); %segment lengths
%     Lmt_p = LMT(sL_p, Xi0); 
%     [strain_f, ~] = Contraction(bpa_i, Lmt_p, [], gemma, Xi3); %Simulated strain, includes loss of usable length due to wrapping
%     unitD_p = UD(bpa_i, L_p);           %unit direction, predicted with updated path
%     F_p = Force(bpa_i, unitD_p, strain_f); %Muscle force, new prediction
%     mA_p = Mom(bpa_i, L_p, unitD_p); %Moment arm vector
%     [strain_p, ~] = Contraction(bpa_i, Lmt_p, [], gemma, []);
%     M_p = Tor(mA_p, F_p, bpa_i.Fm, strain_p); %New torque prediction, using new moment arm, new Force vector, Fmax, and strain_p
    %updated code
    kspr{j} = Spr(klaus(j)); %Calculate the spring rate
    Funit{j} = computeForceVector(klaus(j));  %Calculate force vector in the hip frame
    strain_Xi0{j} = Contraction(klaus(j), [],[],Xi0); %*Calculate contraction from constant length offset
    [L_p{j}, gemma{j}] = Lok(klaus(j), Xi1, Xi2,kspr{j},Funit{j}, strain_Xi0{j}, Xi0); %Calculate deformed bracket geometry and stretched tendon with stiffnesses and length offset.
    unitD_p{j} = UD(klaus(j), L_p{j}); %Calculate force unit vector with deformed geometry
    sL_p{j} = seg(klaus(j), L_p{j}); %Calculate segment lengths with deformed geometry
    Lmt_p{j} = LMT(sL_p{j}, Xi0);   %Calculate musculotendon length with deformed geometry and constant length offset
    strain_p{j} = Contraction(klaus(j), Lmt_p{j},gemma{j},[]); %*Calculate contraction with bracket deformation, tendon stretch, and constant length offset
    F_p{j} = Force(klaus(j), unitD_p{j}, strain_p{j}); %Calculate force vector in body frame
    mA_p{j} = Mom(klaus(j), L_p{j}, unitD_p{j}); %Calculate moment arm
    M_p{j} = Tor(klaus(j), mA_p{j}, F_p{j}, strain_p{j}); %Torque
    %Original code, for reference
%     kspr{j} = Spr(klaus(j));
%     [L_p{j}, gemma{j}] = Lok(klaus(j), Xi1, Xi2,kspr{j});
%     unitD_p{j} = UD(klaus(j), L_p{j});
%     sL_p{j} = seg(klaus(j), L_p{j});
%     Lmt_p{j} = LMT(sL_p{j}, Xi0);
%     strain_p{j} = Contraction(klaus(j), Lmt_p{j},gemma{j});
%     F_p{j} = Force(klaus(j), unitD_p{j}, strain_p{j});
%     mA_p{j} = Mom(klaus(j), L_p{j}, unitD_p{j});
%     M_p{j} = Tor(mA_p{j}, F_p{j}, klaus(j).Fm, strain_p{j});
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

%% -------------Force unit direction ---------------
function F_unit = computeForceVector(klass)
%Calculate the force unit direction from muscle origin (hip frame) to the next
%real point. This takes into account if there are any additional via
%points between muscle origin and muscle insertion. It also takes into
%account if a homogenous transform+ation matrix needs to be used to
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
function contraction = Contraction(klass,Lmt_p,gema,X0)
            rest = klass.rest;      %resting length
            tendon = klass.ten;     %artificial tendon length
            fitting = klass.fitn;   %fitting length
            
            if isempty(Lmt_p)
                Lmt = klass.Lmt;
            else
                Lmt = Lmt_p;    %minus Xi0 is used in LMT function
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
function [LOC, gema] = Lok(klass,X1,X2,kSpr,Funit,strain_predef,X0)
% Inputs:
%   bpa class info
%   X1, X2 stiffness
%   kSpr, tendon stiffness
%   Funit, force unit direction in the hip frame
%   strain_predef – N×1 strain vector (e.g., from Xi0 offset effect)
%   X0, constant length offset
    L = klass.Loc;          %Location of wrapping and attachment points
    rest = klass.rest;      %resting length
    Fm = klass.Fm;          %maximum isometric force
    P = klass.P;            %BPA pressure
    D = klass.dBPA;         %BPA diameter
    KMAX = (rest - klass.Kmax)/rest; %maximum contracted length (meters)
    N = size(L,3);
    
    % Compute Force
    relstrain = strain_predef / KMAX;  %Relative strain
    FF = zeros(size(strain_predef));
    idx = strain_predef > -0.02 & relstrain < 1;
    FF (idx)= festo4(D, relstrain(idx), P) * Fm; %Force magnitude
    F = Funit * FF;  % N×3, already in hip frame
    
    pA = L(1,:,1);                                  %Distance from hip origin to muscle insertion
            switch klass.dBPA
                case 20
%                   Pbr = [-0.8100  -20.222   31.66]/1000;       %from hip origin to bracket bolt closest to the origin of the Bifemsh_Pam
                    Pbr = [9.48  -36.2   30.76]/1000;       %from hip origin to bracket bolt pattern centroid
                case 10
                    Pbr = [-19 22 27.6]/1000;       %from hip origin centroid of bracket cantilever 
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
%                         Fh(ii,:) = -RowVecTrans(T(:,:,ii),F(ii,:));               %Force vector represented in the hip frame, opposite direction
%                         Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:));            %Force vector in the hip frame represented in the bracket frame
                    Fbrh(ii,:) = RowVecTrans(Thbr\eye(4),F(ii,:));            %Force vector in the hip frame represented in the bracket frame
            end
            [epsilon, delta, beta, gema] = fortz(klass,F,X1,X2,kSpr,X0);  %strain from force divided by tensile stiffness
            deflection = [epsilon, delta, beta];    %bracket movement
            pbrAnew = [norm(pbrhA),0,0]+deflection; %New point A, represented in the bracket frame
            
            LOC = L;
            for ii = 1:N                          %Repeat for each orientation 
                pAnew = RowVecTrans(Thbr, pbrAnew(ii,:)); %New point A in the hip frame
                LOC(1,:,ii) = pAnew;      %Update location matrix
            end

end

%% Force and length reduction due to tendon
function [e_axial, e_bendY, e_bendZ, e_cable] = fortz(klass,Fbr,X1,X2,kSpr,X0)
% e_axial, bracket axial elongation
% e_bendY, bracket bending displacement y - direction
% e_bendZ, bracket bending displacement z - direction
% e_cable, tendon cable stretch
% total length change
    
    % Initialize outputs
    [e_axial, e_bendY, e_bendZ, e_cable] = deal(zeros(N,1));
    
    if isinf(X1) && isinf(X2) && isinf(kSpr)
        return;
    end
    
    D = klass.dBPA;         %BPA diameter
    if isempty(X0)
        X0 = 0;
    end
    rest = klass.rest;      %resting length
    tendon = klass.ten;     %tendon length
    fitn = klass.fitn;    %fitting length
    mL = klass.Lmt - X0 - tendon - 2*fitn;   %musculotendon length
    mif = klass.Fm;         %maximum force
    kmax = klass.Kmax;      %maximum contracted length
    KMAX = (rest-kmax)/rest; %turn it into a percentage
    P = klass.P;            %pressure
    
    N = size(Fbr,1);
    % Normalize force vectors safely
    norms = vecnorm(Fbr, 2, 2);
    valid = norms > 1e-3 & all(~isnan(Fbr), 2);
    u_hat_all = normalize(Fbr);
    
    % Vectorized k_b computation
    K_bracket = diag([X1, X2, X2]);       %project bracket stiffness onto force direction
    u_hat = permute(u_hat_all, [3, 2, 1]);  % [1x3xN]
    K_rep = repmat(K_bracket, [1, 1, N]);   % [3x3xN]
    k_b = pagemtimes(pagemtimes(u_hat, K_rep), permute(u_hat, [2, 1, 3]));
    k_b = reshape(k_b, [N, 1]);
    k_eff = 1 ./ (1 ./ k_b + 1 ./ kSpr);  % Nx1

%     ub = rest*KMAX/2;         %upper bound of fzero
    
    % Parallel root solve
    for i = 1:N
        if ~valid(i)
            continue;
        end
        
        % Per-instance constants
        keff = k_eff(i);
        unit_vec = u_hat_all(i, :);
        Lm = mL(i);
        
%         r = 0.0001;     %initial guess

        % Solve for r (deflection) using fzero
        % First, compute contraction from Xi0 if deflection is zero 
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
                r = fzero(relfun, [0, .1]);
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
%             r_bracket = unit_vec * e_bkt;
            e_cable(i) = F_mag/kSpr;
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
        
        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
function F = Force(klass, unitD_p, strain)
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
           scalarForce(scalarForce < 0) = 0;            
            
            F = scalarForce.*unitD_p;

end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
function Mz = Tor(klass, mA_p, F_p, strain_p)  
            Mz = zeros(size(F_p));
           
            switch klass.dBPA
                case 20
                    ss = -.03;       %maximum allowable strain
                case 10
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
        
        springrate = mult*Aeff*E/L;        
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
    valid = norms > 1e-3 & all(~isnan(v), 2);
    vhat = zeros(N, 3);
    vhat(valid, :) = v(valid, :) ./ norms(valid);
end


end