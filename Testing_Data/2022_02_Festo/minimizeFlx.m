%% Optimize predicted torque for extensors.
function [f, varargout] = minimizeFlx(Xi0,Xi1,Xi2)

%% 46cm length
%kf = knee flexor, kf(1) = pinned joint, kf(2) = biomimetic;
%ke = knee extensor, same as above
%kf.L := lengths = [42 46 48] cm
%example: kf(1).L(2).Mz z-axis torque for pinned knee, flexor, 46cm length
%'exp' suffix means experimental
%'_h' suffix means hybrid
%'_p' suffix means prime, as in the new prediction values

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
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

%% 48 cm length
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
                  'Lmt_p',[],'mA_p',[],'M_p',[]);
clear Bifemsh_Pam phiD Ma G Angle Torque InflatedLength ICRtoMuscle TorqueHand A

%% What to optimize based on selection
M_opt = zeros(length(kf(1).L(2).Mexp)); %Optimized torque prediction

%% RMSE, fvu, and Max Residual
h = cell(1,2);
for j = 1:2
    klaus(j) = kf(1).L(j+1);
    bpa(j) = klaus(j);
    bpa(j).Lmt_p = LMT(klaus(j));
    bpa(j).mA_p = Mom(klaus(j));    
    bpa(j).M_p = Tor(klaus(j));
    h{j} = SSE(klaus(j));
end

f = h{1};
g = h{2};
varargout{1} = g;
varargout{2} = bpa;

%% Nested functions, modified from MonoPamExplicit
        %% ------------- Location  ------------------------
        function LOC = Lok(klass)
            L = klass.Loc;      %Location (wrapping, attachment points)
            C = klass.CP;       %Cross point (moves from one frame to another)
%             T = klass.Tk;       %Transformation matrix
            F = klass.F;        %Force vector
            unitD = klass.unitD;
            pkI = [-0.043, -.1115, 0];     %from knee ICR to flexor insertion bracket (where it starts to cantilever)
%             pkI = [0.2575, -104.25, 0];    %from knee ICR to extensor insertion bracket (where it starts to cantilever)
            RkI = [1, 0, 0;                 %Rotation matrix (no rotation)
                   0, 1, 0;
                   0, 0, 1];
            TkI = RpToTrans(RkI, pkI');    %Transformation matrix from knee to flexor bracket frame
%             TkI = RpToTrans(RkI, pkeI');    %Transformation matrix from knee to extensor bracket frame
                    
            LOC = zeros(size(L));
            pB = zeros(size(L,3),3);
            pBI = zeros(size(L,3),3);
            uDkI = zeros(size(L,3),3);
            thetaBI = zeros(size(L,3),1);
            Rot = zeros(3,3,size(L,3));
            yprime = zeros(size(L,3),3);
            epsilon = zeros(size(L,3),3);
            delta = zeros(size(L,3),3);
            pBInew = zeros(size(L,3),3);
            pBnew = zeros(size(L,3),3);
            for ii = 1:size(L, 3)                          %Repeat for each orientation
                for i = 1:size(L, 1)                      %Repeat for all muscle segments
                    LOC(i,:,ii) = L(i, :,ii);
                    if i == C
                        pB(ii,:) = L(C, :,ii);
                        pBI(ii,:) = RowVecTrans(TkI,pB(ii,:));
                        uDkI(ii,:) = RowVecTrans(TkI,unitD(ii,:));
                        thetaBI(ii) = acos(dot(pBI(ii,:),[1,0,0])/(norm(pBI(ii,:))));   %angle between pBI and x axis
                        Rot(:,:,ii) = [cos(thetaBI(ii)) -sin(thetaBI(ii)) 0; ...
                                       sin(thetaBI(ii)) cos(thetaBI(ii)) 0; ...
                                       0    0   1];
                        yprime(ii,:) = Rot(:,:,ii)*[0 1 0]';
                        epsilon(ii,:) = ((dot(F(ii,:),pBI(ii,:))/norm(pBI(ii,:))^2)*pBI(ii,:))./Xi1;  %strain from tensile stiffness
                        delta(ii,:) = ((dot(F(ii,:),yprime(ii,:))/norm(yprime(ii,:))^2)*yprime(ii,:))./Xi2;    %deflection bending stiffness
                        pBInew(ii,:) = pBI(ii,:)+delta(ii,:)+epsilon(ii,:);
                        pBnew(ii,:) = RowVecTrans(TkI\eye(4), pBInew(ii,:));
                        LOC(i,:,ii) = pBnew(ii,:);
                    end

                end
            end
        end

%% ------------- Segment Lengths ------------------------
        function SL = seg(klass)
            L = Lok(klass);
            C = klass.CP;
            T = klass.Tk;
            SL = zeros(size(T, 3), size(L, 1) - 1);
            
            for ii = 1:size(T, 3)                          %Repeat for each orientation
                for i = 1:size(L, 1)-1                      %Repeat for all muscle segments
                    pointA = L(i, :,ii);
                    pointB = L(i+1, :,ii);
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
            L = Lok(klass);
            T = klass.Tk;
            Lmt = zeros(size(T, 3), 1);
            segLengths = seg(klass);
            
            for ii = 1:size(Lmt, 1)                          %Repeat for each orientation
                for i = 1:size(L, 1, 1)-1                      %Repeat for all muscle segments
                    Lmt(ii, 1) = Lmt(ii, 1) + segLengths(ii, i);
                end
            end
            Lmt = Lmt - Xi0;
        end            

        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
        function unitD = UD(klass)
            L = Lok(klass);
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
        function mA = Mom(klass)
            T = klass.Tk;
            L = Lok(klass);
            C = klass.CP;
            unitD = UD(klass);
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L(C, :, i);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
                %mA(i, :) = cross(pointB, unitD(i, :));
            end
        end        
        
        %% -------------- Contraction of the PAM --------------------------
        function contraction = Contraction(klass)
            Lmt = LMT(klass);      %musculotendon length, this calls optimvar Xi0
            rest = klass.rest;
            tendon = klass.ten;
            fitting = klass.fitn;
            
            contraction = (rest-(Lmt-tendon-2*fitting))/rest;    %(minus Xi0 is above)
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
            unitD = UD(klass);
            contract = Contraction(klass);
            rest = klass.rest;
            pres = klass.P;
            kmax = klass.Kmax;  
            KMAX = (rest-kmax)/rest; %turn it into a percentage 
            maxF = klass.Fm;
            
           rel = contract./KMAX;                    %relative strain        
           relPres = pres/620;                      %relative pressure
           
           if dia == 10
                load FestoLookup.mat f_10
                Fn = f_10(rel,relPres);
           elseif dia == 20
               load FestoLookup.mat f20
               Fn = f20(rel,relPres);
           elseif dia == 40
               load FestoLookup.mat f40
               Fn = f40(rel,relPres);
           end 
           scalarForce = Fn.*maxF;

            for i = 1:size(unitD, 1)
                if scalarForce(i) < 0
                    scalarForce(i) = 0;
                end
                if scalarForce(i) > maxF
                    scalarForce(i) = NaN;
                end
            end
            
            
            F = scalarForce.*unitD;

        end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
        function Mz = Tor(klass)
            mA = Mom(klass);
            F = Force(klass);
            Mz = cross(mA,F,2);
            
            for i = 1:size(mA, 1)
                Mz(i, :) = cross(mA(i, :), F(i, :));
            end

        end    

%% Subfunctions
function t = SSE(klass)
     Mpredict1 = Tor(klass);
     Mpredict2 = griddedInterpolant(klass.Ak,Mpredict1(:,3));
     M_opt = Mpredict2(klass.Aexp);
     [RMSE, fvu, maxResid] = Go_OfF(klass.Mexp,M_opt);
     t = [RMSE, fvu, maxResid];
end




end   