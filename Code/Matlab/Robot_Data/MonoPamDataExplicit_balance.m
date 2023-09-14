% Pam Data
% Author: Connor Morrow
% Date: 6/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPamDataExplicit_balance < handle
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which row corresponds with a location where the muscle crosses into a new reference frame
        Diameter                    %Diameter of the BPA
        TransformationMat           %Contains a transformation matrix to change the 
        RestingL                    %Resting Length of the muscle
        Kmax                        %Length of BPA at maximum contraction
        FittingLength               %Length of each end cap (center of hole to bottom port)
        TendonL                     %Length of tendon, if any
        Pressure                    %Pressure of BPA
        SEstiff                     %Series elastic stiffness 
    end
    
    %Dependent properties are those that are calculated by the explicit
    %properties. Matlab will not calculate these until it is queried in the
    %main script
    properties (Dependent)   
        SegmentLengths
        LongestSegment
        MuscleLength
        Contraction
        LengthCheck
        UnitDirection
        MomentArm
        Fmax
        Force
        Fbal
        Torque
    end
    
   
    methods
        %% ------------- Muscle Data Constructor -----------------
        %Constructor Function. By calling 'MuscleData' and entering the
        %muscle information, we construct an object for that muscle.
        function PD = MonoPamDataExplicit(name, location, cross, diameter, t, rest, kmax, tendon, fitn, pres, stiff)
            if nargin == 10
                PD.Name = name;
                PD.Location = location;
                PD.Cross = cross;
                PD.Diameter = diameter;
                PD.TransformationMat = t;
                PD.RestingL = rest;
                PD.Kmax = kmax;
                PD.TendonL = tendon;
                PD.FittingLength = fitn;
                PD.Pressure = pres;
                PD.SEstiff = stiff;
            elseif nargin == 7
                PD.Name = name;
                PD.Location = location;
                PD.Cross = cross;
                PD.Diameter = diameter;
                PD.TransformationMat = t;
                PD.RestingL = rest;
                PD.Kmax = kmax;
                PD.TendonL = 0;
                PD.FittingLength = 0.0254;
                PD.Pressure = 620;
                PD.SEstiff = 1.5834e+06;      % Stiffness of 1/8" steel wire rope is 1.5834e+06
            else
                fprintf('Invalid number of arguments\n')
            end
        end
        
        
        %% ------------- Segment Lengths ------------------------
        function segLengths = get.SegmentLengths(obj)
            L = obj.Location;
            C = obj.Cross;
            T = obj.TransformationMat;
            segLengths = zeros(size(T, 3), size(L, 1) - 1);
            
          
            for ii = 1:size(T, 3)                          %Repeat for each orientation
                for i = 1:size(L, 1)-1                      %Repeat for all muscle segments
                    pointA = L(i, :,ii);
                    pointB = L(i+1, :,ii);
                    if i+1 == C
                        pointB = RowVecTrans(T(:, :, ii), pointB);
                    end
                    segLengths(ii, i) = norm(pointA - pointB);
                end
            end
        end
        
        %% -------------- Longest Segment Calculation ----------------
        function longestSeg = get.LongestSegment(obj)
            L = obj.Location;
            segLengths = obj.SegmentLengths;
            
            % Calculate which muscle segment is the longest on average.
            % This will be where the Pam resides.
            avgSegL = zeros(size(L, 1, 1) - 1);
            for i = 1:size(segLengths, 2)
                avgSegL(i) = mean(segLengths(:, i));
            end
            
            longestSegPointer = 1;
            if size(avgSegL, 1) > 1
                for i = 1:size(avgSegL, 1) - 1
                    if avgSegL(i + 1) > avgSegL(i)
                        longestSegPointer = i + 1;                        
                    end
                end
            end
            longestSeg = segLengths(:, longestSegPointer); 
        end
        
        %% ------------- Muscle Length ------------------------
        %Function that calculates the muscle length, based
        function mL = get.MuscleLength(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            mL = zeros(size(T, 3), 1);
            segLengths = obj.SegmentLengths;
            
            for ii = 1:size(mL, 1)                          %Repeat for each orientation
                for i = 1:size(L, 1, 1)-1                      %Repeat for all muscle segments
                    mL(ii, 1) = mL(ii, 1) + segLengths(ii, i);
                end
            end
        end            

        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
        function unitD = get.UnitDirection(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            C = obj.Cross;
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
        function mA = get.MomentArm(obj)
            T = obj.TransformationMat;
            L = obj.Location;
            C = obj.Cross;
            unitD = obj.UnitDirection;
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L(C, :, i);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
                %mA(i, :) = cross(pointB, unitD(i, :));
            end
        end
        
        %% -------------- Resting PAM Length --------------------------

        %% -------------- Contraction of the PAM --------------------------
        function contraction = get.Contraction(obj)
            mL = obj.MuscleLength;
            rest = obj.RestingL;
            tendon = obj.TendonL;
            fitting = obj.FittingLength;
            
            contraction = (rest-(mL-tendon-2.*fitting))./rest;
%             contraction = zeros(length(mL), 1);
%             for i = 1:length(mL)
%                 contraction(i) = (rest-(mL(i,1)-tendon-2*fitting))/rest;
% %                   contraction(i) = ((rest+tendon+2*fitting)-mL(i,1))/(rest+tendon+2*fitting);
%             end
        end
        
        %% -------------- Length Check --------------------------
        function lengthCheck = get.LengthCheck(obj)
            contraction = obj.Contraction;
            maxContractPercent = 0.25;          %Contracting to 75% of length
            minContractPercent = -0.1;          %Elongating to 110% of length
            restingPamLength = obj.RestingL;
            
            if restingPamLength < 0
                lengthCheck = 'Unusable';
            else
                if max(contraction) <= maxContractPercent
                    if min(contraction) >= minContractPercent
                        lengthCheck = 'Usable';
                    else
                        lengthCheck = 'Unusable';
                    end
                else
                    lengthCheck = 'Unusable';
                end
            end
        end

        %% -------------- Maximum Force --------------------------
        %Calculate the direction of the forced applied by the muscle
        function maxF = get.Fmax(obj)
        %Inputs:
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %Outputs:
        %maxF == Maximum Force, N, produced by BPA at 0% contraction and
        %           620 kPa
        
            dia = obj.Diameter;
            rest = obj.RestingL;

            if dia == 10    
                maxF = maxBPAforce(rest,620);
            elseif dia ==20
                maxF = 1500;
            elseif dia ==40
                maxF = 6000;
            else
                disp('Wrong size diameter BPA')
            end
        end

        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
        function F = get.Force(obj)
        %Inputs:
        %Lmt == muscle-tendon length, scalar
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %pres == measured pressure
        %kmax == maximum contraction length
        %Outputs:
        %F == Force, N           
            dia = obj.Diameter;
            unitD = obj.UnitDirection;
            contract = obj.Contraction;
            mL = obj.MuscleLength;
            rest = obj.RestingL;
            fitting = obj.FittingLength;
            pres = obj.Pressure;
            kmax = obj.Kmax;  
            KMAX = (rest-kmax)/rest; %turn it into a percentage 
            maxF = obj.Fmax;
            stiff = obj.SEstiff;
            tendon =  obj.TendonL;   %Length of artificial tendon and air fittings

           sF = zeros(size(obj.UnitDirection));
            
            %This function will solve for a new muscle length that
            %balances tendon force (fT) and the muscle force (fM) 
            function balanceF = muscleF(nmL)
                k = (rest-(mL-nmL-tendon-2.*fitting))./rest;  %strain 
                rel = k./KMAX;                             %relative strain  
                
               if dia == 10
                Fz = cell(1,4);
                [Fz{1}, Fz{2}, Fz{3}, Fz{4}] = normF10(rel, pres);
                Fn =Fz{3};
                fM = Fn.*maxF;
               elseif dia ~= 10
                fM = festo4(dia, rel, pres); 
               end

                fT = stiff*nmL;

                %Balance of the muscle forces. Solving to find when it
                %becomes equal to 0
                balanceF = fM - fT;
            end
            
          %Repeat the force calculation for every rotation of the joint
           for i = 1:size(unitD, 1)
                ML = obj.MuscleLength(i);              %Muscle-Tendon Length, which is the full calculated length between points in OpenSim

                %Set Function solver parameters
                options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);

                %Determine the normalized muscle length that solves the force
                %equations
                snmL = fsolve(@muscleF, ML, options);

                %With the solved muscle length, we can determine the scalar muscle
                %force by plugging it back into one of the force equations (fT)
                sF = stiff.*snmL;
                
           end

            for i = 1:size(unitD, 1)
                if sF(i) < 0
                    sF(i) = 0;
                end
                if sF(i) > maxF
                    sF(i) = NaN;
                end
            end
            
            SF = diag(sF);
            F = SF*unitD;

        end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
        function tor = get.Torque(obj)
            mA = obj.MomentArm;
            F = obj.Force;
            tor = cross(mA,F,2);
            
%             for i = 1:size(mA, 1)
%                 tor(i, :) = cross(mA(i, :), F(i, :));
%             end
        end    
    end
end