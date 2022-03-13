% Pam Data
% Author: Connor Morrow
% Date: 11/16/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPamDataPinnedExtensor < handle
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which row corresponds with a location where the muscle crosses into a new reference frame
        Diameter                    %Diameter of the BPA
        TransformationMat           %Contains a transformation matrix to change the 
        FittingLength
        TendonL
    end
    
    %Dependent properties are those that are calculated by the explicit
    %properties. Matlab will not calculate these until it is queried in the
    %main script
    properties (Dependent)   
        SegmentLengths
        LongestSegment
        MuscleLength
        RestingL
        Contraction
        LengthCheck
        UnitDirection
        MomentArm
        Force
        Torque
    end
    
   
    methods
        %% ------------- Muscle Data Constructor -----------------
        %Constructor Function. By calling 'MuscleData' and entering the
        %muscle information, we construct an object for that muscle.
        function PD = MonoPamDataPinnedExtensor(name, location, cross, diameter, t)
            if nargin > 0
                PD.Name = name;
                PD.Location = location;
                PD.Cross = cross;
                PD.Diameter = diameter;
                PD.TransformationMat = t;
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
            

%             for j = 1:size(T,3)
%                 if angle(j)>80
%                     L = obj.Location;
%                 elseif j>= 40 && j<= 80
%                     L(:,:,j) = [L(1,:);
%                          L(2,:);
%                          L(3,:);
%                          L(3,:);
%                          L(5,:);];
%                else 
%                     L(:,:,j) = [L(1,:);
%                          L(2,:);
%                          L(2,:);
%                          L(2,:);
%                          L(5,:);];
%                 end      
%             end
            
            for ii = 1:size(T, 3)                          %Repeat for each orientation
                for i = 1:size(L, 1,ii)-1                      %Repeat for all muscle segments
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
        function restingPamLength = get.RestingL(obj)
            
            restingPamLength = 0.460;
            fittingLength = 0.0254;
            tendonLength = 0;
            
            mL = obj.MuscleLength;
            dia = obj.Diameter;
            longestSeg = obj.LongestSegment;
           
            %Calculate the Pam end cap fitting length (estimates currently)
%             if dia == 20
%                 fittingLength = 0.02275;
%             elseif dia == 40
%                 fittingLength = 0.028;
%             else
%                 fittingLength = 0.022;
%             end
% 
            obj.FittingLength = fittingLength;
            
%             restingPamLength = max(mL) - 2*fittingLength-tendonLength;
%             
%             tendonLength = max(mL) - restingPamLength - 2*fittingLength;
%             if tendonLength < 0.08
%                 tendonLength = 0.08;
%             end
%             
%             %If there is only one muscle segment, make the segment length
%             %include the tendon length in the calculation for the resting
%             %PAM length
%             if size(obj.Location, 1) < 3
%                 restingPamLength = restingPamLength - tendonLength;   %The resting Pam length, needs to be accomodate the largest muscle length
%             end
            
            obj.TendonL = tendonLength;
        end
        
        %% -------------- Contraction of the PAM --------------------------
        function contraction = get.Contraction(obj)
            mL = obj.MuscleLength;
            restingPamLength = obj.RestingL;
            
            contraction = zeros(length(mL), 1);
            for i = 1:length(mL)
                contraction(i) = (max(mL) - mL(i))/restingPamLength;
            end
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
        
        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
        function F = get.Force(obj)
        %Inputs:
        %Lmt == muscle-tendon length, scalar
        %rest == resting length of artificial muscle, "size" from Size function
        %dia == diameter of Festo tube, from Size function
        %Outputs:
        %F == Force, N           
            dia = obj.Diameter;
            unitD = obj.UnitDirection;
            contract = obj.Contraction;
            contraction = obj.Contraction;
            mL = obj.MuscleLength;
            rest = obj.RestingL;
            long = max(mL);
            load ForceStrainTable.mat RelativeStrain Force2
            tendon = long - rest;   %Length of artificial tendon and air fittings

            k = (rest-(mL-tendon))/rest; %current strain
            act = [15.31; 18.28; 18.94; 19.50; 27.27; 28.09]/100; %Resting actuator lengths (Hunt 2017)
            strain = [0.1491; 0.1618; 0.1633; 0.1680; 0.1692; 0.1750]; %Max strain for these lengths (Hunt 2017)

            if rest >= max(act)
                kmax = max(strain);                 %maximum strain
            elseif rest <= min(act)
                kmax = min(strain);                 %maximum strain
            else
                kmax = interp1q(act,strain,rest);   %maximum strain
            end

            rel = k/kmax; %relative strain
    
            for i = 1:size(unitD, 1)
                if rel(i) >= 0 && rel(i) <=1
                    scalarForce(i) = interp1(RelativeStrain, Force2, rel(i));
                else
                    scalarForce(i) = 0;
                end
            end

            %If diameter is not 10 mm, then upscale force
            if dia == 20
                scalarForce = (1500/630)*scalarForce;
            end

            if dia == 40
                scalarForce = (6000/630)*scalarForce;
            end

            
            
%             if dia == 20
%                 x = [0, 0.07, 0.11, 0.15, 0.25]';
%                 y = [1400, 800, 600, 400, 0]';
%                 BPAFit = fit(x, y, 'poly2');
%             elseif dia == 40
%                 x = [0, 0.06, 0.12, 0.15, 0.25]';
%                 y = [6000, 3500, 2000, 1500, 0]';
%                 BPAFit = fit(x, y, 'poly2');
%             else
%                 x = [0, 0.1, 0.17, 0.25]';
%                 y = [630, 300, 150, 0]';
%                 BPAFit = fit(x, y, 'exp2');
%             end

%             scalarForce = BPAFit(contract);
            
            F = zeros(size(unitD));
            for i = 1:size(unitD, 1)
                if scalarForce(i) <= 0
                    scalarForce(i) = 0;
                end
                F(i, :) = unitD(i, :)*scalarForce(i);
            end
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
            tor = zeros(size(mA));
            
            for i = 1:size(mA, 1)
                tor(i, :) = cross(mA(i, :), F(i, :));
            end
        end    
    end
end