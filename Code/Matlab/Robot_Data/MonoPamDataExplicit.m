% Pam Data
% Author: Connor Morrow
% Date: 11/16/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPamDataExplicit < handle
    
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
        Force
        Torque
    end
    
   
    methods
        %% ------------- Muscle Data Constructor -----------------
        %Constructor Function. By calling 'MuscleData' and entering the
        %muscle information, we construct an object for that muscle.
        function PD = MonoPamDataExplicit(name, location, cross, diameter, t, rest, kmax, tendon, fit, pres)
            if nargin == 10
                PD.Name = name;
                PD.Location = location;
                PD.Cross = cross;
                PD.Diameter = diameter;
                PD.TransformationMat = t;
                PD.RestingL = rest;
                PD.Kmax = kmax;
                PD.TendonL = tendon;
                PD.FittingLength = fit;
                PD.Pressure = pres;
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
%         function restingPamLength = get.RestingL(obj)
%             
%             restingPamLength = 0.415;
%             fittingLength = 0.0254;
%             tendonLength = 0;
%             
%             mL = obj.MuscleLength;
%             dia = obj.Diameter;
%             longestSeg = obj.LongestSegment;
%            
%             %Calculate the Pam end cap fitting length (estimates currently)
% %             if dia == 20
% %                 fittingLength = 0.02275;
% %             elseif dia == 40
% %                 fittingLength = 0.028;
% %             else
% %                 fittingLength = 0.022;
% %             end
% % 
%             obj.FittingLength = fittingLength;
%             
% %             restingPamLength = max(mL) - 2*fittingLength-tendonLength;
% %             
% %             tendonLength = max(mL) - restingPamLength - 2*fittingLength;
% %             if tendonLength < 0.08
% %                 tendonLength = 0.08;
% %             end
% %             
% %             %If there is only one muscle segment, make the segment length
% %             %include the tendon length in the calculation for the resting
% %             %PAM length
% %             if size(obj.Location, 1) < 3
% %                 restingPamLength = restingPamLength - tendonLength;   %The resting Pam length, needs to be accomodate the largest muscle length
% %             end
%             
%             obj.TendonL = tendonLength;
%         end
        
        %% -------------- Contraction of the PAM --------------------------
        function contraction = get.Contraction(obj)
            mL = obj.MuscleLength;
            rest = obj.RestingL;
            tendon = obj.TendonL;
            fitting = obj.FittingLength;
            
            contraction = zeros(length(mL), 1);
            for i = 1:length(mL)
%                 contraction(i) = (rest-(mL(i,1)-tendon-2*fitting))/rest;
                  contraction(i) = ((rest+tendon+2*fitting)-mL(i,1))/(rest+tendon+2*fitting);
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
            kmax = (rest-kmax)/rest; %turn it into a percentage 
            
           load ForceStrainTable.mat ForceStrain
           tendon =  obj.TendonL;   %Length of artificial tendon and air fittings

           X = linspace(0,620,20); %Pressure for interpolation
           Y = linspace(0,1,30);   %Relative strain range for interpolation
           
           k = zeros(size(unitD,1),1);
           rel = zeros(size(unitD,1),1);
           scalarForce = zeros(size(unitD,1),1);
           for i = 1:size(unitD, 1)
              if dia == 10
                k(i) = (rest-(mL(i,1)-tendon-2*fitting))/rest; %current strain 
                rel(i) = k(i)/kmax; %relative strain
                if contract(i) < 0 && contract(i) >= -0.03
                    scalarForce(i) = interp1([-0.03 -0.015 0], [630 630 510], contract(i),'linear');
                elseif rel(i) >= 0 && rel(i) <= 1
                    scalarForce(i) = interp2(X, Y, ForceStrain, pres, rel(i), 'linear');
                elseif rel(i) > 1
                    scalarForce(i) = 0;
                else
                    scalarForce(i) = NaN;
                end
             
              else %If diameter is not 10 mm, then use Festo Lookup table
                scalarForce(i) = festo4(dia, pres, contract(i));
              end
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