% Bi Pam Data
% Author: Connor Morrow
% Date: 11/23/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef BiPamData
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which row corresponds with a location where the muscle crosses into a new reference frame
        Diameter                    %Diameter of the BPA
        TransformationMat           %Contains a transformation matrix to change the 
        MuscleLength                %Total length of the muscle
        LongestSegment              %Longest contiguous section of a muscle. Where the BPA will reside
        UnitDirection               %The Unit direction of the force about a joint
        MomentArm                   %The moment arm from a joint to the force causing rotation about it
        RestingL
        Contraction
        TendonL
        LengthCheck
        Force                       %Calculate the maximum amount of force that the BPA can output at different orientations
        Torque
    end
    
    %Dependent properties are those that are calculated by the explicit
    %properties. Matlab will not calculate these until it is queried in the
    %main script
    properties (Dependent)        
    end
    
    methods
        %% ------------- Muscle Data Constructor -----------------
        %Constructor Function. By calling 'MuscleData' and entering the
        %muscle information, we construct an object for that muscle.
        function PD = MonoPamData(name, location, cross, diameter, t)
            if nargin > 0
                PD.Name = name;
                PD.Location = location;
                PD.Cross = cross;
                PD.Diameter = diameter;
                PD.TransformationMat = t;
                [PD.MuscleLength, PD.LongestSegment] = computeMuscleLength(PD);
                PD.UnitDirection = computeUnitDirection(PD);
                PD.MomentArm = computeMomentArm(PD);
                [PD.RestingL, PD.Contraction, PD.TendonL, PD.LengthCheck] = computePamLength(PD);
                PD.Force = computeForce(PD);
                PD.Torque = computeTorque(PD);
            else
                fprintf('Invalid number of arguments\n')
            end
        end

        %% ------------- Muscle Length ------------------------
        %Function that calculates the muscle length, based
        function [mL, longestSeg] = computeMuscleLength(obj)
            L = obj.Location;
            C = obj.Cross;
            T = obj.TransformationMat;
            
            mL = zeros(size(T, 3), 1);
            segLengths = zeros(size(T, 3), size(L, 1) - 1);
            
            for iii = 1:size(mL, 2)                         %Repeat for all iterations of the second joint
                for ii = 1:size(mL, 1)                      %Repeat for all iterations of the first joint
                    currentCross = 1;
                    for i = 1:size(L, 1)-1                  %Repeat for each muscle segment
                        pointA = L(i, :);
                        pointB = L(i+1, :);
                        if i+1 == C(currentCross)
                            if currentCross == 1
                                if C(currentCross) == C(currentCross + 1)
                                    pointB = RowVecTrans(T(:, :, ii, currentCross)*T(:, :, iii, currentCross + 1), pointB);
                                else
                                    pointB = RowVecTrans(T(:, :, ii, currentCross), pointB);
                                    currentCross = currentCross + 1;
                                end   
                            else
                                pointB = RowVecTrans(T(:, :, iii, currentCross), pointB);
                                currentCross = currentCross + 1;
                                if currentCross > length(C)
                                    currentCross = currentCross - 1;
                                end
                            end
                        end
                        segLengths(iii, i, ii) = norm(pointA - pointB);
                        mL(ii, iii) = mL(ii, iii) + segLengths(iii, i, ii);    
                    end                   
                end
            end
            
            % Calculate which muscle segment is the longest on average.
            % This will be where the Pam resides.
            avgSegL = zeros(size(L, 1) - 1);
            for i = 1:size(segLengths, 2)
                avgSegL(i) = mean(segLengths(:, i, :));
            end
            
            longestSegPointer = 1;
            if size(avgSegL, 1) > 1
                for i = 1:size(avgSegL, 1)
                    if avgSegL(i + 1) > avgSegL(i)
                        longestSegPointer = i + 1;
                    end
                end
            end
            longestSeg = segLengths(:, longestSegPointer, :);
        end
        
        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
        function unitD = computeUnitDirection(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            C = obj.Cross;
            direction = zeros(size(T, 3), 3);
            unitD = zeros(size(direction));
            
            for iii = 1:size(T, 3)                              %Rotation for second joint
                for ii = 1:size(T, 3)                           %Rotation for first joint
                    for i = 1:size(C, 2)
                        pointA = L(C(i)-1, :);
                        pointB = L(C(i), :);
                        if i > 1
                            if C(i - 1) == C(i)
                                direction(ii, :, iii, i) = RowVecTrans((T(:, :, ii, i)*T(:, :, iii, i))\eye(4), pointA) - pointB;
                            else
                                direction(ii, :, iii, i) = RowVecTrans(T(:, :, iii, i)\eye(4), pointA) - pointB;
                            end
                        else
                            if C(i) == C(i+1)
                                direction(ii, :, iii, i) = RowVecTrans(T(:, :, ii, i)\eye(4), pointA) - RowVecTrans(T(:, :, iii, i+1), pointB);
                            else
                                direction(ii, :, iii, i) = RowVecTrans(T(:, :, ii, i)\eye(4), pointA) - pointB;
                            end
                        end
                        unitD(ii, :, iii, i) = direction(ii, :, iii, i)/norm(direction(ii, :, iii, i));
                    end
                end
            end
        end
        
        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
        %For every ViaPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
        function mA = computeMomentArm(obj)
            T = obj.TransformationMat;
            L = obj.Location;
            C = obj.Cross;
            unitD = obj.UnitDirection;
            mA = zeros(size(T, 3), 3, size(T, 4), 2);
            
            for iii = 1:size(T, 3)                          %Rotation for second joint
                for ii = 1:size(T, 3)                       %Rotation for first joint
                    for i = 1:size(C, 2)                    %Repeat for each crossing point
                        pointB = L(C(i), :);
                        u = unitD(ii, :, iii, i);  
                        if i > 1
                            mA(ii, :, iii, i) = pointB - u*dot(u, pointB);
                        else
                            if C(i + 1) == C(i)
                                pointBp = RowVecTrans(T(:, :, iii, i+1), pointB);
                                mA(ii, :, iii, i) = pointBp - u*dot(u, pointBp);
                            else
                                mA(ii, :, iii, i) = pointB - u*dot(u, pointB);
                            end
                        end
                    end
                end
            end
        end
        
        %% -------------- Pam Length --------------------------
        %This function calculates the length of the Pam that will be used
        %for this muscle. The BPA will be placed in the longest segment of
        %the muscle. The rest of the muscle segments will be considered
        %tendons for attaching to the robot.
        %This also checks to make sure that the BPA fully contracted length
        %is short enough to achieve the full joint rotation.
        function [restingPamLength, contraction, tendonLength, lengthCheck] = computePamLength(obj)
            dia = obj.Diameter;
            longestSeg = obj.LongestSegment;
            mL = obj.MuscleLength;
            contractPercent = 0.25;                     %Maximum contraction percentage of a BPA
           
            %Calculate the Pam end cap fitting length (estimates currently)
            if dia == 20
                fittingLength = 0.025;
            elseif dia == 40
                fittingLength = 0.05;
            else
                fittingLength = 0.0125;
            end
            
            restingPamLength = max(longestSeg) - 2*fittingLength;   %The resting Pam length, needs to be accomodate the largest muscle length
            tendonLength = max(mL) - restingPamLength - 2*fittingLength;
            
            % Make sure there is at least 4 cm of tendon for connections
            if tendonLength < 0.04
                tendonLength = 0.04;
                restingPamLength = restingPamLength - tendonLength;
            end

            contraction = zeros(length(mL), 1);
            for i = 1:length(mL)
                contraction(i) = (mL(i) - tendonLength - 2*fittingLength)/restingPamLength;
            end
            
            if min(contraction) >= (1 - contractPercent)
                lengthCheck = 'Usable';
            else
                lengthCheck = 'Unusable';
            end
        end
        
        %% -------------- Force --------------------------
        %Calculate the directin of the forced applied by the muscle
        function F = computeForce(obj)
            dia = obj.Diameter;
            unitD = obj.UnitDirection;
            
            if dia == 20
                maxForce = 1500;
            elseif dia == 40
                maxForce = 6000;
            else
                maxForce = 630;
            end
            
            F = unitD*maxForce;
        end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
        function tor = computeTorque(obj)
            mA = obj.MomentArm;
            F = obj.Force;
            T = obj.TransformationMat;
            C = obj.Cross;
            tor = zeros(size(mA));
            
            for iii = size(T, 3)                    %Rotation for the second joint
                for ii = size(T, 3)                 %Rotation for the first joint
                    for i = 1:size(C, 2)            %Repeat for eaching crossing point
                        tor(ii, :, iii, i) = cross(mA(ii, :, iii, i), F(ii, :, iii, i));
                    end
                end
            end            
        end    
    end
end