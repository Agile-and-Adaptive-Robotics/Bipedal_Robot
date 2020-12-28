% Bi Pam Data
% Author: Connor Morrow
% Date: 11/23/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef BiPamData < handle
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which row corresponds with a location where the muscle crosses into a new reference frame
        Diameter                    %Diameter of the BPA
        TransformationMat           %Contains a transformation matrix to change the 
        FittingLength
    end
    
    %Dependent properties are those that are calculated by the explicit
    %properties. Matlab will not calculate these until it is queried in the
    %main script
    properties (Dependent)  
        SegmentLengths
        LongestSegment
        MuscleLength
        RestingL
        TendonL
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
        function PD = BiPamData(name, location, cross, diameter, t)
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
        
        %% ------------- Segment Lengths ----------------------------------
        function segLengths = get.SegmentLengths(obj)
            L = obj.Location;
            C = obj.Cross;
            T = obj.TransformationMat;
            segLengths = zeros(size(T, 3), size(T, 3), size(L, 1)-1);
            
            for iii = 1:size(T, 4)                         %Repeat for all iterations of the second joint
                for ii = 1:size(T, 3)                      %Repeat for all iterations of the first joint
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
                        segLengths(ii, iii, i) = norm(pointA - pointB);  
                    end                   
                end
            end
        end
        
        %% ------------- Longest Segment Calculation ----------------------
        function longestSeg = get.LongestSegment(obj)
            L = obj.Location;
            segLengths = obj.SegmentLengths;
            
            % Calculate which muscle segment is the longest on average.
            % This will be where the Pam resides.
            avgSegL = zeros(size(L, 1) - 1);
            for i = 1:size(segLengths, 3)
                avgSegL(i) = mean(mean(segLengths(:, :, i)));
            end
            
            longestSegPointer = 1;
            if size(avgSegL, 1) > 1
                for i = 1:size(avgSegL, 1)-1
                    if avgSegL(i + 1) > avgSegL(i)
                        longestSegPointer = i + 1;
                    end
                end
            end
            longestSeg = segLengths(:, :, longestSegPointer);
        end
            
        %% ------------- Muscle Length ------------------------------------
        %Function that calculates the muscle length, based
        function mL = get.MuscleLength(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            
            mL = zeros(size(T, 3), size(T, 3));
            segLengths = obj.SegmentLengths;
            
            for iii = 1:size(mL, 2)
                for ii = 1:size(mL, 1)
                    for i = 1:size(L, 1) - 1
                        mL(ii, iii) = mL(ii, iii) + segLengths(ii, iii, i);
                    end
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
        function mA = get.MomentArm(obj)
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
        
        %% -------------- Resting PAM Length ------------------------------
        function restingPamLength = get.RestingL(obj)
            dia = obj.Diameter;
            longestSeg = obj.LongestSegment;
            
            %Calculate the Pam end cap fitting length (estimates currently)
            if dia == 20
                fittingLength = 0.025;
            elseif dia == 40
                fittingLength = 0.05;
            else
                fittingLength = 0.0125;
            end
            
            obj.FittingLength = fittingLength;
            restingPamLength = max(max(longestSeg)) - 2*fittingLength;   %The resting Pam length, needs to be accomodate the largest muscle length
        end
        
        %% -------------- Tendon Length -----------------------------------
        function tendonLength = get.TendonL(obj)
            mL = obj.MuscleLength;
            restingPamLength = obj.RestingL;
            fittingLength = obj.FittingLength;
            
            tendonLength = max(max(mL)) - restingPamLength - 2*fittingLength;
            
            % Make sure there is at least 4 cm of tendon for connections
            if tendonLength < 0.04
                tendonLength = 0.04;
                obj.RestingL = restingPamLength - tendonLength;
            end
        end
        
        %% -------------- Contraction of the PAM --------------------------
        function contraction = get.Contraction(obj)
            mL = obj.MuscleLength;
            restingPamLength = obj.RestingL;
            
            contraction = zeros(size(mL));
            for i = 1:size(mL, 1)
                for ii = 1:size(mL, 2)
                    contraction(i, ii) = (max(max(mL)) - mL(i, ii))/restingPamLength;
                end
            end
        end
        
        %% -------------- Length Check ------------------------------------
        function lengthCheck = get.LengthCheck(obj)
            contraction = obj.Contraction;
            contractPercent = 0.25;
            
            if min(min(contraction)) <= contractPercent
                lengthCheck = 'Usable';
            else
                lengthCheck = 'Unusable';
            end
        end
        
        %% -------------- Force --------------------------
        %Calculate the direction of the forced applied by the muscle
        function F = get.Force(obj)
            dia = obj.Diameter;
            unitD = obj.UnitDirection;
            contract = obj.Contraction;
            C = obj.Cross;
            
            if dia == 20
                x = [0, 0.07, 0.11, 0.15, 0.25]';
                y = [1400, 800, 600, 400, 0]';
                BPAFit = fit(x, y, 'poly3');
            elseif dia == 40
                x = [0, 0.06, 0.12, 0.15, 0.25]';
                y = [6000, 3500, 2000, 1500, 0]';
                BPAFit = fit(x, y, 'poly3');
            else
                x = [0, 0.1, 0.17, 0.25]';
                y = [630, 300, 150, 0]';
                BPAFit = fit(x, y, 'exp2');
            end

            scalarForce = BPAFit(contract);
            
            pos = 1;
            F = zeros(size(unitD));
            for iii = 1:size(unitD, 3)                  %Repeat for second joint rotation
                for ii = 1:size(unitD, 1)               %Repeat for first joint rotation
                    for i = 1:size(C, 2)                %Repeat for eaching crossing point
                        if scalarForce(pos) <= 0
                            scalarForce(pos) = 0;
                        end
                        F(ii, :, iii, i) = unitD(ii, :, iii, i)*scalarForce(pos);
                    end
                    pos = pos + 1;
                end
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
            T = obj.TransformationMat;
            C = obj.Cross;
            tor = zeros(size(mA));
            
            for iii = 1:size(T, 3)                    %Rotation for the second joint
                for ii = 1:size(T, 3)                 %Rotation for the first joint
                    for i = 1:size(C, 2)            %Repeat for eaching crossing point
                        tor(ii, :, iii, i) = cross(mA(ii, :, iii, i), F(ii, :, iii, i));
                    end
                end
            end            
        end    
    end
end