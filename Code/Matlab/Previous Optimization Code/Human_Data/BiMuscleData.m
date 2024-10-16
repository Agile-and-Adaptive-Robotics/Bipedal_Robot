% Muscle Data
% Author: Connor Morrow
% Date: 11/12/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification.
% Specifically for Biarticular Muscles

%Creation of a reusable class for categorizing muscle information and
%calculating parameters.

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html



classdef BiMuscleData
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which column corresponds with a wrapping point in the Muscle Location
        MIF                         %Max Isometric Force
        TransformationMat           %Contains a transformation matrix to change the 
        MuscleLength
        UnitDirection               %Unit direction of the force about the joint.
        MomentArm
        Force                       %Force Generated in the muscle
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
        function MD = BiMuscleData(name, location, cross, mif, t)
            if nargin > 0
                MD.Name = name;
                MD.Location = location;
                MD.Cross = cross;
                MD.MIF = mif;
                MD.TransformationMat = t;
                MD.MuscleLength = computeMuscleLength(MD);
                MD.UnitDirection = computeUnitDirection(MD);
                MD.Force = computeForce(MD);
                MD.MomentArm = computeMomentArm(MD);
                MD.Torque = computeTorque(MD);
            else
                fprintf('Invalid number of arguments\n')
            end
        end
        
        %% ------------- Muscle Length ------------------------
%         %Function that calculates the muscle length, based
        function mL = computeMuscleLength(obj)
            L = obj.Location;
            C = obj.Cross;
            T = obj.TransformationMat;
            
            mL = zeros(size(T, 3), size(T, 3));

            
            for iii = 1:size(mL, 2)                         %Repeat for all iterations of the second joint
                for ii = 1:size(mL, 1)                      %Repeat for all iterations of the first joint
                    currentCross = 1;
                    for i = 1:size(L, 1)-1
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
                        mL(ii, iii) = mL(ii, iii) + norm(pointA - pointB);    
                    end                   
                end
            end
        end
        
        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
        function unitD = computeUnitDirection(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            C = obj.Cross;
            direction = zeros(size(T, 3), 3, size(T, 4), 2);
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
            
            for iii = 1:size(T, 3)
                for ii = 1:size(T, 3)
                    for i = 1:size(C, 2)
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
        
        %% -------------- Force --------------------------
        %Calculate the directin of the forced applied by the muscle
        function F = computeForce(obj)
            mif = obj.MIF;
            unitD = obj.UnitDirection;
            
            F = unitD*mif;
        end
        
        %% -------------- Torque about a Joint --------------
        function tor = computeTorque(obj)
            mA = obj.MomentArm;
            F = obj.Force;
            T = obj.TransformationMat;
            C = obj.Cross;
            tor = zeros(size(mA));
            
            for iii = 1:size(T, 3)
                for ii = 1:size(T, 3)
                    for i = 1:size(C, 2)
                        tor(ii, :, iii, i) = cross(mA(ii, :, iii, i), F(ii, :, iii, i));
                    end
                end
            end              
        end    
    end
end
