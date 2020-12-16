% Muscle Data
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Creation of a reusable class for categorizing muscle information and
%calculating parameters.

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html



classdef MonoMuscleData
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which column corresponds with a wrapping point in the Muscle Location
        MIF                         %Max Isometric Force
        TendonSlackLength           %Length of the tendon before it starts producing force
        PennationAngle              %Pennation angle of the muscle attaching to the tendon
        OptFiberLength              %Optimal fiber length of the muscle (where it generates the maximum force)
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
        function MD = MonoMuscleData(name, location, cross, mif, tsl, pennation, ofl, t)
            if nargin > 0
                MD.Name = name;
                MD.Location = location;
                MD.Cross = cross;
                MD.MIF = mif;
                MD.TendonSlackLength = tsl;
                MD.PennationAngle = pennation;
                MD.OptFiberLength = ofl; 
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
            mL = zeros(size(T, 3), 1);
            
            for ii = 1:size(mL, 1)                          %Repeat for each orientation
                for i = 1:size(L, 1)-1                      %Repeat for all muscle segments
                    pointA = L(i, :);
                    pointB = L(i+1, :);
                    if i+1 == C
                        pointB = RowVecTrans(T(:, :, ii), pointB);
                    end
                    mL(ii, 1) = mL(ii, 1) + norm(pointA - pointB);
                end
            end
        end
        
        %% -------------- Force Unit Direction ----------------
        %Calculate the unit direction of the muscle force about the joint.
        function unitD = computeUnitDirection(obj)
            L = obj.Location;
            T = obj.TransformationMat;
            C = obj.Cross;
            direction = zeros(size(T, 3), 3);
            unitD = zeros(size(direction));
            
            for i = 1:size(T, 3)
                pointA = L(C-1, :);
                pointB = L(C, :);
                direction(i, :) = RowVecTrans(T(:, :, i)\eye(4), pointA) - pointB;
                unitD(i, :) = direction(i, :)/norm(direction(i, :));
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
            mA = zeros(size(T, 3), 3);
            
            for i = 1:size(T, 3)
                pointB = L(C, :);
                mA(i, :) = pointB - unitD(i, :)*dot(unitD(i, :), pointB);
            end
        end
        
        %% -------------- Force --------------------------
        %Calculating muscle force, the algorithm comes from Hoy 1990,
        %Thelen 2003, and Millard 2013. Estimates the actual muscle length
        %along the segment length, then attempts to balance the force of
        %the tendon with the force of the pennated muscle to determine the
        %actual muscle length. Having a concrete muscle length provides the
        %actual force.
        
        function F = computeForce(obj)
            T = obj.TransformationMat;
            miF = obj.MIF;
            tsL = obj.TendonSlackLength;
            alphaP = obj.PennationAngle;
            ofL = obj.OptFiberLength;
            F = zeros(size(obj.UnitDirection));
            
            %This function will solve for a nominal muscle length that
            %balances tendon force (fT) and the muscle force ( [fPE + fL]*cos(pennation) ) 
            function balanceF = muscleF(nmL)
                %Values from OpenSim, slightly different than Thelen and Hoy
                gamma = 0.5;                    %Shape factor for force length curve of contractile element
                epsilon = 0.6;                  %Passive muscle strain due to maximum isometric force
                kPE = 4;                        %Exponential shape facor for passive muscle element

                ntsL = tsL/ofL;             %Normalized Tendon Slack Length
                nmtL = mtL/ofL;             %Normalized Musculo-Tendon Length

                %Force of the passive elastic element in the muscle
                %(Thelen eqn 3)
                fPE = (exp(kPE*(nmL - 1)/epsilon) - 1)/(exp(kPE) - 1);

                %Force of the active contractile element of the muscle
                %(Thelen eqn 4)
                fL = exp(-(nmL - 1)^2/gamma);

                %Calculation of the square root of the current pennation
                %angle (Hoy eqn 7
                cosAlpha = sqrt(1 - (sin(alphaP)/nmL)^2);

                %Force of the elastic tendon as it stretches (Hoy eqn 7)
                fT = 37.5/ntsL*(nmtL - nmL*cosAlpha - ntsL);

                %Balance of the muscle forces. Solving to find when it
                %becomes equal to 0
                balanceF = (fL + fPE)*cosAlpha - fT;
            end
            
            %Repeat the force calculation for every rotation of the joint
            for i = 1:size(T, 3)
                unitD = obj.UnitDirection(i, :);        %Direction of the force
                mtL = obj.MuscleLength(i);              %Muscle-Tendon Length, which is the full calculated length between points in OpenSim

                %Estimate the nominal muscle length, without the tendon. If it
                %is thought to be all tendon, the solver will crash
                if mtL == tsL
                    nmL0 = (mtL*1.01 - tsL)/ofL;
                else
                    nmL0 = (mtL - tsL)/ofL;
                end

                %Set Function solver parameters
                options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);

                %Determine the normalized muscle length that solves the force
                %equations
                snmL = fsolve(@muscleF, nmL0, options);

                %With the solved muscle length, we can determine the scalar muscle
                %force by plugging it back into one of the force equations (fT)
                sF = miF*37.5/(tsL/ofL)*(mtL/ofL - snmL*sqrt(1 - (sin(alphaP)/snmL)^2) - tsL/ofL);
                
                %If the force is less than 0, set it to 0 as a muscle can
                %only be in tension
                if sF < 0
                    sF = 0;
                end

                %With the scalar muscle force, we can multiply it by the unit
                %direction of the force to calculate the vectorized version
                F(i, :) = sF*unitD;
            end
        end
        
        %% -------------- Torque about a Joint --------------
        function tor = computeTorque(obj)
            mA = obj.MomentArm;
            F = obj.Force;
            tor = zeros(size(mA));
            
            for i = 1:size(mA, 1)
                tor(i, :) = cross(mA(i, :), F(i, :));
            end
        end    
    end
end
