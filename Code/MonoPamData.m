% Pam Data
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef MonoPamData
    
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
                PD.Force = festo(PD);
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
            segLenghts = zeros(size(T, 3), size(L, 1) - 1);
            
            for ii = 1:size(mL, 1)                          %Repeat for each orientation
                for i = 1:size(L, 1)-1                      %Repeat for all muscle segments
                    pointA = L(i, :);
                    pointB = L(i+1, :);
                    if i+1 == C
                        pointB = RowVecTrans(T(:, :, ii), pointB);
                    end
                    segLength(ii, i) = norm(pointA - pointB);
                    mL(ii, 1) = mL(ii, 1) + segLength;
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
            mif = obj.MIF;
            unitD = obj.UnitDirection;
            
            F = unitD*mif;
        end
        
        %% ---------------------- Torque --------------
        %Calculate torque by multiplying the the force along the 
        %Useful information
        % i -> Index for Crossing Points/Joints
        % ii -> Index for every degree of motion
        % iii -> Index for axes of interest to observe Torque about
        function tor = computeTorque(obj)
            tor = zeros(length(obj.CrossPoints), size(obj.TransformationMat, 4), size(obj.Axis, 2));
            for ii = 1:size(obj.TransformationMat, 4)                            %Repeat Calculation for every degree of motion
                for iii = 1:size(obj.Axis, 2)                                    %Repeat Calculation for every axis of interest to observe torque about
                    tor(:, ii, iii) = obj.MomentArm(:, ii, iii)*obj.Force(ii);
                end
            end
        end    
 
    end
end

%% ---------------- General Functions -----------------
%Function that rotates a vector by a transformation matrix
function v = VecTrans(T, p)
    r = [p; 1];
    s = T*r;
    v = [s(1); s(2); s(3)];
end

%Function for doing cross products and only looking at the axis of
%interest
function c = CrossProd(A, B, n, Rot)
    v = cross(A, B);
    
    NoRotate = 0;
    if n > 9                %When n is 10, 20, or 30, it indicates that the resulting cross product should not be rotated by the rotation matrix
        n = n/10;           %divide the axis of interest by 10 to get 1, 2, or 3 as an axis
        NoRotate = 1;              %indicate that the vector should not be rotated by the rotation matrix
    end
    
    if n == 0           %if the axis of interest is listed as 0, do not perform a rotation
        Rot = eye(3);
    end
    
    switch nargin
        case 3
            c = v(n);
        case 4
            if NoRotate == 1
                c = v(n);
            else
                u = Rot*v;
                c = u(n);
            end
        otherwise
            c = 'Error calculating moment arm';
            error(c)
    end
end

%Function for sizing Festo Muscles based on how much force and length
%change is necessary
function [ID, size, longest]= Size(length, mif)
    %function Size for sizing of Festo artificial muscles
    %Inputs:
    %Length is the musculotendon length, either a column vector or a matrix
    %mif is Maximum Isometric Force
    %Outputs:
    %ID is inside diameter of the specified festo muscle
    %size is the length of the festo muscle
    %long is the longest artificial musculo-tendon length
    %Variables:
    %Fc is a percentage of max theoretical festo muscle force
    %x is length of fitting from end to attachment point.

    %Notes:
    %1) If error calculating ID based on force, then two muscles in parallel
    %should be considered.
    %2) "Size Calc error1" if no change between shortest and largest lengths
    %for give dof(s). 
    %3) "Size Calc error2" if error calculating size because delta is too large 
    %compared to smallest length, then consider disregarding muscle or changing 
    %attachment points(e.g. short interior groin and hip muscles).

    long = max(length);
    short = min(length);
    delta = long - short; 
    F1 = 630;             %Maximum theoretical force for 10mm I.D. Festo muscle
    F2 = 1500;            %Maximum theoretical force for 20mm I.D. Festo muscle
    F3 = 6000;            %Maximum theoretical force for 40mm I.D. Festo muscle
    Fc = 0.9;          
    x = 0.0125;           %Length of air fittings

    if mif < Fc*F1
        ID = 10;
    elseif mif >= Fc*F1 && mif < Fc*F2
        ID = 20;
    elseif mif >= Fc*F2 && mif <= Fc*F3
        ID = 40;
    else
        ID = 'Use additional muscles';
    end

    opt_max = 1.09;  %Maximum allowable length change for optimal performance (Festo)
    %opt_min = 1;
    max_def = 1.25;  %Maximum allowable length change (Festo)


    size = delta/(1 - 1/opt_max);
    if delta == 0
        size = 'Size calc error1';
        return
    end    

    if size+2*x > long
        for def = opt_max:0.001:max_def
            if def <= max_def
                size = delta/(1 - 1/def);
            elseif size+2*x <= long
                def;
                break
    %         elseif def == max_def
    %             def
    %             size = 'Size calc error2';
            end
        end
    end

    if nargout > 2
        longest = long;
    end

end

%Function that generates the force of the Festo Muscle
function F = festo(obj)
    %Inputs:
    %Lmt == muscle-tendon length, scalar
    %rest == resting length of artificial muscle, "size" from Size function
    %dia == diameter of Festo tube, from Size function
    %Outputs:
    %F == Force, N
    
    Lmt = obj.MuscleLength;
    rest = obj.RestingL;
    dia = obj.Diameter;
    long = obj.LongestL;
    
    load ForceStrainTable.mat RelativeStrain Force
    tendon = long - rest;   %Length of artificial tendon and air fittings
    act = [15.31; 18.28; 18.94; 19.50; 27.27; 28.09]; %Resting actuator lengths (Hunt 2017)
    strain = [0.1491; 0.1618; 0.1633; 0.1680; 0.1692; 0.1750]; %Max strain for these lengths (Hunt 2017)
    
    F = zeros(length(Lmt),1);
    for i = 1:length(Lmt)
        k = (rest-(Lmt(i)-tendon))/rest; %current strain


        if rest >= max(act)
            kmax = max(strain);                 %maximum strain
        elseif rest <= min(act)
            kmax = min(strain);                 %maximum strain
        else
            kmax = interp1q(act,strain,rest);   %maximum strain
        end

        rel = k/kmax; %relative strain;

        if rel >= 0 && rel <=1
            F(i) = interp1(RelativeStrain, Force, rel);
        else
            F(i) = 0;
        end

        %If diameter is not 10 mm, then upscale force
        if dia == 20
            F(i) = (1500/630)*F(i);
        elseif dia == 40
            F(i) = (6000/630)*F(i);
        end
    end
    %Reference:
    %Hunt, Alexander J., Alexander Graber-Tilton, and Roger D. Quinn. "Modeling length effects of braided pneumatic actuators."
    %In ASME 2017 International Design Engineering Technical Conferences and Computers and Information in Engineering Conference, pp. V05AT08A008-V05AT08A008.
    %American Society of Mechanical Engineers, 2017.
end