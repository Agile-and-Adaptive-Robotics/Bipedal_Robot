% Pam Data
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates PAM muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

classdef PamData
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Muscle                      %Name of the muscle
        Location
        CrossPoints                 %Designates which column corresponds with a location where the muscle crosses into a new reference frame
        MIF                         %Max Isometric Force
        TransformationMat           %Contains a transformation matrix to change the 
        Axis                        %Contains the axes that the muscle will cause the joint to rotate about. 1>x, 2>y, 3>z
        MuscleLength                %Total length of the muscle
        MomentArm
        Diameter
        RestingL
        LongestL
        Force        
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
        function PD = PamData(muscle, location, points, mif, t, axis)
            if nargin > 0
                PD.Muscle = muscle;
                PD.Location = location;
                PD.CrossPoints = points;
                PD.MIF = mif;
                PD.TransformationMat = t;
                PD.Axis = axis;
                PD.MuscleLength = computeMuscleLength(PD);
                PD.MomentArm = computeMomentArm(PD);
                [PD.Diameter, PD.RestingL, PD.LongestL] = Size(PD.MuscleLength, PD.MIF);
                PD.Force = festo(PD);
                PD.Torque = computeTorque(PD);
            else
                fprintf('Invalid number of arguments\n')
            end
        end

        %% ------------- Muscle Length ------------------------
%         %Function that calculates the muscle length, based
        function mL = computeMuscleLength(obj)
            mL = zeros(1, size(obj.TransformationMat, 4));     %Initialize the muscle length for each transformation matrix to start at 0
            for iii = 1:size(obj.TransformationMat, 4)
                L = obj.Location;                   %Rename the location matrix to something simpler
                j = size(L, 2);                      %Determine the number of via points along the muscle
                T = obj.TransformationMat(:, :, :, iii);          %Rename the Transformation Tensor to something simpler

                cross = obj.CrossPoints(1);                 %The first crossing point
                ii = 1;                             %Index for the wrapping points
                t1 = eye(4);                        %Create an identity transformation matrix
                t2 = eye(4);                        %Create an identity transformation matrix
                for i = 1:j-1
                    while i + 1 == cross                %Begin changing the transformation matrices if we've reached a crossing point
                        t2 = t2*T(:, :, ii);           %Store the first Transformation matrix in t2
                        ii = ii+1;                  %Update the pointer for the transformation matrix
                        if ii <= length(obj.CrossPoints)
                            cross = obj.CrossPoints(ii);        %Set the next wrapping point, if ii isn't currently greater than the length of the via points vector
                        else
                            cross = 0;                   %Cause the while loop to exit when the last wrapping point is used.
                        end
                    end

                    %Calculate the euclidean distance between the two via
                    %points, including the transformation matrix if necessary
                    mL(iii) = mL(iii) + norm(VecTrans(t1, L(:, i))-VecTrans(t2, L(:, i+1)));
    %                 t1 = t2;                        %If the forward point was transformed, store the transformation for the next iteration
                    t2 = t1;                        %Replace the transformtaion matrix with an identity matrix
                end         
            end
        end
        
        %% -------------- Moment Arm --------------------------
        %Calculate the moment arm about a joint
        %For every CrossPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
        function mA = computeMomentArm(obj)
            mA = zeros(length(obj.CrossPoints), size(obj.TransformationMat, 4), size(obj.Axis, 2));
            for ii = 1:size(obj.TransformationMat, 4)       %Repeat calculation for every degree of motion we are observing
                for i = 1:length(obj.CrossPoints)           %Repeat calculation for every joint the muscle will cross over. 
                    T = obj.TransformationMat(:, :, :, ii);
                    L = obj.Location;
                    cross = obj.CrossPoints;

                    %Setting up transformation matrix for the point before
                    %the cross over
                    if i == 1
                        t1 = T(:, :, i);                %if this is the first crossing point, the transformation matrix has to be the first one listed. We do not need to check if the muscle has passed over multiple joints
                    elseif cross(i) == cross(i-1)      
                        t1 = T(:, :, i-1)*T(:, :, i);   %if the new crossing point is the same as the last one, that means that the muscle is spanning two joints. Therefore, we need two transformation matrices to bring the point forward into the current frame
                    else
                        t1 = T(:, :, i);                %if the new crossing point is not the first one or it doesn't span multiple joints, the transformation matrix is good as is.
                    end

                    %Setting up transformation matrix for the point after
                    %the cross over
                    if i == length(cross)                %If this is the last crossing point, then the second cross over point is already in the correct frame
                        t2 = eye(4);
                    elseif cross(i) == cross(i+1)       %If the next crossing point is the same as the current one, then the muscle spans two joints. We will need the next transformation matrix to pull the point into the current frame
                        t2 = T(:, :, i+1);
                    else
                        t2 = eye(4);                    %Otherwise, the second point is in the correct frame already.
                    end

                    direction = VecTrans(t1\eye(4), L(:, cross(i)-1))-VecTrans(t2, L(:, cross(i)));     %Calculate the direction from the previous crossing point to the next
                    unitDirection = direction/norm(direction);                                          %Calculate the unit direction of the direciton vector
                    for iii = 1:size(obj.Axis, 2)       %Repeat the moment arm calculation for every axis of interest 
                        if obj.Axis(i, iii) > 0         %Do not calculate the moment arm if the axis is listed as 0, which can happen if we are interested in multiple axes for one joint, but only one axis for the next joint that the muscles crosses
                            mA(i, ii, iii) = CrossProd(VecTrans(t2, L(:, cross(i))), unitDirection, obj.Axis(iii), T(1:3, 1:3, i)); %Cross the distance vector to the
                        end
                    end
                                    
                end
            end
        end
        
        %% ---------------------- Torque --------------
        function tor = computeTorque(obj)
            tor = zeros(length(obj.CrossPoints), size(obj.TransformationMat, 4), size(obj.Axis, 2));
            for i = 1:size(obj.TransformationMat, 4)                            %Repeat Calculation for every axis
                for ii = 1:length(obj.Axis)
                    tor(:, i, ii) = obj.MomentArm(:, i, ii)*obj.Force(i);
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
            if size+2*x > long && def <= max_def
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