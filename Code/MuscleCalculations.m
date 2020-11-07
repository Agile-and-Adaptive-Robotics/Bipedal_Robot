% Muscle Data
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Creation of a reusable class for categorizing muscle information and
%calculating parameters.

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html



classdef MuscleData
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Name                        %Name of the muscle
        Location
        Cross                       %Designates which column corresponds with a wrapping point in the Muscle Location
        MIF                         %Max Isometric Force
        TransformationMat           %Contains a transformation matrix to change the 
        Axis                        %Contains the axes that the muscle will cause the joint to rotate about. 1>x, 2>y, 3>z
        MuscleLength
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
        function MD = MuscleData(name, location, cross, mif, t)
            if nargin > 0
                MD.Name = name;
                MD.Location = location;
                MD.Cross = cross;
                MD.MIF = mif;
                MD.TransformationMat = t;
                MD.MuscleLength = computeMuscleLength(MD);
%                 MD.MomentArm = computeMomentArm(MD);
%                 MD.Torque = computeTorque(MD);
            else
                fprintf('Invalid number of arguments\n')
            end
        end
        
        %% ------------- Muscle Length ------------------------
%         %Function that calculates the muscle length, based
        function mL = computeMuscleLength(obj)
            L = obj.Location;
            j = size(L, 2);
            T = obj.TransformationMat;
            C = obj.Cross
            mL = zeros(1, size(T, 3));     %Initialize the muscle length for each transformation matrix to start at 0
            
            for ii = 1:size(T, 3)
                for i = 1:size(L, 1)-1
                    pointA = L(i, :);
                    pointB = L(i+1, :);
                    if i+1 >= C
                        pointA = RowVecTrans(T(:, :, ii), pointA);
                    end
                    if i+1 >= C
                        L(i+1) = RowVecTrans(T(:, :, ii), pointB);
                    end
                    
                    mL(ii) = mL(ii) + norm(pointA - pointB);
                end
            end
            
            
            
            for iii = 1:size(obj.TransformationMat, 4)
                L = obj.Location;                   %Rename the location matrix to something simpler
                j = size(L, 2);                      %Determine the number of via points along the muscle
                T = obj.TransformationMat(:, :, :, iii);          %Rename the Transformation Tensor to something simpler

                wrap = obj.ViaPoints(1);                 %The first wrapping point
                ii = 1;                             %Index for the wrapping points
                t1 = eye(4);                        %Create an identity transformation matrix
                t2 = eye(4);                        %Create an identity transformation matrix
                for i = 1:j-1
                    while i + 1 == wrap                %Begin changing the transformation matrices if we've reached a viapoint
                        t2 = t2*T(:, :, ii);           %Store the first Transformation matrix in t2
                        ii = ii+1;                  %Update the pointer for the transformation matrix
                        if ii <= length(obj.ViaPoints)
                            wrap = obj.ViaPoints(ii);        %Set the next wrapping point, if ii isn't currently greater than the length of the via points vector
                        else
                            wrap = 0;                   %Cause the while loop to exit when the last wrapping point is used.
                        end
                    end

                    %Connor Note: The following algorithm matches Ben's, however, I
                    %believe that it might be incorrect, depending on where the
                    %transformation matrices are measured from.

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
        %For every ViaPoint, calculate the moment arm of the muscle about
        %the joint it crosses over
        function mA = computeMomentArm(obj)
            L = obj.Location;
            wrap = obj.ViaPoints;
            for ii = 1:size(obj.TransformationMat, 4)
                T = obj.TransformationMat(:, :, :, ii);
                for i = 1:length(obj.ViaPoints)
                    t1 = T(:, :, i);

                    %This is a mess to conform to how Ben's code works. Will need
                    %to review to make sure that it is accurate. It currently
                    %is set up to emulate the Extensor Digitorum Longus and how
                    %the transformation matrices transform the vectors. 
                    if i+1 > length(wrap)
                        t2 = eye(4);
                    else
                        if i == 1
                            if wrap(i) == wrap(i+1)
                                t2 = T(:, :, i+1);
                            else
                                t2 = eye(4);
                            end
                        else
                            if wrap(i) == wrap(i-1)
                                t1 = T(:, :, i-1)*T(:, :, i);
                                t2 = eye(4);
                            end
                        end
                    end

                    direction = VecTrans(t1\eye(4), L(:, wrap(i)-1))-VecTrans(t2, L(:, wrap(i)));     %Calculate the direction from the previous wrapping point to the next
                    unitDirection = direction/norm(direction);                                          %Calculate the unit direction of the direciton vector
                    
                    if length(obj.Axis) > length(wrap)                      %Check if the number of axes to investigate are greater than the number of wrapping points. In some of 
                                                                            %Ben's code, there might only be one wrapping point, but multiple axes for the moment arm. In order parts
                                                                            %of the code, there are multiple wrapping points, but only one axis
                        for j = 1:length(obj.Axis)
                            mA(j, ii) = CrossProd(VecTrans(t2, L(:, wrap(i))), unitDirection, obj.Axis(j), T(1:3, 1:3, i));
                        end
                    else
                        mA(i, ii) = CrossProd(VecTrans(t2, L(:, wrap(i))), unitDirection, obj.Axis(i), T(1:3, 1:3, i));
                    end
                    
                end
            end
        end
        
        
        %% ---------------------- Torque about a Joint --------------
        function tor = computeTorque(obj)
            for i = 1:size(obj.TransformationMat, 4)
                tor(:, i) = obj.MomentArm(:, i)*obj.Force(i);
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
    switch nargin
        case 3
            c = v(n);
        case 4
            u = Rot*v;
            c = u(n);
        otherwise
            c = 'Error calculating moment arm';
            error(c)
    end
end



