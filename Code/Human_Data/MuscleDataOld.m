% Muscle Data
% Author: Connor Morrow
% Date: 1/14/2020
% Description: This script allows for creating reusable classes, which 
% categorize and calculates muscle information. This will be used in 
% determining muscle placement, optimization, and torque verification

%Creation of a reusable class for categorizing muscle information and
%calculating parameters.

%Refer to https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html



classdef MuscleDataOld
    
    %% ------------Public Properties---------------------------
    %List of explicit properties for the muscles
    properties
        Muscle                      %Name of the muscle
        Location
        ViaPoints                   %Designates which column corresponds with a wrapping point in the Muscle Location
        MIF                         %Max Isometric Force
        OFL                         %Optimal Fiber Length
        TSL                         %Tendon Slack Length
        PA                          %Pennation Angle
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
        function MD = MuscleData(muscle, location, viapoints, mif, ofl, tsl, pa, t, axis)
            if nargin > 0
                MD.Muscle = muscle;
                MD.Location = location;
                MD.ViaPoints = viapoints;
                MD.MIF = mif;
                MD.OFL = ofl;
                MD.TSL = tsl;
                MD.PA = pa;
                MD.TransformationMat = t;
                MD.Axis = axis;
                MD.MuscleLength = computeMuscleLength(MD);
                MD.MomentArm = computeMomentArm(MD);
                MD.Force = computeForce(MD);
                MD.Torque = computeTorque(MD);
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
        
        %% ---------------------- Muscle Force --------------
        % Reminder: Clean this code up if possible. Looks ugly after the
        % copy paste.
        function f = computeForce(obj)
            for i = 1:length(obj.MuscleLength)
                Lmt = obj.MuscleLength(i);
                mif = obj.MIF;
                ofl = obj.OFL;
                tsl = obj.TSL;
                pa = obj.PA;

                kPE = 4;    %passive exponential shape factor (OpenSim Kshape Pasive)
                eom = 0.6;  %passive muscle strain at max isometric force (OpenSim FmaxMuscleStrain)
                y = 0.5;    %active shape factor (OpenSim KshapeActive)
                


                if Lmt ~= tsl
                    x0 = (Lmt-tsl)/ofl; %Initial guess
                else
                    x0 = (Lmt*1.001-tsl)/ofl; %Initial guess, avoid computational crash
                end
                options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);
                Lma = fsolve(@myfunc,x0,options);
                
                f(i) = mif*(37.5/(tsl/ofl))*((Lmt/ofl)-(Lma)*((1-(sin(pa)/(Lma))^2)^(1/2))-(tsl/ofl)); %mif at Lmt
                if f < 0
                       f(i) = 0;  %Muscle force is tensile only (Millard 2013)
                end
            end
            
            function Fbal = myfunc(Lmn)

                Fpe = (exp(kPE*((Lmn)-1)/eom)-1)/(exp(kPE)-1); %Passive force-length curve, normalized (Thelen 2003)
                fL = exp(-(((Lmn)-1).^2/y));  %Active force-length curve, normalized (Thelen 2003)
                cosa = ((1-(sin(pa)/(Lmn)).^2)^(1/2)); %cosine of pennation angle (Hoy 1990)
                fT = (37.5/(tsl/ofl))*((Lmt/ofl)-(Lmn)*cosa-(tsl/ofl)); %Normalized tendon force (Hoy 1990)
                Fbal = (fL+Fpe)*cosa-fT;

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



