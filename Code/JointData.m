classdef JointData
    
    properties
        Joint                   %Name of the Joint
        RotAxis                 %Axis of rotation for joint
        MaxTheta                %Maximum rotation for joint
        MinTheta                %Minimum rotation for the joint
        Home                    %Home position vector for the joint
        Divisions               %Number of points between min and max theta
        Theta                   %Current angle of the joint
        RotMat                  %Rotation matrix for the joint
        TransformationMat       %Transformation matrix for the join
    end
    
    properties (Dependent)

    end
    
    methods
        %% ----------------- Joint Data Constructor
        function JD = JointData(joint, rotaxis, maxtheta, mintheta, home, divisions)
            if nargin > 0
                JD.Joint = joint;
                JD.RotAxis = rotaxis;
                JD.MaxTheta = maxtheta;
                JD.MinTheta = mintheta;
                JD.Home = home;
                JD.Divisions = divisions;
                JD.Theta = computeTheta(JD);
                JD.RotMat = computeRotMat(JD);
                JD.TransformationMat = computeTransformationMat(JD);
            else
                fprintf('Invalid number of arguments\n')
            end
        end
        
        %% --------------- Create the Linearly Space Theta values
        function theta = computeTheta(obj)
            theta = linspace(obj.MinTheta, obj.MaxTheta, obj.Divisions);
        end
        
        
        %% ------------------ Create Rotation Matrix
        %Generate a rotation matrix for the joint. Based off of Modern
        %Robotics, pp. 74
        function rotmat = computeRotMat(obj)
            theta = obj.Theta;
            ra = obj.RotAxis;
            rotmat = zeros(3, 3, length(theta));
            
            for i = 1:length(theta)
                rotmat(1, 1, i) = cos(theta(i)) + ra(1)^2*(1-cos(theta(i))); 
                rotmat(1, 2, i) = ra(1)*ra(2)*(1-cos(theta(i)))-ra(3)*sin(theta(i));
                rotmat(1, 3, i) = ra(1)*ra(3)*(1-cos(theta(i))) + ra(2)*sin(theta(i));

                rotmat(2, 1, i) = ra(1)*ra(2)*(1-cos(theta(i))) + ra(3)*sin(theta(i));
                rotmat(2, 2, i) = cos(theta(i)) + ra(2)^2*(1-cos(theta(i)));
                rotmat(2, 3, i) = ra(2)*ra(3)*(1-cos(theta(i))) - ra(1)*sin(theta(i));

                rotmat(3, 1, i) = ra(1)*ra(3)*(1-cos(theta(i))) - ra(2)*sin(theta(i));
                rotmat(3, 2, i) = ra(2)*ra(3)*(1-cos(theta(i))) + ra(1)*sin(theta(i));
                rotmat(3, 3, i) = cos(theta(i)) + ra(3)^2*(1-cos(theta(i)));
            end
            
        end
        
        %% ------------Create Transformation Matrix
        function tm = computeTransformationMat(obj)
            theta = obj.Theta;
            tm = zeros(4, 4, length(theta));
            rotmat = obj.RotMat;
            h = obj.Home;
            
            for i = 1:length(theta)
                tm(:, :, i) = [rotmat(:, :, i), h; 0 0 0 1];
            end
        end      
           
    end
end
