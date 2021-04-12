%% PAM Placement Cost Function, for optimization about the knee
% This is the cost function for trying to determine the optimal place to
% attachment pneumatic artificial muscles on a humanoid robot. This
% function calculates the error between the robot generated torque and the
% human torque. The more positive the number, the better the robot is at
% creating equivalent or greater torque values than the human. In order to
% capture the negative torque contribution that some muscles create, the
% function switches between subtracting the human torque from the robot
% torque, and subtracting the robot torque from the human torque. This is
% because we want the robot to be able to produce more positive torque when
% the human torque is positive, and we want the robot to produce more
% negative torque when the human torque is negative. this cannot be done
% with just taking the absolute value, because if the robot is producing
% positive torque and the human is producing negative torque, the absolute
% value will not reflect that we wish the robot torque to be negative. 

function C = costFunctionKnee(TorqueH, TorqueR)
    % Set up the weights that will be used for the cost function. In the
    % future this might be an input variable.
    Gt = 1;     %Gain for the torque component
    Ga = 10;     %Gain for the angle between vectors
    
    % Sum the differences in torque between the human and robot
    Ct = 0;
    for i = 1:size(TorqueH, 1)                          % First joint rotation OR rotation about the x axis for multiple DoF joint
        for ii = 1:size(TorqueH, 2)                     % x, y, z torque
            if ii == 1
                Gt = 0.25;
            elseif ii == 2
                Gt = 0.25;
            else 
                Gt = 5;
            end
            for iii = 1:size(TorqueH, 3)                % Second joint rotation OR rotation about the y axis for multiple DoF joint
                for iiii = 1:size(TorqueH, 4)           % Torque about different joints OR rotation about the z axis for multiple DoF joint
                    if TorqueH(i, ii, iii, iiii) >= 0
                        Ct = Ct + Gt*(TorqueR(i, ii, iii, iiii) - TorqueH(i, ii, iii, iiii));
                    else
                        Ct = Ct + Gt*(TorqueH(i, ii, iii, iiii) - TorqueR(i, ii, iii, iiii));
                    end
                end
            end
        end
    end
    
    Ct = Ct/(i*ii*iii*iiii);
    
    % The difference in angle between the human torque vector and the PAM
    % torque vector
    
    Ca = 0;
    for i = 1:size(TorqueH, 1)
        for ii = 1:size(TorqueH, 3)
            for iii = 1:size(TorqueH, 4)
                uvecH = TorqueH(i, :, ii, iii)/norm(TorqueH(i, :, ii, iii));
                
                %Sometimes the BPA can't produce any force due to high
                %contraction. We will set it equal to negative the human
                %vector to maximize the penalty. Consider changing later
                if norm(TorqueR(i, :, ii, iii)) == 0
                    uvecR = -uvecH;
                else
                    uvecR = TorqueR(i, :, ii, iii)/norm(TorqueR(i, :, ii, iii));
                end
                Ca = Ca + dot(uvecH, uvecR) - 1;      %Dot product ranges from 1 (aligned) to -1 (opposite direction). Using the -1, it shifts the range to 0-> -2
            end
        end
    end
    
    Ca = Ca/(i*ii*iii);
    
    % Sum the costs
    C = -Ct - Ga*Ca;
end