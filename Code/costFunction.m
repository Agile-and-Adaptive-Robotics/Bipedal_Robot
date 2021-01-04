%% PAM Placement Cost Function
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

function C = costFunction(TorqueH, TorqueR)

    % Sum the differences in torque between the human and robot
    C = 0;
    for i = 1:size(TorqueH, 1)                          % First joint rotation OR rotation about the x axis for multiple DoF joint
        for ii = 1:size(TorqueH, 2)                     % x, y, z torque
            for iii = 1:size(TorqueH, 3)                % Second joint rotation OR rotation about the y axis for multiple DoF joint
                for iiii = 1:size(TorqueH, 4)           % Torque about different joints OR rotation about the z axis for multiple DoF joint
                    if TorqueH(i, ii, iii, iiii) >= 0
                        C = C + TorqueR(i, ii, iii, iiii) - TorqueH(i, ii, iii, iiii);
                    else
                        C = C + TorqueH(i, ii, iii, iiii) - TorqueR(i, ii, iii, iiii);
                    end
                end
            end
        end
    end
    
    C = C/(i*ii*iii*iiii);
end