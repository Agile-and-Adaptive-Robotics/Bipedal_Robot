%% Calculating the Cost of the torque error
% The first part of the cost function will look at how the torque between
% the human muscles and the PAM placements differ. We want the robot
% muscles to be able to develop equal to or greater than the amount of
% torque that the human muscles can generate. 

% For the first test, let's look specifically at the x axis torque

PosxTorqueHxzRotation = abs(xTorqueHxzRotation);
PosxTorqueRxzRotation = abs(xTorqueRxzRotation);

PosxTorqueExzRotation = PosxTorqueRxzRotation - PosxTorqueHxzRotation;

figure
surf(mTheta*180/pi, mGamma*180/pi, PosxTorqueHxzRotation, 'EdgeColor', 'none')
% surf(mTheta*180/pi, mGamma*180/pi, PosxTorqueExzRotation, 'EdgeColor', 'none')
title('AddMag X Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

%% Calculate the total scalar error between the human and robot torque

eC(1) = 0;
for i = 1:size(TorqueH, 1)
    for ii = 1:size(TorqueH, 2)
        for iii = 1:size(TorqueH, 3)
            for iiii = 1:size(TorqueH, 4)
                if TorqueH(i, ii, iii, iiii) >= 0
                    eC(1) = eC(1) + TorqueR(i, ii, iii, iiii) - TorqueH(i, ii, iii, iiii);
                else
                    eC(1) = eC(1) + TorqueH(i, ii, iii, iiii) - TorqueR(i, ii, iii, iiii);
                end
            end
        end
    end
end

% Normalize the error to the amount of data points that it was summed over
eC(1) = eC(1)/(i*ii*iii*iiii);
