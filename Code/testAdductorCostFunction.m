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