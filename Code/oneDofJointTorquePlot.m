%One degree of freedom joint torque calculation
% Author: Connor Morrow
% Date: 1/12/2021
% Description: This script plots the torque of human muscles, PAMs, and the
% error between the two when actuating about a 1DoF joint

%% X axis torque for adduction/abduction and flexion/extension
xTorqueEzRotation = size(xTorqueHzRotation);
yTorqueEzRotation = size(yTorqueHzRotation);
zTorqueEzRotation = size(zTorqueHzRotation);

for i = 1:length(xTorqueRzRotation)
    if xTorqueHzRotation(i) >= 0
        xTorqueEzRotation(i) = xTorqueRzRotation(i) - xTorqueHzRotation(i);
    else
        xTorqueEzRotation(i) = xTorqueHzRotation(i) - xTorqueRzRotation(i);
    end
    
    if yTorqueHzRotation(i) >= 0
        yTorqueEzRotation(i) = yTorqueRzRotation(i) - yTorqueHzRotation(i);
    else
        yTorqueEzRotation(i) = yTorqueHzRotation(i) - yTorqueRzRotation(i);
    end
    
    if zTorqueHzRotation(i) >= 0
        zTorqueEzRotation(i) = zTorqueRzRotation(i) - zTorqueHzRotation(i);
    else
        zTorqueEzRotation(i) = zTorqueHzRotation(i) - zTorqueRzRotation(i);
    end
end

figure
hold on
sgtitle('Adductor Magnus Torque through Abduction/Adduction and Flexion/Extension')

subplot(3, 2, 1)
plot(phi*180/pi, xTorqueHzRotation, phi*180/pi, xTorqueRzRotation)
title('Muscle and PAM X Torque')
legend('Human', 'PAM')

subplot(3, 2, 2)
plot(phi*180/pi, xTorqueEzRotation)
title('Adjusted Error X Torque')

subplot(3, 2, 3)
plot(phi*180/pi, yTorqueHzRotation, phi*180/pi, yTorqueRzRotation)
title('Muscle and PAM Y Torque')
legend('Human', 'PAM')

subplot(3, 2, 4)
plot(phi*180/pi, yTorqueEzRotation)
title('Adjusted Error Y Torque')

subplot(3, 2, 5)
plot(phi*180/pi, zTorqueHzRotation, phi*180/pi, zTorqueRzRotation)
title('Muscle and PAM Z Torque')
legend('Human', 'PAM')

subplot(3, 2, 6)
plot(phi*180/pi, zTorqueEzRotation)
title('Adjusted Error Z Torque')

hold off