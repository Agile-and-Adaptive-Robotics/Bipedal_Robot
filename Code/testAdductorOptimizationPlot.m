%% X axis torque for adduction/abduction and flexion/extension


% %% Plotting Error
% 
% xTorqueExzRotation = xTorqueRxzRotation - xTorqueHxzRotation;
% 
% yTorqueExzRotation = yTorqueRxzRotation - yTorqueHxzRotation;
% 
% zTorqueExzRotation = zTorqueRxzRotation - zTorqueHxzRotation;

for i = 1:length(xTorqueRxzRotation)
    if xTorqueHxzRotation(i) >= 0
        xTorqueExzRotation(i) = xTorqueRxzRotation(i) - xTorqueHxzRotation(i);
    else
        xTorqueExzRotation(i) = xTorqueHxzRotation(i) - xTorqueRxzRotation(i);
    end
    
    if yTorqueHxzRotation(i) >= 0
        yTorqueExzRotation(i) = yTorqueRxzRotation(i) - yTorqueHxzRotation(i);
    else
        yTorqueExzRotation(i) = yTorqueHxzRotation(i) - yTorqueRxzRotation(i);
    end
    
    if zTorqueHxzRotation(i) >= 0
        zTorqueExzRotation(i) = zTorqueRxzRotation(i) - zTorqueHxzRotation(i);
    else
        zTorqueExzRotation(i) = zTorqueHxzRotation(i) - zTorqueRxzRotation(i);
    end
end

figure
hold on
sgtitle('Adductor Magnus Torque through Abduction/Adduction and Flexion/Extension')
subplot(3, 3, 1)
surf(mTheta*180/pi, mGamma*180/pi, xTorqueHxzRotation, 'EdgeColor', 'none')
title('Muscle X Torque')

subplot(3, 3, 2)
surf(mTheta*180/pi, mGamma*180/pi, xTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM X Torque')

subplot(3, 3, 3)
surf(mTheta*180/pi, mGamma*180/pi, xTorqueExzRotation, 'EdgeColor', 'none')
title('Adjusted Error X Torque')

subplot(3, 3, 4)
surf(mTheta*180/pi, mGamma*180/pi, yTorqueHxzRotation, 'EdgeColor', 'none')
title('Muscle Y Torque')

subplot(3, 3, 5)
surf(mTheta*180/pi, mGamma*180/pi, yTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM Y Torque')

subplot(3, 3, 6)
surf(mTheta*180/pi, mGamma*180/pi, yTorqueExzRotation, 'EdgeColor', 'none')
title('Adjusted Error Y Torque')

subplot(3, 3, 7)
surf(mTheta*180/pi, mGamma*180/pi, zTorqueHxzRotation, 'EdgeColor', 'none')
title('Muscle Z Torque')

subplot(3, 3, 8)
surf(mTheta*180/pi, mGamma*180/pi, zTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM Z Torque')

subplot(3, 3, 9)
surf(mTheta*180/pi, mGamma*180/pi, zTorqueExzRotation, 'EdgeColor', 'none')
title('Adjusted Error Z Torque')

hold off