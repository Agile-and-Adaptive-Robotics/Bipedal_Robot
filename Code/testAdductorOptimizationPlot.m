%% X axis torque for adduction/abduction and flexion/extension
figure
surf(mTheta*180/pi, mGamma*180/pi, xTorqueHxzRotation, 'EdgeColor', 'none')
title('AddMag X Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

figure
surf(mTheta*180/pi, mGamma*180/pi, xTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM AddMag X Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

%% Y axis torque for adduction/abduction and flexion/extension
figure
surf(mTheta*180/pi, mGamma*180/pi, yTorqueHxzRotation, 'EdgeColor', 'none')
title('AddMag Y Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

figure
surf(mTheta*180/pi, mGamma*180/pi, yTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM AddMag Y Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

%% Z axis torque for adduction/abduction and flexion/extension
figure
surf(mTheta*180/pi, mGamma*180/pi, zTorqueHxzRotation, 'EdgeColor', 'none')
title('AddMag Z Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

figure
surf(mTheta*180/pi, mGamma*180/pi, zTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM AddMag Z Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

%% Plotting Error

xTorqueExzRotation = xTorqueRxzRotation - xTorqueHxzRotation;

yTorqueExzRotation = yTorqueRxzRotation - yTorqueHxzRotation;

zTorqueExzRotation = zTorqueRxzRotation - zTorqueHxzRotation;

figure
surf(mTheta*180/pi, mGamma*180/pi, xTorqueExzRotation, 'EdgeColor', 'none')
title('Error AddMag X Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

figure
surf(mTheta*180/pi, mGamma*180/pi, yTorqueExzRotation, 'EdgeColor', 'none')
title('Error AddMag Y Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

figure
surf(mTheta*180/pi, mGamma*180/pi, zTorqueExzRotation, 'EdgeColor', 'none')
title('Error AddMag Z Tor during Abd/Add and Flex/Ext')
xlabel('Abduction/Adduction, deg')
ylabel('Flexion/Extension, deg')
zlabel('Torque, N*m')
colorbar
caxis([-40, 150])

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
title('Error X Torque')

subplot(3, 3, 4)
surf(mTheta*180/pi, mGamma*180/pi, yTorqueHxzRotation, 'EdgeColor', 'none')
title('Muscle Y Torque')

subplot(3, 3, 5)
surf(mTheta*180/pi, mGamma*180/pi, yTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM Y Torque')

subplot(3, 3, 6)
surf(mTheta*180/pi, mGamma*180/pi, yTorqueExzRotation, 'EdgeColor', 'none')
title('Error Y Torque')

subplot(3, 3, 7)
surf(mTheta*180/pi, mGamma*180/pi, zTorqueHxzRotation, 'EdgeColor', 'none')
title('Muscle Z Torque')

subplot(3, 3, 8)
surf(mTheta*180/pi, mGamma*180/pi, zTorqueRxzRotation, 'EdgeColor', 'none')
title('PAM Z Torque')

subplot(3, 3, 9)
surf(mTheta*180/pi, mGamma*180/pi, zTorqueExzRotation, 'EdgeColor', 'none')
title('Error Z Torque')

hold off