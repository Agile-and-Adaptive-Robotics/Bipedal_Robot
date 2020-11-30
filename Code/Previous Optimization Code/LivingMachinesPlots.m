%% Living Machines Plots
close all



figure
hold on
plot(C, 'k')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')
xlim([0, length(C)])
ylim([0, 800])
set(gca, 'FontSize', 12)
hold off

%Individual Cost Components

figure
subplot(2, 2, 1)
plot(eC, 'k')
title('Error Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')
subplot(2, 2, 2)
plot(mC, 'k')
title('Muscle Length Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')
subplot(2, 2, 3)
plot(disC, 'k')
title('Distance Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')
subplot(2, 2, 4)
plot(dC, 'k')
title('Diameter Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')

figure
subplot(2, 4, [1 2 5 6])
plot(C, 'k')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')
xlim([0, length(C)])
ylim([0, 800])
set(gca, 'FontSize', 12)

subplot(2, 4, 3)
plot(eC, 'k')
title('Error Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')

subplot(2, 4, 4)
plot(mC, 'k')
title('Muscle Length Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')

subplot(2, 4, 7)
plot(disC, 'k')
title('Distance Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')

subplot(2, 4, 8)
plot(dC, 'k')
title('Diameter Component')
xlabel('Iterations', 'FontWeight', 'Bold')
ylabel('Cost Value', 'FontWeight', 'Bold')


%Error plots
% figure
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error1, 'EdgeColor', 'none')
% title(strcat(RobotTitle1, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% zlabel('Torque, N*m')
% colorbar; caxis(caxisRange)
% 
% figure
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error2, 'EdgeColor', 'none')
% title(strcat(RobotTitle2, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% zlabel('Torque, N*m')
% colorbar; caxis(caxisRange)
% 
% figure
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError1, 'EdgeColor', 'none')
% title(strcat('New ', RobotTitle1, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% zlabel('Torque, N*m')
% colorbar; caxis(caxisRange)
% 
% figure
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError2, 'EdgeColor', 'none')
% title(strcat('New ', RobotTitle2, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% zlabel('Torque, N*m')
% colorbar; caxis(caxisRange)

caxisRange = [-200 200];
RobotAxis1Label = ['Flexion/Extension,' char(176)];
RobotAxis2Label = ['Lateral Bending,' char(176)];

% figure
% subplot(2, 3, 1)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error1, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('Previous Model X Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 3, 2)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error2, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('Previous Model Y Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 3, 3)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error3, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('Previous Model Z Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 3, 4)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError1, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('New Model X Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 3, 5)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError2, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('New Model Y Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 3, 6)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError3, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('New Model Z Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 



% figure
% subplot(2, 2, 1)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error1, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('Previous Model X Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 2, 2)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error3, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('Previous Model Z Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 2, 3)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError1, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('New Model X Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
% 
% subplot(2, 2, 4)
% surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError3, 'EdgeColor', 'none')
% zlabel('Torque, N*m')
% caxis(caxisRange)
% zlim([-200 200])
% title('New Model Z Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

figure
subplot(2, 2, 1)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, abs(Error1), 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('Previous Model X Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

subplot(2, 2, 2)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, abs(Error3), 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('Previous Model Z Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

subplot(2, 2, 3)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, abs(NewError1), 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('New Model X Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

subplot(2, 2, 4)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, abs(NewError3), 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('New Model Z Axis Torque Difference'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

%Figure showing the human torques and the robot torques
figure
subplot(2, 2, 1)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, HumanTorque1, 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('Human X Axis Torque'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

subplot(2, 2, 2)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, HumanTorque3, 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('Human Z Axis Torque'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

subplot(2, 2, 3)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewRobotTorque1, 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('New Model X Axis Torque'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

subplot(2, 2, 4)
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewRobotTorque3, 'EdgeColor', 'none')
zlabel('Torque, N*m')
caxis(caxisRange)
zlim([0 250])
title('New Model Z Axis Torque'); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)

