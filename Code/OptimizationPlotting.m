%% -------------------- Plotting ----------------------
%Generate Plots to see how the torques changed after optimization, and
%whether they are better/more in sync with the human model 

%Plot how the cost function changes over the iterations

%Plots for the new robot configuration
figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewRobotTorque1, 'EdgeColor', 'none')
title(RobotTitle1); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewRobotTorque2, 'EdgeColor', 'none')
title(RobotTitle2); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('RobotTorque3', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewRobotTorque3, 'EdgeColor', 'none')
    title(RobotTitle3); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('RobotTorque4', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewRobotTorque4, 'EdgeColor', 'none')
    title(RobotTitle4); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

%Plots for the error from the first robot configuration to the new robot
%configuration
Error1 = HumanTorque1 - RobotTorque1;
Error2 = HumanTorque2 - RobotTorque2;

if exist('HumanTorque3', 'var') == 1
    Error3 = HumanTorque3 - RobotTorque3;
    if exist('HumanTorque4', 'var') == 1
        Error4 = HumanTorque4 - RobotTorque4;
    end
end

NewError1 = HumanTorque1 - NewRobotTorque1;
NewError2 = HumanTorque2 - NewRobotTorque2;

if exist('HumanTorque3', 'var') == 1
    NewError3 = HumanTorque3 - NewRobotTorque3;
    if exist('HumanTorque4', 'var') == 1
        NewError4 = HumanTorque4 - NewRobotTorque4;
    end
end

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error1, 'EdgeColor', 'none')
title(strcat(RobotTitle1, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error2, 'EdgeColor', 'none')
title(strcat(RobotTitle2, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('HumanTorque3', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error3, 'EdgeColor', 'none')
    title(strcat(RobotTitle3, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('HumanTorque4', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, Error4, 'EdgeColor', 'none')
    title(strcat(RobotTitle4, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end


figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError1, 'EdgeColor', 'none')
title(strcat('New ', RobotTitle1, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError2, 'EdgeColor', 'none')
title(strcat('New ', RobotTitle2, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('HumanTorque3', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError3, 'EdgeColor', 'none')
    title(strcat('New ', RobotTitle3, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('HumanTorque4', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, NewError4, 'EdgeColor', 'none')
    title(strcat('New ', RobotTitle4, ' Error')); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end




%% ------ Final plot: Bones and Muscle Locations ----------
%Figure that prints the model of the robot and the muscle

if Muscle1.Diameter == 40
    LW = 4;
elseif Muscle1.Diameter == 20
    LW = 2;
else
    LW = 1;
end

%Original Location
L1 = Muscle1.Location;
L1(:, Cross - 1) = LocationTracker1(:, 1);
L1(:, Cross) = LocationTracker2(:, 1);
figure
hold on
plot3(0, 0, 0, 'o', 'color', 'r')
plot3(Joint1a.Home(1), -Joint1a.Home(2), Joint1a.Home(3), 'o', 'color', 'r')
plot3(Spine(1, :), -Spine(3, :), Spine(2, :), '.', 'color', 'b')
plot3(Sacrum(1, :), -Sacrum(3, :), Sacrum(2, :), '.', 'color', 'b')
plot3(Pelvis(1, :), -Pelvis(3, :), Pelvis(2, :), '.', 'color', 'b')
plot3(L1(1, :), -L1(3, :), L1(2, :), 'color', 'r', 'LineWidth', LW)
axis(axisLimits)
title('Original Muscle Location')
hold off

%New best location
L2 = Muscle1.Location;
L2(:, Cross - 1) = StartingLocation1;
L2(:, Cross) = StartingLocation2;
figure
hold on
plot3(0, 0, 0, 'o', 'color', 'r')
plot3(Joint1a.Home(1), -Joint1a.Home(2), Joint1a.Home(3), 'o', 'color', 'r')
plot3(Spine(1, :), -Spine(3, :), Spine(2, :), '.', 'color', 'b')
plot3(Sacrum(1, :), -Sacrum(3, :), Sacrum(2, :), '.', 'color', 'b')
plot3(Pelvis(1, :), -Pelvis(3, :), Pelvis(2, :), '.', 'color', 'b')
plot3(L2(1, :), -L2(3, :), L2(2, :), 'color', 'r', 'LineWidth', LW)
axis(axisLimits)
title('Current Optimization Solution for Muscle Location')
hold off