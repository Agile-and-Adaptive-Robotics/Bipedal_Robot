%% -------------------- Plotting ----------------------
%Generate Plots to see how the torques changed after optimization, and
%whether they are better/more in sync with the human model 

%Plot how the cost function changes over the iterations

%Plots for the new robot configuration
figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, OriginalRobotTorque1, 'EdgeColor', 'none')
title(RobotTitle1); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

figure
surf(RobotAxis1*180/pi, RobotAxis2*180/pi, OriginalRobotTorque2, 'EdgeColor', 'none')
title(RobotTitle2); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
zlabel('Torque, N*m')
colorbar; caxis(caxisRange)

if exist('RobotTorque3', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, OriginalRobotTorque3, 'EdgeColor', 'none')
    title(RobotTitle3); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

if exist('RobotTorque4', 'var') == 1
    figure
    surf(RobotAxis1*180/pi, RobotAxis2*180/pi, OriginalRobotTorque4, 'EdgeColor', 'none')
    title(RobotTitle4); xlabel(RobotAxis1Label); ylabel(RobotAxis2Label)
    zlabel('Torque, N*m')
    colorbar; caxis(caxisRange)
end

%Plots for the error from the first robot configuration to the new robot
%configuration
Error1 = HumanTorque1 - OriginalRobotTorque1;
Error2 = HumanTorque2 - OriginalRobotTorque2;

if exist('HumanTorque3', 'var') == 1
    Error3 = HumanTorque3 - OriginalRobotTorque3;
    if exist('HumanTorque4', 'var') == 1
        Error4 = HumanTorque4 - OriginalRobotTorque4;
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
clear L
for m = 1:MuscleNum
    L{m} = Muscles{m}.Location;
    L{m}(:, Muscles{m}.CrossPoints - 1) = LocationTracker1{m}(:, 1);
    L{m}(:, Muscles{m}.CrossPoints) = LocationTracker2{m}(:, 1);
end

figure
hold on
plot3(0, 0, 0, 'o', 'color', 'r')
plot3(Joint1a.Home(1), -Joint1a.Home(2), Joint1a.Home(3), 'o', 'color', 'r')
plot3(Spine(1, :), -Spine(3, :), Spine(2, :), '.', 'color', 'b')
plot3(Sacrum(1, :), -Sacrum(3, :), Sacrum(2, :), '.', 'color', 'b')
plot3(Pelvis(1, :), -Pelvis(3, :), Pelvis(2, :), '.', 'color', 'b')
for m = 1:MuscleNum
    plot3(L{m}(1, :), -L{m}(3, :), L{m}(2, :), 'color', 'r', 'LineWidth', LW)
end
axis(axisLimits)
title('Original Muscle Location')
hold off

[a, b] = min(C);

%New best location
for m = 1:MuscleNum
    Lb{m}(:, :) = Muscles{m}.Location;
    Lb{m}(:, Muscles{m}.CrossPoints - 1) = LocationTracker1{m}(:, b);
    Lb{m}(:, Muscles{m}.CrossPoints) = LocationTracker2{m}(:, b);
end
figure
hold on
plot3(0, 0, 0, 'o', 'color', 'r')
plot3(Joint1a.Home(1), -Joint1a.Home(2), Joint1a.Home(3), 'o', 'color', 'r')
plot3(Spine(1, :), -Spine(3, :), Spine(2, :), '.', 'color', 'b')
plot3(Sacrum(1, :), -Sacrum(3, :), Sacrum(2, :), '.', 'color', 'b')
plot3(Pelvis(1, :), -Pelvis(3, :), Pelvis(2, :), '.', 'color', 'b')
for m = 1:MuscleNum
    plot3(Lb{m}(1, :), -Lb{m}(3, :), Lb{m}(2, :), 'color', 'r', 'LineWidth', LW)
end
axis(axisLimits)
title('Current Optimization Solution for Muscle Location')
hold off

figure
hold on
plot3(0, 0, 0, 'o', 'color', 'r')
plot3(Joint1a.Home(1), -Joint1a.Home(2), Joint1a.Home(3), 'o', 'color', 'r')
plot3(Spine(1, :), -Spine(3, :), Spine(2, :), '.', 'color', 'b')
plot3(Sacrum(1, :), -Sacrum(3, :), Sacrum(2, :), '.', 'color', 'b')
plot3(Pelvis(1, :), -Pelvis(3, :), Pelvis(2, :), '.', 'color', 'b')
for m = 1:MuscleNum
    plot3(L{m}(1, :), -L{m}(3, :), L{m}(2, :), 'color', 'r', 'LineWidth', LW)
end
for m = 1:MuscleNum
    plot3(Lb{m}(1, :), -Lb{m}(3, :), Lb{m}(2, :), 'color', 'g', 'LineWidth', LW)
end
axis(axisLimits)
title('Current Optimization Solution for Muscle Location')
hold off