%% This script moves the point 
Location = Bifemsh_Pam.Location;

[val, pos] = min(abs(phi));

u = Bifemsh_Pam.UnitDirection(pos, :);

move = 0.05*u

newL1 = Location(1, :) - move

newL2 = Location(2, :) + move

newL = [newL1; newL2];

if isequal(Bones{1}, 'Pelvis')
    offset1 = 0;
elseif isequal(Bones{1}, 'Femur')
    offset1 = Hip;
end

if isequal(Bones{2}, 'Femur')
    offset2 = Hip;
elseif isequal(Bones{2}, 'Tibia')
    offset2 = Knee;
end


for i = 1:size(RMuscleLocation, 2)
    newL(1:RMuscleCross{i}-1, :) = newL(1:RMuscleCross{i}-1, :) + offset1;
    newL(RMuscleCross{i}:end, :) = newL(RMuscleCross{i}:end, :) + offset1 + offset2;
    newL = newL*RotationM;
end


figure
hold on
%Bone Plotting
plot3(Spine(:, 1), Spine(:, 2), Spine(:, 3), '.', 'color', 'b');
plot3(Sacrum(:, 1), Sacrum(:, 2), Sacrum(:, 3), '.', 'color', 'b');
plot3(Pelvis(:, 1), Pelvis(:, 2), Pelvis(:, 3), '.', 'color', 'b');
plot3(Femur(:, 1), Femur(:, 2), Femur(:, 3), '.', 'color', 'b');
plot3(Tibia(:, 1), Tibia(:, 2), Tibia(:, 3), '.', 'color', 'b');
plot3(Talus(:, 1), Talus(:, 2), Talus(:, 3), '.', 'color', 'b');
plot3(Calcaneus(:, 1), Calcaneus(:, 2), Calcaneus(:, 3), '.', 'color', 'b');
plot3(Toes(:, 1), Toes(:, 2), Toes(:, 3), '.', 'color', 'b');

%MusclePlotting
for i = 1:size(HMuscleLocation, 2)
    plot3(HMuscleLocation{i}(:, 1), HMuscleLocation{i}(:, 2), HMuscleLocation{i}(:, 3), '.-', 'color', 'r', 'linewidth', 3);
end
for i = 1:size(RMuscleLocation, 2)
        plot3(RMuscleLocation{i}(:, 1), RMuscleLocation{i}(:, 2), RMuscleLocation{i}(:, 3), '.-', 'color', 'g', 'linewidth', 3);
end
plot3(newL(:, 1), newL(:, 2), newL(:, 3), '.-', 'color', 'k', 'linewidth', 2)
axis(axisLimits)
hold off

A