function d = AttachmentDistance(bone, Location)

n = 0;                      %value that will store the number of columns or the number of bone geometric points
for i = 1:size(bone, 2)
    n = n + size(bone{i}, 2);
end


%Find the closest bone point
d = 1000;                               %Value that will be updated with the closest distance from the location
boneIndex = 0;                          %updates to reflect which bonesegment the closest point is in
bonePoint = zeros(3,1);                 %updates with the bone point closest to the attacment point
for i = 1:size(bone, 2)                 %Iterate for every bone segment
    for ii = 1:size(bone{i}, 2)         %Iterate for every point in that bone segment
        for iii = 1:size(Location, 2)   %Iterate for every attachment point location
            check = norm(bone{i}(:, ii) - Location(:, iii));
            if check(1) < check(2)
                boneIndex = i;
                bonePoint
        end
    end
end
