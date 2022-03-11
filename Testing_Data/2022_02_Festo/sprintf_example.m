tests = 4;
runsperseries = [10 15 30 12];

for i = 1:tests
    for j = 1:runsperseries(i)
                file_name = sprintf('FlxTest%0.0f_%0.0f.mat', j,i)
    end
end