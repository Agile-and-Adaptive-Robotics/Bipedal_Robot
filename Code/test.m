% This is a script for testing junk related to the bipedal optimization

for i = 1:10
        for ii = 1:15
            if ii == 2*i
                myBreak = 1;
                break
            end
        end
        if myBreak == 1
            break
        end
end

ep = 10;
while ep > 10^-10 || i == 1
    ep = ep*0.1;
end

z = 2+2;