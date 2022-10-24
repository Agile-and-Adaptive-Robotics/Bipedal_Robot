%% This function finds force of a 10mm BPA at a given pressure, relative strain, and resting length

function F = bpaForce10(restingL, relStrain, Pressure)

rest = restingL;
k = relStrain;
P = Pressure;

maxForce = maxBPAforce(rest);

Fn = cell(1,4);
[Fn{1}, Fn{2}, Fn{3}, Fn{4}] = normF10(k, P);


F = zeros([size(Fn{1}),length(Fn)]);

for i = 1:length(Fn)
            A = Fn{i};
            F(:,:,i) = maxForce.*A;
end

for i = 1:length(Fn)
    for r = 1:size(Fn{1},1)
        for c = 1:size(Fn{1},2)
            if F(r,c,i) < 0
                F(r,c,i) = NaN;
            elseif F(r,c,i) > 1.1*maxForce
                F(r,c,i) = NaN;
            else
            end
        end
    end
end

end