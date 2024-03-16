%% This function finds force of a 10mm BPA at a given pressure, relative strain, and resting length

function F = bpaForce10(restingL, relStrain, Pressure)
rest = restingL;
k = relStrain;
P = Pressure/620;

maxF = maxBPAforce(rest);

load FestoLookup.mat f_10
Fn = f_10(k,P);
F = Fn.*maxF;

        for i = 1:size(F, 1)
            if F(i) < 0
                    F(i) = 0;
            end
            if F(i) > maxF
                    F(i) = NaN;
            end
        end

end