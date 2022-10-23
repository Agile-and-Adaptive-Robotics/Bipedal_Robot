%% This function finds force of a 10mm BPA at a given pressure, relative strain, and resting length

function F = bpaForce10(restingL, relStrain, Pressure)

rest = restingL;
k = relStrain;
P = Pressure;

maxForce = maxBPAforce(rest);

Fn = normF10(k, P);

for i = 1:length(Fn)
    F{i} = maxForce.*Fn{i};
end

end