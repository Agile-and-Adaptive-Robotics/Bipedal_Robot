function F = festo(Lmt, rest, dia, long)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%Outputs:
%F == Force, N

load ForceStrainTable.mat RelativeStrain Force
tendon = long - rest;   %Length of artificial tendon and air fittings

k = (rest-(Lmt-tendon))/rest; %current strain
act = [15.31; 18.28; 18.94; 19.50; 27.27; 28.09]; %Resting actuator lengths (Hunt 2017)
strain = [0.1491; 0.1618; 0.1633; 0.1680; 0.1692; 0.1750]; %Max strain for these lengths (Hunt 2017)

if rest >= max(act)
    kmax = max(strain);                 %maximum strain
elseif rest <= min(act)
    kmax = min(strain);                 %maximum strain
else
    kmax = interp1q(act,strain,rest);   %maximum strain
end

rel = k/kmax; %relative strain;

if rel >= 0 && rel <=1
    F = interp1(RelativeStrain, Force, rel);
else
    F = 0;
end

%If diameter is not 10 mm, then upscale force
if dia == 20
    F = (1500/630)*F;
end

if dia == 40
    F = (6000/630)*F;
end

end

%Reference:
%Hunt, Alexander J., Alexander Graber-Tilton, and Roger D. Quinn. "Modeling length effects of braided pneumatic actuators."
%In ASME 2017 International Design Engineering Technical Conferences and Computers and Information in Engineering Conference, pp. V05AT08A008-V05AT08A008.
%American Society of Mechanical Engineers, 2017.