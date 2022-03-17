function F = festo3(Lmt, rest, dia, long, pres, kmax)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%long == longest musculotendon length
%pres == pressure in kPa
%kmax == maximum measured contraction
%Outputs:
%F == Force, N

load ForceStrainTable.mat RelativeStrain ForceStrain
tendon = long - rest;   %Length of artificial tendon and air fittings

k = (rest-(Lmt-tendon))/rest; %current strain

rel = k/kmax; %relative strain;

pressure = linspace(0,620,20);
column = linspace(1,20,20);
col = interp1(pressure,column,pres); %interpolate force strain column from pressure

for i=1:size(Lmt, 2)
 if rel(i) >= 0 && rel(i) <=1
    F(i) = interp2(ForceStrain, col, rel(i), 'linear',  0);
 else
    F(i) = 0;
 end
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