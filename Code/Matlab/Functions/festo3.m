function F = festo3(Lmt, rest, dia, pres, kmax)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%long == longest musculotendon length
%pres == pressure in kPa
%kmax == maximum measured contraction length. Input as length, will
%convert to percent in the code
%Outputs:
%F == Force, N

load ForceStrainTable.mat ForceStrain
tendon = 0;   %Length of artificial tendon and air fittings
fitting = 0.0254; %End cap length
            
kmax = (rest-kmax)/rest; %Convert maximum contraction from length into percent

X = linspace(0,620,20); %Pressure for interpolation
Y = linspace(0,1,30);   %Relative strain range for interpolation
           
k = zeros(size(Lmt,1),1);
rel = zeros(size(Lmt,1),1);
F = zeros(size(Lmt,1),1);
            for i = 1:size(Lmt, 1)
                k(i,1) = (rest-(Lmt(i,1)))/rest; %current strain
                rel(i,1) = k(i,1)/kmax; %relative strain
                if rel(i,1) >= 0 && rel(i,1) <=1
                    F(i,1) = interp2(X, Y, ForceStrain, pres, rel(i), 'linear',  0);
                else
                    F(i,1) = 0;
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