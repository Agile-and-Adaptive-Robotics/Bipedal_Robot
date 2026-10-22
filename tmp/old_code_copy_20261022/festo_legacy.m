function F = festo_legacy(Lmt, rest, dia, long)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, from Size function
%dia == diameter of Festo tube, from Size function
%Outputs:
%F == Force, N

A0 = 254.3; %kPa
A1 = 192.0; %kPa
A2 = 2.0265;
A3 = -0.461;
A4 = -0.331; %1/mN
A5 = 1.230; %kPa/mN
A6 = 15.6; %kPa
S = 0;
P = 620; %Pressure, kPa
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

    function ForceBal = myfunc(F)
        
        ForceBal = A0 + A1 *tan(A2*(k/(A4*F+kmax)+A3)) + A5*F + A6*S - P;
             
    end 

x0 = 250; %Initial guess, N
options = optimoptions('fsolve','FunctionTolerance',0.01);
F = fsolve(@myfunc,x0,options); %Force, N

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