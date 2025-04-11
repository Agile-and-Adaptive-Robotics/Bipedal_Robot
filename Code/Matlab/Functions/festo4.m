%% Create a force lookup table based on festo data
function F = festo4(dia, rel, pres)
%Inputs:
%dia == diameter of Festo tube, from Size function
%pres == pressure in kPa
%rel == relative strain
%Outputs:
%F == Force, % of maximum

clear X Y

P = pres/620;           % Normalize pressure
% Fmax20 = 1500;          % Max 20mm BPA force
% Fmax40 = 6000;          % Max 20mm BPA force

    if dia == 10
        load FestoLookup.mat f_10
        F = f_10(rel,P);
    elseif dia == 20
        load FestoLookup.mat f20
       	F = f20(rel,P);
    elseif dia == 40
        load FestoLookup.mat f40
        F = f40(rel,P);
    end
    
%     F(rel > 1) = 0; %No force if shorter than shortest length
%     F(F > 1.05) = NaN; %If force is greater than 5% of it's maximum, return NaN.
%     F = Fn.*maxF;

  end
    
