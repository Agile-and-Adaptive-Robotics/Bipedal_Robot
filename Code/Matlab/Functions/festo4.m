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
    
    for i = length(F)
        if F(i) > 1
            F(i) = NaN;
        elseif F(i)<0
            F(i) = 0;
        else
        end
    end
%     F = Fn.*maxF;

  end
    
