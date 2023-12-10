%% Create a force lookup table based on festo data
function F = festo4(dia, rel, pres)
%Inputs:
%dia == diameter of Festo tube, from Size function
%pres == pressure in kPa
%rel == relative strain
%Outputs:
%F == Force, N

load FestoLookup.mat f40 f20
clear X Y

%rel = contract/0.25;    % Normalize contraction
P = pres/600;           % Normalize pressure
Fmax20 = 1500;          % Max 20mm BPA force
Fmax40 = 6000;          % Max 20mm BPA force

% X = linspace(0,600,7); %Pressure for interpolation
% Y1 = linspace(-0.05,0.25,31);   %Relative strain range for interpolation
% Y2 = linspace(-0.04,0.25,30);

  if dia == 20
            F = f20(rel, P).*Fmax20;
            for i = 1:length(F)
                if F > 1500
                    F = NaN;
                end
            end
  elseif dia == 40
            F = f40(rel, P).*Fmax40;
            for i = 1:length(F)
                if F > 6000
                    F = NaN;
                end
            end
  else

            F = f10(rel,pres);    % This uses the equation fit from Dr. Hunt's lookup table. We now have more accurate ways to calculate force.
        for i = length(F)
            if F(i) < 0
                F(i) = 0;
            end
        end
  end
    
