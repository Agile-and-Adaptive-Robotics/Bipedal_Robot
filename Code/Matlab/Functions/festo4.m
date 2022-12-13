%% Create a force lookup table based on festo data
function F = festo4(dia, pres, contract)
%Inputs:
%dia == diameter of Festo tube, from Size function
%pres == pressure in kPa
%contract == percent of original length contracted, e.g. a 400 mm BPA that
%is inflacted to a length of 300 mm would be 25% (contract = .25)
%Outputs:
%F == Force, N

load FestoLookup.mat f40 f20
clear X Y

rel = contract/0.25;    % Normalize contraction
P = pres/600;           % Normalize pressure
Fmax20 = 1500;          % Max 20mm BPA force
Fmax40 = 6000;          % Max 20mm BPA force

% X = linspace(0,600,7); %Pressure for interpolation
% Y1 = linspace(-0.05,0.25,31);   %Relative strain range for interpolation
% Y2 = linspace(-0.04,0.25,30);

  if dia == 20
            F = f20(rel, P)*Fmax20;
            if F > 1500
                F = 1500;
            end
  elseif dia == 40
            F = f40(rel,P)*Fmax40;
            if F > 6000
                F = 6000;
            end
  else
            x1 = [ -.03   -.02  -.01      0]';
            z1 = [741.9    613   523  458.2]';
            z2 = [759.2  629.3 539.3  473.3]';
            z3 = [785.1  653.5 562.6  495.9]';
            x = [x1; x1; x1];
            y = [580*ones(length(x1),1);
                 600*ones(length(x1),1);
                 630*ones(length(x1),1)];
            z = [z1; z2; z3]; 
            BPAFit = fit([x, y],z,'linearinterp','Normalize','on');
            F = BPAFit(contract,pres);
  end

    if F < 0
        F = 0;
    end