%% Create a force lookup table based on festo data
function F = festo4(dia, pres, contract)
%Inputs:
%dia == diameter of Festo tube, from Size function
%pres == pressure in kPa
%contract == percent of original length contracted, e.g. a 400 mm BPA that
%is inflacted to a length of 300 mm would be 25% (contract = .25)
%Outputs:
%F == Force, N

load ForceStrainTable.mat FestoLookup40 FestoLookup20
clear X Y

X = linspace(0,600,7); %Pressure for interpolation
Y = linspace(-0.05,0.25,31);   %Relative strain range for interpolation


  if dia == 20
            F = interp2(Y, X, FestoLookup20, contract, pres, 'linear',  0);
  elseif dia == 40
            F = interp2(Y, X, FestoLookup40, contract, pres, 'linear',  0);
  else
            x = [0, 0.1, 0.17, 0.25]';
            y = [630, 300, 150, 0]';
            BPAFit = fit(x, y, 'linearinterp');
            F = BPAFit(contract);
  end