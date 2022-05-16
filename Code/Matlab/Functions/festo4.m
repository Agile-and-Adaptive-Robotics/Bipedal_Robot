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

X = linspace(0,600,7); %Pressure for interpolation
Y1 = linspace(-0.05,0.25,31);   %Relative strain range for interpolation
Y2 = linspace(-0.04,0.25,30);

  if dia == 20
            F = f20(contract, pres);
            if F > 1500
                F = 1500;
            end
  elseif dia == 40
            F = f40(contract, pres);
            if F > 6000
                F = 6000;
            end
  else
            x = [0,   0.1, 0.17, 0.224]';
            z = [630, 300, 150,  0]';
            BPAFit = fit(x, z, 'linearinterp');
            F = BPAFit(contract);
  end

    if F < 0
        F = 0;
    end