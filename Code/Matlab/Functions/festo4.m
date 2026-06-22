%% Create a force lookup table based on festo data
function F = festo4(dia, rel, pres)
%Inputs:
%dia == diameter of Festo tube, from Size function
%pres == pressure in kPa
%rel == relative strain
%Outputs:
%F == Force, % of maximum

persistent f_10 f20 f40

P = pres/620;           % Normalize pressure
% Fmax20 = 1500;          % Max 20mm BPA force
% Fmax40 = 6000;          % Max 20mm BPA force

     if isempty(f20)
        S = load('FestoLookup.mat', 'f_10', 'f20', 'f40');
        f_10 = S.f_10;
        f20  = S.f20;
        f40  = S.f40;
    end
    
    switch dia
        case 10
            F = f_10(rel, P);
        case 20
            F = f20(rel, P);
        case 40
            F = f40(rel, P);
        otherwise
            error('Unsupported BPA diameter: %g', dia);
    end

    F(rel > 1) = 0; %No force if shorter than shortest length
%     F(F > 1.05) = NaN; %If force is greater than 5% of it's maximum, return NaN.
%     F = Fn.*maxF;

  end
    
