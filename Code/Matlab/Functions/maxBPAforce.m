%% Equation for maximum force in 10, 20, and 40 mm BPAs given a resting length
%%Ben has five 10 mm ID BPAs with force measurements at maximum pressure and
%many relative strains. He also has many BPAs measured at zero strain and
%various pressures. 
%Lawrence has data on five BPAs with high resolution of pressure data but
%only at a few different amounts of relative strain.
%Mo tested many 20 mm BPAs at zero strain.
%40 mm BPA limit is specified from Festo tool.

%Inputs:
%   - Resting length (in meters)
%   - Diameter (as a string?): if no diameter given, use 10 mm;
%   - Pressure (in kPa); if no pressure given, use 620;
%Outputs:
%   - Force (Newtons)

function z = maxBPAforce(restingLength, varargin)
p = inputParser;

checkLength = @(x) (x > 0) && isnumeric(x);
addRequired(p,'restingLength',checkLength);

defaultDia = '10';
validDia = {'10','20','40'};
checkDia = @(x) any(validatestring(x,validDia));
addOptional(p,'diameter',defaultDia,checkDia);

defaultPres = 620;
checkPres = @(y) (y > 0) && isnumeric(y);
addOptional(p,'pressure',defaultPres,checkPres);
parse(p, restingLength, varargin{:})

x = restingLength;
y = p.Results.pressure;
dia = p.Results.diameter;

    switch dia
        case '10'
        %Parameters for F_{max10} equation
        a1 = 0.4895;                    %N\per\kPa
        a2 = 0.03068;                   %\per\kPa\per\m
        % a = [a1 a2];
        % a_ub = [];      %Upper bound on coefficient confidence interval
        % a_lb = [];      %Lower bound on coefficient confidence interval                   
        z = y.*(a1.*atan(a2.*(x-0.0075).*y));
        
        case '20'
        %Parameters for F_{max20} equation
        a1 = 1.4877;                    %N\per\kPa
        a2 = 0.0248;                   %\per\kPa\per\m
        % a = [a1 a2];
        % a_ub = [];      %Upper bound on coefficient confidence interval
        % a_lb = [];      %Lower bound on coefficient confidence interval                   
        z = y.*(a1.*atan(a2.*(x-0.0075).*y));
        
        case '40'
        %From Festo tool. Irrespective of length.
        z = 6398.4;
        
        otherwise
            error('Error in maxBPAforce function');
    end
        
        
    
end