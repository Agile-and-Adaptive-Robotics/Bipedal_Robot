%%Hunt Equation is Dr. Hunt's equation for the Pressure-force-length
%%relationship of Festo 10mm BPAs
% Input:
%  epsilon = contraction
%  epsilon_max = maximum contraction
%  P = Pressure
%  S = Inflating, Deflating
% Output:
%  F = Force
 
 function A = HuntEq(epsilon, epsilon_max, P, S)
 % parameters
 a6 = 15.5970188831175;
 a5 = 1.22902680741856;
 a4 = -0.000331227432105755;
 %a3 = -0.4609;
 a3 = -0.460864147695817;
 a2 = 2.02665622940754;
 a1 = 191.919433464452;
 a0 = 254.327839024804;
 %S = 0;
 %a = [a0 a1 a2 a3 a4 a5 a6];
 A = NaN(length(epsilon),length(P)+1);         %Make lookup table 
 A(:,1) = epsilon/epsilon_max;                 %Make first column of Lookup table relative strain

fI = (1-epsilon-epsilon_max)*P;                          %Create a matrix of initial guesses for force


%options = optimoptions('fsolve','Display','none','FunctionTolerance',0.0001);
options = optimset(@fzero);
%options.PlotFcns = @optimplotfval;
options.FunValCheck = 'on';

    for i = 1:length(epsilon)
        for j = 2:length(P)+1
            P1 = P(j-1);
            epsilon1 = epsilon(i);
  %          f0 = fI(i,j-1)-0.1*P1;                        %initial guess, Newtons
            f0 = 1;
            [g fval exitflag output] = fzero(@(F)createLookup(F, P1, epsilon1, epsilon_max, a0, a1, a2, a3, a4, a5, a6, S),f0,options)
            if g < 0
                g = 0;
            end
            A(i,j) = g;
   
        end
    end

% disp(A/1000)

        function Balance = createLookup(F, P1, epsilon1, epsilon_max, a0, a1, a2, a3, a4, a5, a6, S)
        

            Balance = -P1 + a0 + a1.*tan( a2.*( epsilon1./(a4.*F + epsilon_max) + a3 ) ) + a5.*F + a6*S;

        end
    

end
