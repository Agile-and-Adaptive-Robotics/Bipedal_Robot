% P = [];
% Torque = [];
% MomentArm = [];
% Force = Torque./MomentArm;

epsilon_max = 1;                    %Maximum relative contraction
epsilon = linspace(0,1,30)';        %Relative contraction
P = linspace(0,620,19);             %Pressure
B = genlookup(epsilon, epsilon_max, P);

 function A = genlookup(epsilon, epsilon_max, P)
 % parameters
 a6 = 15.6;
 a5 = 1.23;
 a4 = -0.331e-6;
 a3 = -0.461;
 a2 = 2.0265;
 a1 = 192;
 a0 = 254.3;
 S = 0;
 %a = [a0 a1 a2 a3 a4 a5 a6];
 A = NaN(length(epsilon),length(P)+1);         %Make lookup table 
 A(:,1) = epsilon;                           %Make first column of Lookup table epsilon

fI = (1-epsilon)*P;                          %Create a matrix of initial guesses for force

options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);
    for i = 1:length(epsilon)
        for j = 2:length(P)+1
            P1 = P(j-1);
            epsilon1 = epsilon(i);
            f0 = fI(i,j-1);                        %initial guess, Newtons
            g = fsolve(@(F)createLookup(F, P1, epsilon1, epsilon_max, a0, a1, a2, a3, a4, a5, a6, S),f0,options);
            if g < 0
                g = 0;
            end
            A(i,j) = g;
        end
    end

disp(A)

        function Balance = createLookup(F, P1, epsilon1, epsilon_max, a0, a1, a2, a3, a4, a5, a6, S)
        

            Balance = -P1 + a0 + a1*tan( a2*( epsilon1/(a4*F + epsilon_max) + a3 ) ) + a5*F + a6*S;

        end
    

end
