% P = [];
% Torque = [];
% MomentArm = [];
% Force = Torque./MomentArm;

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

epsilon_max = 1;                    %Maximum relative contraction
epsilon = linspace(0,1,30)';        %Relative contraction
P = linspace(0,620,20);             %Pressure

A = NaN(length(epsilon),length(P));         %Make lookup table 
A(:,1) = epsilon;                           %Make first column of Lookup table epsilon
options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);

for i = 1:length(epsilon)
    for j = 2:length(P)
        f0 = 10;                %initial guess, Newtons
        P1 = P(j);
        epsilon1 = epsilon(i);
        vars = [P1, epsilon1];
        g = fsolve(@createLookup,f0,options);
        if g < 0
            g = 0;
        end
        A(i,j) = g;
    end
end

disp(A)

    function Balance = createLookup(F)
        

        Balance = P1 - a0 + a1*tan( a2*( epsilon1/(a4*F + epsilon_max) + a3 ) ) + a5*F + a6*S;

    end
    


    