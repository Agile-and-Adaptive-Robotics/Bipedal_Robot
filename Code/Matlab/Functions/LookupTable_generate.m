row = 30;
col = 19;

epsilon_max = 0.18;                    %Maximum relative contraction
epsilon = linspace(0,epsilon_max,row)';%Relative contraction
P = linspace(0,620,col);             %Pressure

%% Isometric Condition
S0 = 0;                             %Not inflating or deflating
B0 = HuntEq(epsilon, epsilon_max, P, S0);

%Only keep first zero value in columns. Set others to NaN. This helps with
%curve fitting.
A0 = B0(:,2:col+1);
for i = 2:row
    for j = 1:col
        if isnan(A0(i-1,j)) == 1
            A0(i,j) = NaN;
        elseif A0(i-1,j) == 0
            A0(i,j) = NaN;
        else 
        end
    end
end

%% Inflating condition 
S1 = 1;                             %Inflating
B1 = HuntEq(epsilon, epsilon_max, P, S1);
A1 = B1(:,2:col+1);
for i = 2:row
    for j = 1:col
        if isnan(A1(i-1,j)) == 1
            A1(i,j) = NaN;
        elseif A1(i-1,j) == 0
            A1(i,j) = NaN;
        else 
        end
    end
end

%% Deflating condition
S_1 = -1;                           %Deflating
B_1 = HuntEq(epsilon, epsilon_max, P, S_1);
A_1 = B_1(:,2:col+1);
for i = 2:row
    for j = 1:col
        if isnan(A_1(i-1,j)) == 1
            A_1(i,j) = NaN;
        elseif A_1(i-1,j) == 0
            A_1(i,j) = NaN;
        else 
        end
    end
end

%% This is now in a separate m-file called HuntEq
%  function A = genlookup(epsilon, epsilon_max, P)
%  % parameters
%  a6 = 15.6;
%  a5 = 1.23;
%  a4 = -0.331e-6;
%  a3 = -0.461;
%  a2 = 2.0265;
%  a1 = 192;
%  a0 = 254.3;
%  S = 0;
%  %a = [a0 a1 a2 a3 a4 a5 a6];
%  A = NaN(length(epsilon),length(P)+1);         %Make lookup table 
%  A(:,1) = epsilon;                           %Make first column of Lookup table epsilon
% 
% fI = (1-epsilon)*P;                          %Create a matrix of initial guesses for force
% 
% 
% options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);
%     for i = 1:length(epsilon)
%         for j = 2:length(P)+1
%             P1 = P(j-1);
%             epsilon1 = epsilon(i);
%             f0 = fI(i,j-1);                        %initial guess, Newtons
%             g = fsolve(@(F)createLookup(F, P1, epsilon1, epsilon_max, a0, a1, a2, a3, a4, a5, a6, S),f0,options);
%             if g < 0
%                 g = 0;
%             end
%             A(i,j) = g;
%         end
%     end
% 
% disp(A)
% 
%         function Balance = createLookup(F, P1, epsilon1, epsilon_max, a0, a1, a2, a3, a4, a5, a6, S)
%         
% 
%             Balance = -P1 + a0 + a1*tan( a2*( epsilon1/(a4*F + epsilon_max) + a3 ) ) + a5*F + a6*S;
% 
%         end
%     
% 
% end
