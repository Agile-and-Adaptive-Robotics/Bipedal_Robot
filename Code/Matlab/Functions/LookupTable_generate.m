row = 30;
col = 19;

epsilon_max = 0.171;                    %Maximum relative contraction
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

%% Normalize all the pressure, strain, and force
%    - scaling everything helps with curve fitting
%    - Fmax seems to be a function of resting length

rel = epsilon./epsilon_max;   %relative strain
Prel = P./max(P);             %relative pressure
Fmax = max(A_1(1,:));         %maximum force

%create force-strain-pressure tables scaled to maximum force
z0 = A0./Fmax;
z1 = A1./Fmax;
z_1 = A_1./Fmax;




