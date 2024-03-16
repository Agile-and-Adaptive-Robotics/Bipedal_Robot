function F = festo3(Lmt, rest, dia, pres, kmax, ten, Fitting)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%long == longest musculotendon length
%pres == pressure in kPa
%kmax == maximum measured contraction length. Input as length, will be
%        convert to percent in the code.
%ten == tendon length
%Fitting == fitting length
%Outputs:
%F == Force, N
if nargin == 5              %Use if comparing measured muscle length to resting
    tendon = 0;
    fitting = 0;
elseif nargin == 6          %Use if measured muscle length and "tendon" where tendon can be a stand in for everything else in the Lmt that isn't active muscle
    tendon = ten;
    fitting = 0;
elseif nargin == 7          %Use if muscle length, tendon, and fitting sizes are known
    tendon = ten;
    fitting = Fitting;
else
    fprintf('Invalid number of arguments\n')
end
           
kmax = (rest-kmax)/rest; %Convert maximum contraction from length into percent
P = pres/620;             %Normalize the pressure

k = zeros(size(Lmt,1),1);
rel = zeros(size(Lmt,1),1);

    for i = 1:size(Lmt, 1)
        k(i,1) = (rest-(Lmt(i,1)-tendon-2*fitting))/rest; %current strain
        rel(i,1) = k(i,1)/kmax; %relative strain 
    end

    if ~isstring(dia)
        dia = num2str(dia);
    end
    maxF = maxBPAforce(rest,dia);

    if dia == 10
        load FestoLookup.mat f_10
        Fn = f_10(rel,P);
    elseif dia == 20
       load FestoLookup.mat f20
       Fn = f20(rel,P);
    elseif dia == 40
        load FestoLookup.mat f40
        Fn = f40(rel,P);
    end

    for i = length(Fn)
        if Fn(i) > 1
            Fn(i) = NaN;
        elseif Fn(i)<0
            Fn(i) = 0;
        else
        end
    end
    
    F = Fn.*maxF;
            
end