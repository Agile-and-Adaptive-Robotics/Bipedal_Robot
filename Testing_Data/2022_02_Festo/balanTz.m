function [x, fval] = balanTz(Angle, Torque, InflatedLength, ICRtoMuscle, pres, class)
    %Lm == length of muscle (BPA), measured
    %Lmt == length of musculotendon, predicted
    %pres == BPA pressure
    %rest == resting length
    %dia == BPA diameter
    %KMAX == maximum contracted length
    
    Dia = num2str(class.Diameter);
    rest = obj.RestingL;
    kmax = class.Kmax;  
    KMAX = (rest-kmax)/rest; %turn it into a percentage
    fitting = class.FittingLength;
    tendon = class.TendonL;
    Lmt = class.MuscleLength;
    contract = class.Contraction;
    unitD = class.UnitDirection;
    MA = class.MomentArm;
    MAz = sqrt(MA(:,1).^2+MA(:,2).^2);
    Fvec = class.Force;
    Fpred = vecnorm(Fvec,2,2);
    Fmax = class.Fmax;
    Tcalmag = vecnorm(class.Torque,2,2);
    Tcalz = class.Torque(:,3);
            
    relpred = contract./KMAX;        %relative strain        

    y = pres;    %pressure
    relPres = pres/620;              %relative pressure
    
    Fmax = maxBPAforce(rest,Dia);
    if dia == 10
        load FestoLookup.mat f_10
        Fn = f_10;
    elseif dia == 20
        load FestoLookup.mat f20
        Fn = f20;
    elseif dia == 40
        load FestoLookup.mat f40
       Fn = f40;
    else
        disp('Error with uninflated diameter')
    end
    
    d = optimvar('d',1,'LowerBound',0);
    kSE = optimvar('y',1,'LowerBound',0);
    
    function Fbal = myfunc(u)

        strain = (rest-(Lmt-u(1)-tendon-2*fitting))/rest;      
        rel = strain./kmax;
        Fbpa = Fmax.*Fn(rel,relPres); %

        Ft = u(2).*u(1); %Series elastic element

        Fbal(1) = Fbpa-Ft;
        Fbal(2) = Lmt - (Lm+tendon+2.*fitting+u(1));

    end 

    x0(1) = 0.01; %Initial guess
    x0(2) = 1000;         %Spring rate, N/m    
    options = optimoptions('fsolve','Display','iter','FunctionTolerance',0.001,'StepTolerance',1*10^-6,'PlotFcn',[@optimplotx, @optimplotfval]);
    
    %x0 is the guess of the Muscle length
    [x, fval] = fsolve(@myfunc,x0,options);
    
    strain2 = (rest-(Lm-x(1)))/rest;   
    rel2 = strain2./kmax;
    f = Fmax*(a0*(exp(-a1.*rel2)-1)+y.*exp(-a3*((rel2).^2)))-x(2).*x(1); %mif at Lmt
    if f < 0
           f = 0;  %Muscle force is tensile only
    end
end

% Hoy, M. G., Zajac, F. E., & Gordon, M. E. (1990). A musculoskeletal model of the human lower extremity: the effect of muscle, tendon, and moment arm on the moment-angle relationship of musculotendon actuators at the hip, knee, and ankle. Journal of biomechanics, 23(2), 157-169.
% 
% Thelen, D. G. (2003). Adjustment of muscle mechanics model parameters to simulate dynamic contractions in older adults. Journal of biomechanical engineering, 125(1), 70-77.
% 
% Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). Flexing computational muscle: modeling and simulation of musculotendon dynamics. Journal of biomechanical engineering, 135(2), 021005.