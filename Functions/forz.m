function f = forz(Lmt, mif, ofl, tsl, pa)
    %Lmt == muscle-tendon length
    %mif == max isometric force at optimum fiber length
    %ofl == optimal fiber length
    %tsl == tendon slack length
    %pa == Pennation angle

    kPE = 4;    %passive exponential shape factor (OpenSim Kshape Pasive)
    eom = 0.6;  %passive muscle strain at max isometric force (OpenSim FmaxMuscleStrain)
    y = 0.5;    %active shape factor (OpenSim KshapeActive)

        function Fbal = myfunc(Lmn)

            Fpe = (exp(kPE*((Lmn)-1)/eom)-1)/(exp(kPE)-1); %Passive force-length curve, normalized (Thelen 2003)
            fL = exp(-(((Lmn)-1).^2/y));  %Active force-length curve, normalized (Thelen 2003)
            cosa = ((1-(sin(pa)/(Lmn)).^2)^(1/2)); %cosine of pennation angle (Hoy 1990)
            fT = (37.5/(tsl/ofl))*((Lmt/ofl)-(Lmn)*cosa-(tsl/ofl)); %Normalized tendon force (Hoy 1990)
            Fbal = (fL+Fpe)*cosa-fT;

        end 

    if Lmt ~= tsl
        x0 = (Lmt-tsl)/ofl; %Initial guess
    else
        x0 = (Lmt*1.001-tsl)/ofl; %Initial guess, avoid computational crash
    end
    options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);
    Lma = fsolve(@myfunc,x0,options);
    f = mif*(37.5/(tsl/ofl))*((Lmt/ofl)-(Lma)*((1-(sin(pa)/(Lma))^2)^(1/2))-(tsl/ofl)); %mif at Lmt
    if f < 0
           f = 0;  %Muscle force is tensile only (Millard 2013)
    end
end

% Hoy, M. G., Zajac, F. E., & Gordon, M. E. (1990). A musculoskeletal model of the human lower extremity: the effect of muscle, tendon, and moment arm on the moment-angle relationship of musculotendon actuators at the hip, knee, and ankle. Journal of biomechanics, 23(2), 157-169.
% 
% Thelen, D. G. (2003). Adjustment of muscle mechanics model parameters to simulate dynamic contractions in older adults. Journal of biomechanical engineering, 125(1), 70-77.
% 
% Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). Flexing computational muscle: modeling and simulation of musculotendon dynamics. Journal of biomechanical engineering, 135(2), 021005.