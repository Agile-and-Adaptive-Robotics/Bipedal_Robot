function C = MuscleLengthCost(C, Muscle, GLength)
    for ii = 1:length(Muscle.MuscleLength)
        C(k) = C(k) + GLength*Muscle.MuscleLength(ii);
    end
end