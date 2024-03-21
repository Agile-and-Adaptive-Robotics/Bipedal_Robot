%% This function finds force of a 10mm BPA at a given pressure, relative strain, and resting length
%-Use this to back calculate scalar force values
%-Might not be appropriate to use an upper bound since we are back
% calculating force
function F = bpaForce10(restingL, relStrain, Pressure)
rest = restingL;
k = relStrain;
P = Pressure/620;

maxF = maxBPAforce(rest);

load FestoLookup.mat f_10
Fn = f_10(k,P);
F = Fn.*maxF;

%No negative force
F(F<0)=0;

%Maximum force limit
%limit of Fmax10 = 476.72 N
%maximum error = 30.9 N
% F(F>507.63)=NaN;          

end