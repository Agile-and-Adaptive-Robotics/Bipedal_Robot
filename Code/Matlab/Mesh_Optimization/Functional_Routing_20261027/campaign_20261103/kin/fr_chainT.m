function T = fr_chainT(kin, q, fromBody, toBody)
%FR_CHAINT Transform mapping points expressed in fromBody into toBody:
% p_to = T(1:3,1:3)*p_from + T(1:3,4).
Tf = fr_bodyT(kin, q, fromBody);
Tt = fr_bodyT(kin, q, toBody);
T = Tt \ Tf;
end
