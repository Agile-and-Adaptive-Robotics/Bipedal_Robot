function T=fr_robotKneeT(k,degrees)
% Map OpenSim tibia coordinates through the registered robot t1 frame.
a=ppval(k.theta,degrees);
c=cos(a); s=sin(a);
T=[c,-s,0,0;s,c,0,0;0,0,1,0;0,0,0,1];
T(1:3,4)=ppval(k.translation,degrees);
T=T*k.adapter;
end
