Sp = 85000;
At = 0.606;

Fi = 0.75*Sp*At*10^(-3)

C = 1/6;

Fbmin = Fi;

Fbmax = Fi + C*10000*10^(-3)

sigmin = Fbmin/At

sigmax = Fbmax/At

sigm = (sigmin + sigmax)/2

siga = (sigmax - sigmin)/2

sigi = Fi/At

Se = 18.6;
Sut = 120;

n = Se*(Sut - sigi)/(Sut*siga + Se*(sigm - sigi))

r = sqrt(2)

Fdd = 12000/(4*r)

Fby = 250 + sin(deg2rad(45))*Fdd

Fbx = sin(deg2rad(45))*Fdd

Fb = sqrt(Fby^2 + Fbx^2)

tau = Fb/(pi/4*0.5^2)

n = 0.577*57*10^3/tau