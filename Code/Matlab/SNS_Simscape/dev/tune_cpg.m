% tune_cpg.m — fast ODE prototype of the half-center CPG to pick gains
% before touching the Simulink build.
% State: Vx (RG ext), Vf (RG flex), A_x, A_f (Adp neurons)
% RG:  Cm*dV/dt = Gm*(Vr-V) + Idrive + gmi*sat(Vp)*(E-V) + gadp*sat(A)*(E-V)
% Adp: Cma*dA/dt = Gma*(Vr-A) + ge*sat(Vrg)*(0-A)
% sat(u) = clip((u-ThrPre)*SlopePre, 0, 1); E = -72 (inh); Vr = -52

Vr = -52; ThrP = -45; SlP = 0.5; E = -72;
sat = @(u) min(max((u - ThrP)*SlP, 0), 1);

p.gmi  = 0.8;   p.gse = 0.6;
p.gadp = 2.0;
p.ge   = 0.12;
p.Gm   = 0.2;  p.Cm = 5;
p.Gma  = 0.05; p.Cma = 30;
p.Ix = 15; p.If = 6;    % strong asymmetry: ext wins the onset race

dt = 0.1;
T = 30000; N = T/dt;
Vx = Vr; Vf = Vr - 0.5; Ax = Vr; Af = Vr;
outVx = zeros(N,1); outVf = zeros(N,1);
for k = 1:N
    Ix  = p.Ix + p.gmi*sat(Vf)*(E-Vx) + p.gadp*sat(Ax)*(E-Vx) + p.gse*sat(Vx)*(0-Vx);
    If  = p.If + p.gmi*sat(Vx)*(E-Vf) + p.gadp*sat(Af)*(E-Vf) + p.gse*sat(Vf)*(0-Vf);
    dVx = (p.Gm*(Vr-Vx) + Ix)/p.Cm;
    dVf = (p.Gm*(Vr-Vf) + If)/p.Cm;
    dAx = (p.Gma*(Vr-Ax) + p.ge*sat(Vx)*(0-Ax))/p.Cma;
    dAf = (p.Gma*(Vr-Af) + p.ge*sat(Vf)*(0-Af))/p.Cma;
    Vx = Vx + dt*dVx; Vf = Vf + dt*dVf; Ax = Ax + dt*dAx; Af = Af + dt*dAf;
    outVx(k) = Vx; outVf(k) = Vf;
end
t = (1:N)'*dt/1000;
d = sign(outVx - outVf); d(d==0) = 1;
sw = sum(abs(diff(d)) > 0);
fprintf('gmi=%.2f gadp=%.2f ge=%.2f Ix=%.1f If=%.1f -> %d switches', p.gmi, p.gadp, p.ge, p.Ix, p.If, sw);
if sw >= 2
    idx = find(abs(diff(d)) > 0);
    fprintf(', period ~ %.2f s\n', 2*mean(diff(t(idx))));
else
    fprintf(' (LATCHED)\n');
end
% print a few samples
for k = round(linspace(1, N, 8))
    fprintf('  t=%5.2f s  Vx=%7.2f  Vf=%7.2f\n', t(k), outVx(k), outVf(k));
end
