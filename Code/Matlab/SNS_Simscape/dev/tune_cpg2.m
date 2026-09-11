% tune_cpg2.m — half-center v2: adaptation current g*S_adp subtracted from RG
Vr = -52; ThrP = -45; SlP = 0.5; E = -72;
satP = @(u) min(max((u - ThrP)*SlP, 0), 1);   % synapse presat (SlopePre 0.5)
satN = @(u) min(max(u + 45, 0), 1);           % RG output S (Thr -45, Slope 1)

p.gmi   = 0.8;    % mutual inhibition
p.gadp  = 15;     % adaptive current (nA at S_adp = 1)
p.ge    = 0.12;   % RG -> Adp excitation (Esyn 0)
p.Gm    = 0.2;  p.Cm = 5;      % RG tau 25 ms
p.Gma   = 0.05; p.Cma = 30;    % Adp tau 600 ms
p.tauz  = 500;    % ms (Cma/Gma)
p.Ix = 4.0; p.If = 3.6;        % drives: free V = -32 / -34 mV

dt = 0.2; T = 30000; N = T/dt;
Vx = Vr; Vf = Vr - 0.5; zx = 0; zf = 0;
% Adp neuron voltage state: dzx/dt = (satN(Vx) - zx)/tauz via RC: use z directly
out = zeros(N, 4);
for k = 1:N
    Ix  = p.Ix + p.gmi*satP(Vf)*(E-Vx) - p.gadp*zx;
    If  = p.If + p.gmi*satP(Vx)*(E-Vf) - p.gadp*zf;
    dVx = (p.Gm*(Vr-Vx) + Ix)/p.Cm;
    dVf = (p.Gm*(Vr-Vf) + If)/p.Cm;
    dzx = (satN(Vx) - zx)/p.tauz;
    dzf = (satN(Vf) - zf)/p.tauz;
    Vx = Vx + dt*dVx; Vf = Vf + dt*dVf;
    zx = min(max(zx + dt*dzx, 0), 1); zf = min(max(zf + dt*dzf, 0), 1);
    out(k,:) = [Vx, Vf, zx, zf];
end
t = (1:N)'*dt/1000;
d = sign(out(:,1) - out(:,2)); d(d==0) = 1;
sw = sum(abs(diff(d)) > 0);
fprintf('gmi=%.1f gadp=%.0f Ix=%.1f If=%.1f -> %d switches', p.gmi, p.gadp, p.Ix, p.If, sw);
if sw >= 2
    idx = find(abs(diff(d)) > 0);
    fprintf(', period ~ %.2f s\n', 2*mean(diff(t(idx))));
else
    fprintf(' (LATCHED)\n');
end
Sx = satN(out(:,1)); Sf = satN(out(:,2));
fprintf('S_ext range %.2f..%.2f, S_flex range %.2f..%.2f\n', min(Sx), max(Sx), min(Sf), max(Sf));
for k = round(linspace(1, N, 10))
    fprintf('  t=%5.2f s  Vx=%7.2f (S=%.2f)  Vf=%7.2f (S=%.2f)\n', t(k), out(k,1), satN(out(k,1)), out(k,2), satN(out(k,2)));
end
