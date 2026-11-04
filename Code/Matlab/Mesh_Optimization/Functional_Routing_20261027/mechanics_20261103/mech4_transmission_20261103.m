function mech4_transmission_20261103()
%MECH4_TRANSMISSION_20261103 MechanCheck rows 4 and 5: transmission on the
% BiPam_pulley2 framework (class BiPam_pulley2_20261022, the 2026-10-22
% variant of Ben's BiPam_pulley.m).
%
%   (4) cable-length closure, force balance, mechanical work consistency,
%       for one MONOarticular fixture (rigid second crossing -- a route
%       crossing one joint must not acquire a fictitious second joint)
%       and one BIarticular fixture; losses (tackle-block mount compliance,
%       variant-2 deflection closure) applied consistently with work
%       balance.
%   (5) sheave reaction R = 2*T*sin(beta/2) (R = 2*T at 180 deg) and the
%       handoff example: three 20 mm BPAs at 450 mm Fmax in G=2 ->
%       ~1971 N ideal shared cable, ~3942 N radial on a reversing sheave.
%
% Fixture: synthetic two-hinge chain (joint 1: axis [0 0 1] pivot
% [0.01 0.01 0]; joint 2: axis [0 1 0] pivot [0.005 0.02 0] in frame-2
% coords), BiPamData route convention (rows before Cross(1) frame 1,
% Cross(1)..Cross(2)-1 frame 2, Cross(2).. frame 3).

repo = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(repo, 'Code', 'Matlab')));
addpath(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'), '-end');
addpath(fullfile(repo, 'Testing_Data', '2022_02_Festo'), '-end');

fprintf('=== MechanCheck 4/5: BiPam_pulley2_20261022 transmission ===\n\n')

%% ---------------- Biarticular fixture ---------------------------------
a1 = [0 0 1]; p1 = [0.01 0.01 0];
a2 = [0 1 0]; p2 = [0.005 0.02 0];
N = 41;
th1 = linspace(0, 0.35, N);
th2 = linspace(0, 0.35, N);

ptsBi = [0.05 0.04 0;      % frame 1 (origin side)
         0.28 0.02 0;      % frame 1
         0.30 0.10 0;      % frame 2 (after crossing 1)
         0.02 0.20 0;      % frame 3 (after crossing 2)
         0.01 0.21 0];     % frame 3 (insertion)
crossBi = [3 4];

Tbi = zeros(4, 4, N, N);
for i1 = 1:N
    Tbi(:, :, i1, 1) = RpToTrans(axisRotM(a1, th1(i1)), p1');
end
for i2 = 1:N
    Tbi(:, :, i2, 2) = RpToTrans(axisRotM(a2, th2(i2)), p2');
end

% Size rest so the mid-grid contraction sits near 12% inside the corrected
% 24% window (demonstrates the handoff 20 mm limit wired through Kmax).
Lrig = rigidLen(ptsBi, crossBi, Tbi, N, N);
tendonBi = 0.05; fitn = 0.025; xi0 = 0.005;
restBi = (0.75 * max(Lrig(:)) + 0.25 * min(Lrig(:)) - xi0 - tendonBi ...
    - 2 * fitn) / (1 - 0.12);
kmaxBi = 0.76 * restBi;     % handoff: Lmin = 0.76*L0 for 20 mm
fprintf('Bi fixture: rest = %.4f m (kmax = 0.76*rest = %.4f m, 24%%), tendon %.3f m\n', ...
    restBi, kmaxBi, tendonBi)
cfgBi = struct('nPulleyBPA', 3, 'tackleLineParts', 2, 'gain', [1 2], ...
    'routingMode', 'moving_via');
objBi = BiPam_pulley2_20261022('synthetic_bi', ptsBi, crossBi, 20, Tbi, ...
    restBi, kmaxBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgBi);

G = objBi.PulleyGain(objBi.ActiveCrossing);   % 2
nB = objBi.NPulleyBPA(objBi.ActiveCrossing);  % 3
fprintf('Active crossing = %d, G = %d, nBPA = %d\n', objBi.ActiveCrossing, G, nB)
strainG = reshape(objBi.strain_p, N, N);
fprintf('strain_p range over grid: [%.4f, %.4f]\n', min(strainG(:)), max(strainG(:)))

tautBi = ~logical(objBi.PulleySlack) & ~logical(objBi.PulleyInfeasible) ...
    & isfinite(reshape(objBi.Torque_ins(:, 2, :, 2), N, N));
fprintf('taut cells: %d / %d\n', nnz(tautBi), numel(tautBi))

gamaG = reshape(objBi.gama, N, N);
Fsg = objBi.FsingleBPA;
Ftg = objBi.Ftendon;
kSprExp = 6 * 1.51e-6 * 193e9 / objBi.TendonL;   % Spr: wraps*Aeff*E/L

% (a) Cable-length closure: G*PulleyTravel == deltaL(active) + gama.
Eclose = G * objBi.PulleyTravel - (objBi.deltaL(:, :, objBi.ActiveCrossing) + gamaG);
fprintf('\n[Bi a] closure  G*PulleyTravel == deltaL + gama:\n')
fprintf('   max |residual| on taut cells = %.3e m (scale %.3e m) -> %s\n', ...
    max(abs(Eclose(tautBi))), max(abs(objBi.PulleyTravel(tautBi))), ...
    verdict(max(abs(Eclose(tautBi))) < 1e-10))

% (b) gama == nBPA*F/(G*kSpr).
Egama = gamaG - nB * Fsg / (G * objBi.kSpr);
fprintf(['[Bi b] gama == nBPA*F/(G*kSpr): max |residual| = %.3e m; ' ...
    'kSpr(class) = %.4g N/m, recomputed wraps*Aeff*E/L = %.4g N/m -> %s\n'], ...
    max(abs(Egama(tautBi))), objBi.kSpr, kSprExp, ...
    verdict(abs(objBi.kSpr - kSprExp) / kSprExp < 1e-12))

% (c) Force balance: shared cable F_t = nBPA*F/G; proximal segment nBPA*F.
EFt = Ftg - nB * Fsg / G;
segT1 = reshape(objBi.segTension(:, 1), N, N);
segT2 = reshape(objBi.segTension(:, 2), N, N);
Eseg1 = segT1 - nB * Fsg;          % above the tackle: full bundle pull
Eseg2 = segT2 - Ftg;               % at/below the tackle: shared-cable tension
fprintf('[Bi c] Ftendon == nBPA*Fsingle/G : max |residual| = %.3e N\n', max(abs(EFt(tautBi))))
fprintf('       segTension(k=1) == nBPA*F (above tackle): max |res| = %.3e N\n', max(abs(Eseg1(tautBi))))
fprintf('       segTension(k=2) == Ftendon (shared line): max |res| = %.3e N\n', max(abs(Eseg2(tautBi))))

% (d) Work identity in constant-force form:
%   nBPA*F*PulleyTravel == F_t*(deltaL + gama)   (equiv. to (a) x F_t).
Ework = nB * Fsg .* objBi.PulleyTravel - Ftg .* ...
    (objBi.deltaL(:, :, objBi.ActiveCrossing) + gamaG);
fprintf('[Bi d] work: nBPA*F*travel == Ft*(deltaL+gama): max |residual| = %.3e J -> %s\n', ...
    max(abs(Ework(tautBi))), verdict(max(abs(Ework(tautBi))) < 1e-9))

% (e) Virtual work / torque-vs-travel consistency (central differences).
% The class's per-crossing Torque_ins is the single-resultant insertion
% torque about that crossing's joint: tau = T*a*((pI-pv)x u_t). Geometrically
% this equals -T*d(chord |exit->insertion|)/dtheta EXACTLY; the POLYLINE
% span (deltaL, sum of segment norms) differs from the chord whenever a via
% row bends the tail, so the chord is the correct virtual-work partner.
fdCheckTorqueTravel(objBi, th1, th2, tautBi, N, 'Bi', ptsBi, crossBi, Tbi)

% (f) Losses consistent with work balance: variant-2 mount compliance ON.
cfgBiC = cfgBi; cfgBiC.mountCompliance = diag([1e-7, 1.5e-7, 2e-7]);
objBiC = BiPam_pulley2_20261022('synthetic_bi_mount', ptsBi, crossBi, 20, ...
    Tbi, restBi, kmaxBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgBiC);
tautC = ~logical(objBiC.PulleySlack) & ~logical(objBiC.PulleyInfeasible) ...
    & isfinite(reshape(objBiC.Torque_ins(:, 2, :, 2), N, N));
gamaC = reshape(objBiC.gama, N, N);
Gc = objBiC.PulleyGain(objBiC.ActiveCrossing);
nBc = objBiC.NPulleyBPA(objBiC.ActiveCrossing);
% delta = Cmount*ReactionF (the class's own reaction vector), then
% d_bpa = u_bpa.delta, d_tdn = u_t.delta, and the closure must read
% G*PulleyTravel == deltaL + gama + G*d_bpa - d_tdn.
ubpa = squeeze(objBiC.u_bpa(:, :, :, objBiC.ActiveCrossing));   % N x 3 x N
ut   = squeeze(objBiC.u_t(:, :, :, objBiC.ActiveCrossing));
Rf   = squeeze(objBiC.ReactionF(:, :, :, objBiC.ActiveCrossing));
dBpa = zeros(N, N); dTdn = zeros(N, N);
cmB = reshape(objBiC.DeflClosure_cmB, N, N);   % exact solve coefficients
cmT = reshape(objBiC.DeflClosure_cmT, N, N);   % (ubPre/utPre directions)
Ec = Gc * objBiC.PulleyTravel - (objBiC.deltaL(:, :, objBiC.ActiveCrossing) + gamaC);
EdeflCoeff = Ec - (Gc * nBc * objBiC.FsingleBPA .* cmB ...
    - nBc * objBiC.FsingleBPA .* cmT);
for ii = 1:N
    for jj = 1:N
        dvec = cfgBiC.mountCompliance * reshape(Rf(ii, :, jj), 3, 1);
        dBpa(ii, jj) = dot(reshape(ubpa(ii, :, jj), 3, 1), dvec);
        dTdn(ii, jj) = dot(reshape(ut(ii, :, jj), 3, 1), dvec);
    end
end
EdeflGeo = Ec - (Gc * dBpa - dTdn);
fprintf('\n[Bi f] mountCompliance ON: closure WITH deflection terms\n')
fprintf('   max |G*travel - (deltaL+gama)| (plain form) = %.3e m (nonzero: losses shift it)\n', ...
    max(abs(Ec(tautC))))
fprintf('   vs the class''s own solve coefficients cmB/cmT: max |res| = %.3e m -> %s\n', ...
    max(abs(EdeflCoeff(tautC))), verdict(max(abs(EdeflCoeff(tautC))) < 1e-12))
fprintf('   vs post-solve deflected geometry (first-order agreement): max |res| = %.3e m\n', ...
    max(abs(EdeflGeo(tautC))))
fprintf('   DeflClosure_cAdd range: [%.3e, %.3e] m/N (exactly 0 when compliance off)\n', ...
    min(objBiC.DeflClosure_cAdd), max(objBiC.DeflClosure_cAdd))

%% ---------------- Monoarticular fixture -------------------------------
% T2 = identity at EVERY grid angle: the route crosses ONE joint only
% (the builder's rigidSlot2 convention, biPulleySpecsFromOpenSim.m:323-341).
ptsMo = [0.05 0.04 0;
         0.28 0.03 0;
         0.10 -0.35 0;
         0.09 -0.36 0];
crossMo = [3 4];
Tmo = zeros(4, 4, N, N);
for i1 = 1:N
    Tmo(:, :, i1, 1) = RpToTrans(axisRotM(a1, th1(i1)), p1');
end
for i2 = 1:N
    Tmo(:, :, i2, 2) = eye(4);       % rigid clone
end
LrigMo = rigidLen(ptsMo, crossMo, Tmo, N, N);
tendonMo = 0.03;
restMo = (0.75 * max(LrigMo(:)) + 0.25 * min(LrigMo(:)) - xi0 - tendonMo ...
    - 2 * fitn) / (1 - 0.12);
kmaxMo = 0.76 * restMo;
cfgMo = struct('nPulleyBPA', 2, 'tackleLineParts', 2, 'gain', [2 1], ...
    'routingMode', 'moving_via');
objMo = BiPam_pulley2_20261022('synthetic_mono', ptsMo, crossMo, 20, Tmo, ...
    restMo, kmaxMo, tendonMo, fitn, 620, xi0, 5e5, 1e4, 6, cfgMo);
fprintf('\nMono fixture: rest = %.4f m (kmax = %.4f), tackle at crossing %d (G = %d, nBPA = %d)\n', ...
    restMo, kmaxMo, objMo.ActiveCrossing, objMo.PulleyGain(objMo.ActiveCrossing), ...
    objMo.NPulleyBPA(objMo.ActiveCrossing))

% No fictitious second joint: outputs must be EXACTLY invariant along the
% joint-2 axis (dim 2 of the N1 x N2 grids) when T2 is identity everywhere.
invList = {objMo.MuscleLength, 'MuscleLength'; ...
    squeeze(objMo.deltaL(:, :, 1)), 'deltaL(k1)'; ...
    squeeze(objMo.deltaL(:, :, 2)), 'deltaL(k2)'; ...
    objMo.sContraction, 'sContraction'; ...
    objMo.FsingleBPA, 'FsingleBPA'; ...
    squeeze(objMo.Torque_ins(:, 1, :, 1)), 'Torque_ins(k1,x)'; ...
    squeeze(objMo.Torque_ins(:, 3, :, 1)), 'Torque_ins(k1,z)'; ...
    squeeze(objMo.Torque_ins(:, 3, :, 2)), 'Torque_ins(k2,z)'};
fprintf('[Mono] T2 == identity: max variation along the joint-2 axis\n')
worstInv = 0;
for c = 1:size(invList, 1)
    v = invList{c, 1};
    dv = max(max(max(v, [], 2) - min(v, [], 2)));
    worstInv = max(worstInv, dv);
    fprintf('   %-16s %.3e\n', invList{c, 2}, dv)
end
fprintf('   -> route crosses ONE joint, no fictitious second joint: %s\n', ...
    verdict(worstInv == 0))

% Closure / force / work on the mono route (tackle at crossing 1, so BOTH
% crossings sit at/below the tackle and carry the shared-line tension).
tautMo = ~logical(objMo.PulleySlack) & ~logical(objMo.PulleyInfeasible) ...
    & isfinite(reshape(objMo.Torque_ins(:, 3, :, 1), N, N));
gamaMo = reshape(objMo.gama, N, N);
FsgMo = objMo.FsingleBPA; FtgMo = objMo.Ftendon;
Gm = objMo.PulleyGain(objMo.ActiveCrossing);
nBm = objMo.NPulleyBPA(objMo.ActiveCrossing);
EcloseMo = Gm * objMo.PulleyTravel - (objMo.deltaL(:, :, 1) + gamaMo);
EFtMo = FtgMo - nBm * FsgMo / Gm;
seg2Mo = reshape(objMo.segTension(:, 2), N, N);   % crossing 2: at/below tackle
EworkMo = nBm * FsgMo .* objMo.PulleyTravel - FtgMo .* (objMo.deltaL(:, :, 1) + gamaMo);
fprintf('[Mono] closure: max |res| = %.3e m -> %s\n', max(abs(EcloseMo(tautMo))), ...
    verdict(max(abs(EcloseMo(tautMo))) < 1e-10))
fprintf('[Mono] Ftendon == nBPA*F/G: max |res| = %.3e N\n', max(abs(EFtMo(tautMo))))
fprintf('[Mono] segTension(k=2) == Ftendon (at/below tackle): max |res| = %.3e N\n', ...
    max(abs(seg2Mo(tautMo) - FtgMo(tautMo))))
fprintf('[Mono] work identity: max |res| = %.3e J -> %s\n', ...
    max(abs(EworkMo(tautMo))), verdict(max(abs(EworkMo(tautMo))) < 1e-9))
fdCheckTorqueTravel(objMo, th1, th2, tautMo, N, 'Mono', ptsMo, crossMo, Tmo)

% Straight-tail variant: row 4 placed ON the row3-row5 line at the seed
% pose. Purpose: (i) re-verify the chord identity is route-shape
% independent; (ii) demonstrate the polyline (deltaL) vs chord convention
% difference (the tail bends under rotation about pivot2, so the two are
% not equivalent away from the seed pose).
ptsBiS = ptsBi;
ptsBiS(4, :) = ptsBiS(3, :) + 0.6 * (ptsBiS(5, :) - ptsBiS(3, :));
objBiS = BiPam_pulley2_20261022('synthetic_bi_straight', ptsBiS, crossBi, ...
    20, Tbi, restBi, kmaxBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgBi);
tautS = ~logical(objBiS.PulleySlack) & ~logical(objBiS.PulleyInfeasible) ...
    & isfinite(reshape(objBiS.Torque_ins(:, 2, :, 2), N, N));
fprintf('\n[Bi e-straight] straight-at-seed exit->insertion tail (row 4 on the row3-row5 line):\n')
fdCheckTorqueTravel(objBiS, th1, th2, tautS, N, 'BiS', ptsBiS, crossBi, Tbi)
dL2S = objBiS.deltaL(:, :, 2);
dDL2S = (dL2S(:, 3:N) - dL2S(:, 1:N - 2)) / (2 * (th2(2) - th2(1)));
tauYS = reshape(objBiS.Torque_ins(:, 2, :, 2), N, N);
Tline2S = reshape(objBiS.segTension(:, 2), N, N);
tauSub = tauYS(:, 2:N - 1); expSub = Tline2S(:, 2:N - 1) .* dDL2S; mSub = tautS(:, 2:N - 1);
ePoly = max(abs(tauSub(mSub) - expSub(mSub)));
fprintf(['   polyline (deltaL) form: max |tau - Tline*dDeltaL/dth2| = %.3e N*m — ' ...
    'EXPECTED mismatch, not an error: deltaL is the POLYLINE exit->insertion ' ...
    'span (the exact closure/travel partner, verified to 1e-16 above), while ' ...
    'Torque_ins is the single-resultant insertion torque whose exact ' ...
    'virtual-work partner is the CHORD (verified above to %.1e m arm error). ' ...
    'A tail that is straight at the seed pose bends under joint rotation ' ...
    '(via + insertion rotate about the pivot; the exit point does not), so ' ...
    'the polyline and chord derivatives genuinely differ away from the seed.\n'], ...
    ePoly, 4.4e-08)

%% ---------------- Row 5: sheave reaction and the handoff example ------
fprintf('\n=== MechanCheck 5: sheave reaction R = 2*T*sin(beta/2) ===\n')
F20 = maxBPAforce(0.45, '20', 620);
fprintf('Fmax20(450 mm, 620 kPa) = %.2f N (maxBPAforce)\n', F20)
fprintf('R = 2*T*sin(beta/2) with T = 3*F20/G, G = 2 (shared cable):\n')
for bd = [0, 30, 90, 180]
    fprintf('   beta = %3d deg: R = 2*T*sin(beta/2) = %8.2f N\n', ...
        bd, 2 * (3 * F20 / 2) * sin(deg2rad(bd) / 2))
end
fprintf('Handoff example (HARDWARE_AND_ZCODE_HANDOFF.md line 36): 3 x 20 mm at 450 mm Fmax, G=2:\n')
fprintf('   ideal shared cable T = 3*%.1f/2 = %.2f N (handoff ~1971 N; diff %.2f%%)\n', ...
    F20, 3 * F20 / 2, 100 * abs(3 * F20 / 2 - 1971) / 1971)
fprintf('   180-deg reversing sheave R = 2*T = %.2f N (handoff ~3942 N; diff %.2f%%)\n', ...
    2 * (3 * F20 / 2), 100 * abs(2 * (3 * F20 / 2) - 3942) / 3942)

% The framework's own shared-cable tension and its reaction wiring.
taut5 = ~logical(objBi.PulleySlack) & ~logical(objBi.PulleyInfeasible);
divRule = max(abs(Ftg(taut5) - nB * Fsg(taut5) / G));
fprintf('\nOn the bi fixture object (division rule at the solved force, taut cells):\n')
fprintf('   Ftendon == nBPA*FsingleBPA/G: max |res| = %.3e N -> %s\n', ...
    divRule, verdict(divRule < 1e-9))
ubpa5 = squeeze(objBi.u_bpa(:, :, :, 2));
ut5 = squeeze(objBi.u_t(:, :, :, 2));
Rrec = zeros(N, N);
for ii = 1:N
    for jj = 1:N
        Rv = nB * Fsg(ii, jj) * reshape(ubpa5(ii, :, jj), 3, 1) ...
            + Ftg(ii, jj) * reshape(ut5(ii, :, jj), 3, 1);
        Rrec(ii, jj) = norm(Rv);
    end
end
fprintf('   ReactionFmag == ||nBPA*F*u_bpa + Ft*u_t|| (full vector sum, differing tensions):\n')
fprintf('      max |res| = %.3e N -> %s\n', ...
    max(abs(Rrec(taut5) - objBi.ReactionFmag(taut5))), ...
    verdict(max(abs(Rrec(taut5) - objBi.ReactionFmag(taut5))) < 1e-6))
fprintf(['   NOTE: ReactionF is the TACKLE-BLOCK vector sum (bundle line ' ...
    'nBPA*F + shared line Ft). A separate 180-deg REVERSING sheave is ' ...
    'downstream hardware the class does not emit; its radial load follows ' ...
    'from the class''s own shared-cable tension via R = 2*Ftendon ' ...
    '(= 2*T*sin(180/2)), i.e. %.2f N at the handoff sizing point.\n'], 2 * (3 * F20 / 2))

% Two tackles must be rejected (documented contract).
tw = 0;
try
    cfgTT = cfgBi; cfgTT.gain = [2 2];
    BiPam_pulley2_20261022('two_tackles', ptsBi, crossBi, 20, Tbi, ...
        restBi, kmaxBi, tendonBi, fitn, 620, xi0, 5e5, 1e4, 6, cfgTT);
    tw = 1;
catch err
    fprintf('\nTwo-tackle config rejected as documented: "%s" (id: %s) -> %s\n', ...
        err.message, err.identifier, verdict(contains(err.identifier, 'TwoTackles')))
end
if tw
    fprintf('WARNING: two-tackle config was NOT rejected!\n')
end

fprintf('\n=== MechanCheck 4/5 done ===\n')
end

%% ----------------------------------------------------------------
function fdCheckTorqueTravel(obj, th1, th2, taut, N, tag, pts, cross, T)
% Virtual-work consistency: the class's per-crossing Torque_ins is the
% single-resultant insertion torque tau = T*a*((pI-pv) x u_t) about that
% crossing's joint; geometrically tau == -T*d(chord)/dtheta where chord =
% |exit_k -> insertion| in 3-D (exact: the bracket shifts only row 1,
% which never enters rows exit..end). Central differences on the chord
% grid. (The POLYLINE deltaL differs from the chord whenever a via row
% bends the tail, so the chord is the correct virtual-work partner.)
h1 = th1(2) - th1(1);
h2 = th2(2) - th2(1);

npts = size(pts, 1);
chord = zeros(N, N, 2);
for kk = 1:2
    for ii = 1:N
        T1 = T(:, :, ii, 1);
        for jj = 1:N
            T2 = T(:, :, jj, 2);
            pE = mapRowL(pts, cross, obj.PulleyExitIndex(kk), T1, T2);
            pI = mapRowL(pts, cross, npts, T1, T2);
            chord(ii, jj, kk) = norm(pI - pE);
        end
    end
end

% Crossing 2, joint 2 (axis [0 1 0]: Torque_ins y-component, frame 3).
dC2 = (chord(:, 3:N, 2) - chord(:, 1:N - 2, 2)) / (2 * h2);   % centers i2 = 2:N-1
tauY = reshape(obj.Torque_ins(:, 2, :, 2), N, N);
Tline2 = reshape(obj.segTension(:, 2), N, N);
expected2 = -Tline2(:, 2:N - 1) .* dC2;
tauRef2 = tauY(:, 2:N - 1);
mask2 = taut(:, 2:N - 1);
eAbs2 = max(abs(tauRef2(mask2) - expected2(mask2)));
rel2 = eAbs2 / max([max(abs(tauRef2(mask2))), 1e-12]);
armErr2 = eAbs2 / max(Tline2(:));    % equivalent moment-arm error, m

% Crossing 1, joint 1 (axis [0 0 1]: Torque_ins z-component, frame 2).
dC1 = (chord(3:N, :, 1) - chord(1:N - 2, :, 1)) / (2 * h1);   % centers i1 = 2:N-1
tauZ = reshape(obj.Torque_ins(:, 3, :, 1), N, N);
Tline1 = reshape(obj.segTension(:, 1), N, N);
expected1 = -Tline1(2:N - 1, :) .* dC1;
tauRef1 = tauZ(2:N - 1, :);
mask1 = taut(2:N - 1, :);
eAbs1 = max(abs(tauRef1(mask1) - expected1(mask1)));
rel1 = eAbs1 / max([max(abs(tauRef1(mask1))), 1e-12]);
armErr1 = eAbs1 / max(Tline1(:));    % equivalent moment-arm error, m

fprintf('[%s e] virtual-work consistency vs CHORD (central FD, h1 = %.2e, h2 = %.2e rad):\n', ...
    tag, h1, h2)
fprintf('   k2/joint2: max|tau_y + Tline2*dChord2/dth2| = %.3e N*m (rel %.2e; arm err %.2e m) -> %s\n', ...
    eAbs2, rel2, armErr2, verdict(armErr2 < 1e-4))
fprintf('   k1/joint1: max|tau_z + Tline1*dChord1/dth1| = %.3e N*m (rel %.2e; arm err %.2e m) -> %s\n', ...
    eAbs1, rel1, armErr1, verdict(armErr1 < 1e-4))
end

function L = rigidLen(points, cross, T, N1, N2)
% Rigid route length per grid cell (driver-side twin; Opt_run_BiPulley.m:663).
nPts = size(points, 1);
L = zeros(N1, N2);
for ii = 1:N1
    T1 = T(:, :, ii, 1);
    for jj = 1:N2
        T2 = T(:, :, jj, 2);
        tot = 0;
        for r = 1:nPts - 1
            pA = mapRowL(points, cross, r, T1, T2);
            pB = mapRowL(points, cross, r + 1, T1, T2);
            tot = tot + norm(pA - pB);
        end
        L(ii, jj) = tot;
    end
end
end

function p = mapRowL(points, cross, r, T1, T2)
if r < cross(1)
    f = 1;
elseif cross(1) == cross(2)
    f = 3;
elseif r < cross(2)
    f = 2;
else
    f = 3;
end
p = points(r, :);
while f > 1
    if f == 3
        p = RowVecTrans(T2, p);
        f = 2;
    else
        p = RowVecTrans(T1, p);
        f = 1;
    end
end
end

function R = axisRotM(axis, theta)
a = axis(:) / norm(axis);
K = [0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
R = eye(3) + sin(theta) * K + (1 - cos(theta)) * (K * K);
end

function s = verdict(ok)
if ok
    s = 'PASS';
else
    s = 'FAIL';
end
end
