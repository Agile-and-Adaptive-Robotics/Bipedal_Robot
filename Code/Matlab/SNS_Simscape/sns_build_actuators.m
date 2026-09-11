%% sns_build_actuators.m — add real BPA actuator blocks (10/20/40 mm) + a biological muscle block to SNS_Library.slx
%
% Built 2026-09-10. Uses BEN'S EQUATIONS verbatim:
%
% Normalized force surface (FestoLookup.mat sfits, exported to
% festo_lookup_coeffs.json; same surface festo4.m evaluates):
%   Fn(rel, Pn) = a0*(exp(-a1*rel) - 1) + Pn*exp(-a3*rel^2)
%   rel  = relative strain  = contraction / KMAX
%   Pn   = pressure / 620 kPa
%   dia     a0              a1              a3
%   10   0.568207874671  4.25442545542   0.55972777762
%   20   0.257852586017  6.4766142989    1.32087718059
%   40   0.122366088343  10.4714342293   2.02328814095
%
% Max isometric force at zero strain, 620 kPa (maxBPAforce.m):
%   10 mm: Fmax = 620*(0.4895*atan(0.03068*(Rest-0.0075)*620))   [N, RMSE 14.7 N]
%   20 mm: Fmax = 620*(1.4877*atan(0.0248*(Rest-0.0075)*620))    [N, RMSE 23.8 N]
%   40 mm: Fmax = 6000 N  (BalanceX3 class value; maxBPAforce says 6398.4
%          from the Festo tool — editable mask param, ask Ben which he wants)
%
% Contraction (MonoPamDataExplicit_balanceX3.get.Contraction):
%   contraction = (Rest - (L - TendonL - 2*FittingL)) / Rest
%   KMAX        = (Rest - Kmax) / Rest
%   F           = Fn * Fmax, then  F(rel >= 1) = 0,  F < 0 -> 0   (festo4.m)
%
% Biological muscle (BioMuscle): OpenSim Thelen2003-style Hill-type, rigid
% tendon:
%   lambda = l_ce/l_opt, l_ce = max(Lmt - l_slack, eps)/cos(alpha)
%   fl  = exp(-(lambda-1)^2 / 0.45)                       (active force-length)
%   fv  = Hill: shortening  fv = b*(1+A)/(s + b) - A,  s = -v_tilde >= 0,
%               b = A*v_max;  lengthening fv = 1.4 - 0.4*exp(-v_tilde/0.2)
%   fpas = (exp(kpe*(lambda-1)/e_pas) - 1)/(exp(kpe) - 1), lambda > 1, else 0
%   F   = Fmax * (Act*fl*fv + fpas) * cos(alpha), clamped [0, 2*Fmax]
%   v_tilde = v_mt / cos(alpha) / l_opt   (v_mt = musculotendon velocity m/s)
% Defaults: A = 0.25 (Hill curvature), v_max = 10 l_opt/s, e_pas = 0.6,
% kpe = 4, alpha = 0 deg. l_slack modeled rigid (tendon slack length only
% offsets fiber length); swap in series elasticity when needed.

lib = 'SNS_Library';
if ~bdIsLoaded(lib)
    load_system(fullfile(fileparts(mfilename('fullpath')), [lib '.slx']));
end
set_param(lib, 'Lock', 'off');   % libraries load read-only by default

%% ---------------- BPA 10 / 20 / 40 mm ----------------
% a0, a1, a3 per diameter (closed-form Festo surfaces, see header)
surf = { ...
    'BPA_10mm', 10,  0.568207874671, 4.25442545542, 0.55972777762, ...
    '(620*(0.4895*atan(0.03068*(Rest-0.0075)*620)))'; ...
    'BPA_20mm', 20,  0.257852586017, 6.4766142989,  1.32087718059, ...
    '(620*(1.4877*atan(0.0248*(Rest-0.0075)*620)))'; ...
    'BPA_40mm', 40,  0.122366088343, 10.4714342293, 2.02328814095, '6000'};

for s = 1:size(surf, 1)
    nm  = surf{s, 1};
    dia = surf{s, 2};
    a0  = surf{s, 3};
    a1  = surf{s, 4};
    a3  = surf{s, 5};
    fmaxDef = surf{s, 6};
    blk = [lib '/' nm];
    if getSimulinkBlockHandle(blk) > 0
        delete_block(blk);
    end
    add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 1280+180*(s-1) 160 1400+180*(s-1)]);
    delete_line(blk, 'In1/1', 'Out1/1');
    delete_block([blk '/In1']);
    delete_block([blk '/Out1']);
    B = @(n, p, v) add_block(p, [blk '/' n], 'Position', v);
    B('P',        'simulink/Sources/In1',                        [25 28 55 42]);   % pressure kPa
    B('L',        'simulink/Sources/In1',                        [25 143 55 157]); % actuator length m
    B('Pn',       'simulink/Math Operations/Gain',               [90 18 120 48]);  % P/620
    B('Pn620',    'simulink/Sources/Constant',                   [90 55 115 80]);  % 620
    B('tfit',     'simulink/Sources/Constant',                   [90 150 135 180]);% TendonL+2*FittingL
    B('Lsub',     'simulink/Math Operations/Sum',                [170 88 200 122]);% L - tfit
    B('con',      'simulink/Math Operations/Sum',                [235 78 265 112]);% Rest - Lmus
    B('Rest_c',   'simulink/Sources/Constant',                   [170 130 205 160]);% Rest
    B('relDiv',   'simulink/Math Operations/Gain',               [300 83 330 113]);% 1/KMAX
    B('KMAX_c',   'simulink/Sources/Constant',                   [300 130 350 160]);% (Rest-Kmax)/Rest
    B('negA1rel', 'simulink/Math Operations/Gain',               [370 55 400 85]); % -a1*rel
    B('exp1',     'simulink/Math Operations/Math Function',      [410 55 440 85]);
    B('a0g',      'simulink/Math Operations/Gain',               [470 55 500 85]); % a0*
    B('minusA0',  'simulink/Math Operations/Sum',                [530 50 560 90]); % - a0
    B('a0_c',     'simulink/Sources/Constant',                   [470 100 500 130]);
    B('rel2',     'simulink/Math Operations/Math Function',      [370 150 400 180]);
    B('negA3',    'simulink/Math Operations/Gain',               [430 150 460 180]);% -a3*rel^2
    B('exp2',     'simulink/Math Operations/Math Function',      [490 150 520 180]);
    B('pTerm',    'simulink/Math Operations/Product',            [550 120 580 150]);% Pn * exp2
    B('FnSum',    'simulink/Math Operations/Sum',                [620 70 650 110]);% a0 term + pTerm
    B('FmaxG',    'simulink/Math Operations/Gain',               [690 78 720 108]);% * Fmax
    B('relGate',  'simulink/Logic and Bit Operations/Compare To Constant', [690 140 730 170]);
    B('Fgate',    'simulink/Math Operations/Product',            [760 80 790 120]);% F * (rel<1)
    B('Fsat',     'simulink/Discontinuities/Saturation',         [820 85 850 115]);
    B('F',        'simulink/Sinks/Out1',                         [890 93 920 107]);

    set_param([blk '/Pn'],    'Gain', '1/620');
    set_param([blk '/Pn620'], 'Value', '620');
    set_param([blk '/tfit'],  'Value', 'TendonL + 2*FittingL');
    set_param([blk '/Lsub'],  'Inputs', '-+');
    set_param([blk '/con'],   'Inputs', '+-');   % Rest - (L - tendon - 2*fitting)
    set_param([blk '/Rest_c'],'Value', 'Rest');
    set_param([blk '/relDiv'],'Gain', '1/(Rest-Kmax)');  % rel = ((Rest-L)/Rest)/KMAX = (Rest-L)/(Rest-Kmax)
    set_param([blk '/KMAX_c'],'Value', '(Rest-Kmax)/Rest');
    set_param([blk '/negA1rel'], 'Gain', sprintf('-%.14g', a1));
    set_param([blk '/exp1'],  'Operator', 'exp');
    set_param([blk '/a0g'],   'Gain', sprintf('%.14g', a0));
    set_param([blk '/minusA0'],'Inputs', '+-');
    set_param([blk '/a0_c'],  'Value', sprintf('%.14g', a0));
    set_param([blk '/rel2'],  'Operator', 'square');
    set_param([blk '/negA3'], 'Gain', sprintf('-%.14g', a3));
    set_param([blk '/exp2'],  'Operator', 'exp');
    set_param([blk '/pTerm'], 'Inputs', '2');
    set_param([blk '/FnSum'], 'Inputs', '++');
    % Fmax baked into the internal gain: block params inside a masked
    % subsystem resolve mask params (Rest), while mask-param DEFAULT strings
    % that reference sibling params do NOT resolve at sim time (learned the
    % hard way 2026-09-10). Double-click into the block to override Fmax.
    set_param([blk '/FmaxG'], 'Gain', fmaxDef);
    set_param([blk '/relGate'], 'relop', '<', 'const', '1');
    set_param([blk '/Fgate'], 'Inputs', '2');
    set_param([blk '/Fsat'],  'UpperLimit', fmaxDef, 'LowerLimit', '0');

    L = @(a, b) add_line(blk, a, b, 'autorouting', 'on');
    L('P/1', 'Pn/1');
    L('Pn/1', 'pTerm/1');
    L('L/1', 'Lsub/2');
    L('tfit/1', 'Lsub/1');
    L('Lsub/1', 'con/2');
    L('Rest_c/1', 'con/1');
    L('con/1', 'relDiv/1');
    L('relDiv/1', 'negA1rel/1');
    L('negA1rel/1', 'exp1/1');
    L('exp1/1', 'a0g/1');
    L('a0g/1', 'minusA0/1');
    L('a0_c/1', 'minusA0/2');
    L('relDiv/1', 'rel2/1');
    L('rel2/1', 'negA3/1');
    L('negA3/1', 'exp2/1');
    L('exp2/1', 'pTerm/2');
    L('minusA0/1', 'FnSum/1');
    L('pTerm/1', 'FnSum/2');
    L('FnSum/1', 'FmaxG/1');
    L('FmaxG/1', 'Fgate/1');
    L('relDiv/1', 'relGate/1');
    L('relGate/1', 'Fgate/2');
    L('Fgate/1', 'Fsat/1');
    L('Fsat/1', 'F/1');

    m = Simulink.Mask.create(blk);
    m.Type = sprintf('SNS BPA %d mm Actuator', dia);
    m.Description = sprintf(['Real Festo BPA actuator, %d mm ID, Ben''s equations: ' ...
        'F = Fmax*[a0*(exp(-a1*rel)-1) + (P/620)*exp(-a3*rel^2)], rel = ' ...
        'contraction/KMAX, contraction = (Rest-(L-TendonL-2*FittingL))/Rest. ' ...
        'Zero force for rel>=1 and F<0 (festo4.m). Fmax = maxBPAforce(Rest,620) ' ...
        'is baked into the internal gain (edit inside the block to override); ' ...
        'Inputs: P [kPa], L [m]. Output: F [N].'], dia);
    m.Display = bpaIconCode(sprintf('BPA %d', dia));
    m.addParameter('Name', 'Rest', 'Type', 'edit', 'Prompt', 'Resting length Rest (m)', 'Value', '0.20');
    m.addParameter('Name', 'Kmax', 'Type', 'edit', 'Prompt', 'Fully-contracted length Kmax (m)', 'Value', '0.165');
    m.addParameter('Name', 'TendonL', 'Type', 'edit', 'Prompt', 'Tendon length (m)', 'Value', '0');
    m.addParameter('Name', 'FittingL', 'Type', 'edit', 'Prompt', 'One end-fitting length (m)', 'Value', '0');
    set_param(blk, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
        'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');
end

%% ---------------- BioMuscle (OpenSim Thelen2003-style) ----------------
blk = [lib '/BioMuscle'];
if getSimulinkBlockHandle(blk) > 0
    delete_block(blk);
end
add_block('simulink/Ports & Subsystems/Subsystem', blk, 'Position', [40 1820 160 1940]);
delete_line(blk, 'In1/1', 'Out1/1');
delete_block([blk '/In1']);
delete_block([blk '/Out1']);
B = @(n, p, v) add_block(p, [blk '/' n], 'Position', v);
B('Act',     'simulink/Sources/In1',                        [25 28 55 42]);   % activation 0..1
B('Lmt',     'simulink/Sources/In1',                        [25 93 55 107]);  % musculotendon length m
B('Vmt',     'simulink/Sources/In1',                        [25 158 55 172]); % musculotendon velocity m/s
% --- fiber length (rigid tendon) ---
B('lslack',  'simulink/Sources/Constant',                   [90 100 130 130]); % l_slack
B('lsub',    'simulink/Math Operations/Sum',                [170 88 200 122]);% Lmt - l_slack
B('cosA',    'simulink/Math Operations/Trigonometric Function', [170 150 200 180]);
B('Adeg',    'simulink/Sources/Constant',                   [100 185 135 215]); % alpha deg
B('d2r',     'simulink/Math Operations/Gain',               [170 190 200 220]);
B('lceDiv',  'simulink/Math Operations/Gain',               [230 90 260 120]); % /(cos a)
B('lmax',    'simulink/Discontinuities/Saturation',         [290 90 320 120]); % lower clamp 1e-3
B('lam',     'simulink/Math Operations/Gain',               [350 90 380 120]); % /l_opt
% --- active force-length ---
B('lm1',     'simulink/Math Operations/Sum',                [410 60 440 90]); % lambda-1
B('one_c',   'simulink/Sources/Constant',                   [410 100 435 120]);
B('flSq',    'simulink/Math Operations/Math Function',      [470 60 500 90]); % square
B('flG',     'simulink/Math Operations/Gain',               [530 60 560 90]); % -1/0.45
B('flExp',   'simulink/Math Operations/Math Function',      [650 60 680 90]); % exp
% --- force-velocity ---
B('cosV',    'simulink/Math Operations/Gain',               [230 158 260 188]); % v/cos(a)
B('lOpt_c',  'simulink/Sources/Constant',                   [290 200 330 230]); % l_opt
B('vDiv',    'simulink/Math Operations/Divide',             [290 155 320 190]); % /l_opt -> v_tilde
B('vNeg',    'simulink/Math Operations/Gain',               [350 158 380 188]); % s = -v_tilde
B('vSwitch', 'simulink/Signal Routing/Switch',              [440 155 470 205]); % s>=0 ? conc : ecc
B('concB',   'simulink/Sources/Constant',                   [400 155 430 185]); % b = A*v_max
B('bSum',    'simulink/Math Operations/Sum',                [470 250 500 280]); % s + b
B('bAb',     'simulink/Sources/Constant',                   [470 300 500 330]); % b*(1+A)
B('concDiv', 'simulink/Math Operations/Divide',             [540 255 570 290]);
B('concA',   'simulink/Math Operations/Sum',                [610 255 640 285]); % ... - A
B('hillA_c', 'simulink/Sources/Constant',                   [540 300 570 330]); % A
% eccentric: 1.4 - 0.4*exp(-v_tilde/0.2)
B('eccT',    'simulink/Math Operations/Gain',               [400 340 430 370]); % -v_tilde/0.2
B('eccExp',  'simulink/Math Operations/Math Function',      [460 340 490 370]);
B('eccA',    'simulink/Math Operations/Gain',               [520 340 550 370]); % -0.4*
B('eccSum',  'simulink/Math Operations/Sum',                [580 340 610 370]); % 1.4 + ...
B('c14_c',   'simulink/Sources/Constant',                   [520 385 555 405]);
% --- passive force-length ---
B('plm1',    'simulink/Math Operations/Sum',                [410 430 440 460]); % lambda-1
B('pG1',     'simulink/Math Operations/Gain',               [470 430 500 460]); % kpe/e_pas
B('pExp',    'simulink/Math Operations/Math Function',      [530 430 560 460]);
B('pSub',    'simulink/Math Operations/Sum',                [590 430 620 460]); % -1
B('pNorm',   'simulink/Math Operations/Gain',               [650 430 680 460]); % 1/(exp(kpe)-1)
B('pSat',    'simulink/Discontinuities/Saturation',         [710 430 740 460]); % [0, inf)
% --- combine: F = Fmax*(Act*fl*fv + fpas)*cos(a) ---
B('flfv',    'simulink/Math Operations/Product',            [740 70 770 100]);
B('actMul',  'simulink/Math Operations/Product',            [800 60 830 90]);
B('pAdd',    'simulink/Math Operations/Sum',                [860 70 890 110]);
B('FmaxG',   'simulink/Math Operations/Gain',               [920 78 950 108]);
B('cosMul',  'simulink/Math Operations/Product',            [980 80 1010 110]);
B('Fsat',    'simulink/Discontinuities/Saturation',         [1040 85 1070 115]);
B('F',       'simulink/Sinks/Out1',                         [1110 93 1140 107]);

set_param([blk '/lsub'],   'Inputs', '+-');   % Lmt - l_slack
set_param([blk '/lslack'], 'Value', 'l_slack');
set_param([blk '/vSwitch'], 'Criteria', 'u2 >= Threshold');
set_param([blk '/lmax'],   'UpperLimit', '10', 'LowerLimit', '1e-3');
set_param([blk '/lam'],    'Gain', '1/l_opt');
set_param([blk '/lm1'],    'Inputs', '+-');
set_param([blk '/flSq'],   'Operator', 'square');
set_param([blk '/flG'],    'Gain', '-1/0.45');
set_param([blk '/flExp'],  'Operator', 'exp');
set_param([blk '/cosA'],   'Operator', 'cos');
set_param([blk '/Adeg'],   'Value', 'alpha');
set_param([blk '/d2r'],    'Gain', 'pi/180');
set_param([blk '/lceDiv'], 'Gain', '1/cos(alpha*pi/180)');
set_param([blk '/cosV'],   'Gain', '1/cos(alpha*pi/180)');
set_param([blk '/lOpt_c'], 'Value', 'l_opt');
set_param([blk '/vNeg'],   'Gain', '-1');
set_param([blk '/concB'],  'Value', 'A_hill*v_max');          % b = A*v_max
set_param([blk '/bSum'],   'Inputs', '++');
set_param([blk '/bAb'],    'Value', 'A_hill*v_max*(1+A_hill)');
set_param([blk '/concA'],  'Inputs', '+-');
set_param([blk '/hillA_c'],'Value', 'A_hill');
set_param([blk '/eccT'],   'Gain', '-1/0.2');
set_param([blk '/eccExp'], 'Operator', 'exp');
set_param([blk '/eccA'],   'Gain', '-0.4');
set_param([blk '/eccSum'], 'Inputs', '++');
set_param([blk '/c14_c'],  'Value', '1.4');
set_param([blk '/pG1'],    'Gain', 'kpe/e_pas');
set_param([blk '/pSub'],   'Inputs', '+-');
set_param([blk '/one_c'],  'Value', '1');
set_param([blk '/pNorm'],  'Gain', '1/(exp(kpe)-1)');
set_param([blk '/pSat'],   'UpperLimit', 'inf', 'LowerLimit', '0');
set_param([blk '/flfv'],   'Inputs', '2');
set_param([blk '/actMul'], 'Inputs', '2');
set_param([blk '/pAdd'],   'Inputs', '++');
set_param([blk '/FmaxG'],  'Gain', 'Fmax');
set_param([blk '/cosMul'], 'Inputs', '2');
set_param([blk '/Fsat'],   'UpperLimit', '2*Fmax', 'LowerLimit', '0');

L = @(a, b) add_line(blk, a, b, 'autorouting', 'on');
L('Lmt/1', 'lsub/1');
L('lslack/1', 'lsub/2');
L('lsub/1', 'lceDiv/1');
L('lceDiv/1', 'lmax/1');
L('lmax/1', 'lam/1');
% active FL
L('lam/1', 'lm1/1');
L('one_c/1', 'lm1/2');
set_param([blk '/plm1'],    'Inputs', '+-');
L('lm1/1', 'flSq/1');
L('flSq/1', 'flG/1');
L('flG/1', 'flExp/1');
% velocity path
L('Adeg/1', 'd2r/1');
L('d2r/1', 'cosA/1');
L('Vmt/1', 'cosV/1');
L('cosV/1', 'vDiv/1');
L('lOpt_c/1', 'vDiv/2');
L('vDiv/1', 'vNeg/1');
L('vNeg/1', 'vSwitch/2');    % control: s >= 0 -> concentric branch
L('vNeg/1', 'eccT/1');
% concentric: fv = b(1+A)/(s+b) - A
L('concB/1', 'bSum/2');
L('vNeg/1', 'bSum/1');
L('bAb/1', 'concDiv/1');
L('bSum/1', 'concDiv/2');
L('concDiv/1', 'concA/1');
L('hillA_c/1', 'concA/2');
L('concA/1', 'vSwitch/1');   % data in: concentric fv
% eccentric: fv = 1.4 - 0.4*exp(-v/0.2)
L('eccT/1', 'eccExp/1');
L('eccExp/1', 'eccA/1');
L('eccA/1', 'eccSum/1');
L('c14_c/1', 'eccSum/2');
L('eccSum/1', 'vSwitch/3');  % data in: eccentric fv
% passive
L('lam/1', 'plm1/1');
L('one_c/1', 'plm1/2');
L('plm1/1', 'pG1/1');
L('pG1/1', 'pExp/1');
L('pExp/1', 'pSub/1');
L('one_c/1', 'pSub/2');
L('pSub/1', 'pNorm/1');
L('pNorm/1', 'pSat/1');
% combine
L('flExp/1', 'flfv/1');
L('vSwitch/1', 'flfv/2');
L('Act/1', 'actMul/1');
L('flfv/1', 'actMul/2');
L('actMul/1', 'pAdd/1');
L('pSat/1', 'pAdd/2');
L('pAdd/1', 'FmaxG/1');
L('FmaxG/1', 'cosMul/1');
L('cosA/1', 'cosMul/2');
L('cosMul/1', 'Fsat/1');
L('Fsat/1', 'F/1');

m = Simulink.Mask.create(blk);
m.Type = 'SNS Biological Muscle (Thelen-style)';
m.Description = ['OpenSim Thelen2003-style Hill-type biological muscle, rigid tendon: ' ...
    'F = Fmax*(Act*fl(lambda)*fv(v) + fpas(lambda))*cos(alpha). Active FL = ' ...
    'exp(-(lambda-1)^2/0.45); Hill force-velocity (A=0.25, v_max in l_opt/s); ' ...
    'exponential passive FL (kpe, e_pas); fiber length = (Lmt-l_slack)/cos(alpha). ' ...
    'Inputs: Act [0..1], Lmt [m], Vmt [m/s]. Output: F [N].'];
m.Display = bpaIconCode('BIO', [0.87 0.62 0.60]);
m.addParameter('Name', 'Fmax', 'Type', 'edit', 'Prompt', 'Max isometric force Fmax (N)', 'Value', '1000');
m.addParameter('Name', 'l_opt', 'Type', 'edit', 'Prompt', 'Optimal fiber length l_opt (m)', 'Value', '0.10');
m.addParameter('Name', 'l_slack', 'Type', 'edit', 'Prompt', 'Tendon slack length (m)', 'Value', '0.05');
m.addParameter('Name', 'alpha', 'Type', 'edit', 'Prompt', 'Pennation angle (deg)', 'Value', '0');
m.addParameter('Name', 'v_max', 'Type', 'edit', 'Prompt', 'Max shortening velocity (l_opt/s)', 'Value', '10');
m.addParameter('Name', 'A_hill', 'Type', 'edit', 'Prompt', 'Hill curvature A (a/Fmax)', 'Value', '0.25');
m.addParameter('Name', 'kpe', 'Type', 'edit', 'Prompt', 'Passive FL exponent kpe', 'Value', '4');
m.addParameter('Name', 'e_pas', 'Type', 'edit', 'Prompt', 'Passive FL strain at Fmax e_pas', 'Value', '0.6');
set_param(blk, 'MaskIconFrame', 'off', 'MaskIconUnits', 'autoscale', ...
    'MaskIconOpaque', 'on', 'MaskIconRotate', 'none');

save_system(lib);
set_param(lib, 'Lock', 'on');   % blocks can only be INSTANCED from a LOCKED library
fprintf('SNS_Library.slx: added BPA_10mm/BPA_20mm/BPA_40mm + BioMuscle (11 blocks total).\n');

%% ---------------- local functions ----------------
function s = bpaIconCode(lbl, fill)
    % Fusiform actuator with heavy outline. BPA = salmon; BIO = darker red.
    if nargin < 2, fill = [0.94 0.76 0.74]; end
    s = strjoin({ ...
        't_ = linspace(0, 2*pi, 73);' ...
        'patch(0.92*cos(t_), 0.42*sin(t_), [0 0 0]);' ...
        sprintf('patch(0.81*cos(t_), 0.33*sin(t_), [%.3g %.3g %.3g]);', fill) ...
        ['disp(''' lbl ''');'] ...
        }, newline);
end
