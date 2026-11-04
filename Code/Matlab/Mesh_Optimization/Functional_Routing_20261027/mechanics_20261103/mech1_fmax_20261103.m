function mech1_fmax_20261103()
%MECH1_FMAX_20261103 MechanCheck row 1: Fmax from maxBPAforce.m at 620 kPa
% and 450 mm resting length, 10 mm and 20 mm, vs the handoff expectation
% (~441 N / ~1314 N in HARDWARE_AND_ZCODE_HANDOFF.md "Load sizing that the
% optimizer must supply"). Also prints the exact equation evaluated and
% the bundle arithmetic used by MechanCheck row 5.
%
% Output: console log (captured by mech1_fmax_20261103.log).

repo = 'D:/GitHub/Bipedal_Robot';
addpath(genpath(fullfile(repo, 'Code', 'Matlab')));
addpath(fullfile(repo, 'Code', 'Matlab', 'Mesh_Optimization'), '-end');
addpath(fullfile(repo, 'Testing_Data', '2022_02_Festo'), '-end');

fprintf('=== MechanCheck 1: maxBPAforce at 620 kPa, L0 = 450 mm ===\n')
fprintf('Source: Code/Matlab/Functions/maxBPAforce.m\n')
fprintf('  10 mm branch (lines 39-48): z = y*(a1*atan(a2*(x-0.0075)*y)), ')
fprintf('a1 = 0.4895 N/kPa, a2 = 0.03068 1/(kPa*m)\n')
fprintf('  20 mm branch (lines 50-59): z = y*(a1*atan(a2*(x-0.0075)*y)), ')
fprintf('a1 = 1.4877 N/kPa, a2 = 0.0248 1/(kPa*m)\n\n')

L0 = 0.45;      % m, resting length
P  = 620;       % kPa

% 10 mm: evaluate the equation term by term, then via the function.
arg10 = 0.03068 * (L0 - 0.0075) * P;
F10_hand = P * (0.4895 * atan(arg10));
F10 = maxBPAforce(L0, '10', 620);
fprintf('10 mm:  atan argument = 0.03068*(0.45-0.0075)*620 = %.6f\n', arg10)
fprintf('        atan(arg) = %.6f rad\n', atan(arg10))
fprintf('        hand-substituted: 620*0.4895*atan(arg) = %.4f N\n', F10_hand)
fprintf('        maxBPAforce(0.45,10,620)            = %.4f N\n', F10)
fprintf('        handoff expects ~441 N -> abs diff %.4f N (%.3f%%)\n\n', ...
    abs(F10 - 441), 100 * abs(F10 - 441) / 441)

% 20 mm
arg20 = 0.0248 * (L0 - 0.0075) * P;
F20_hand = P * (1.4877 * atan(arg20));
F20 = maxBPAforce(L0, '20', 620);
fprintf('20 mm:  atan argument = 0.0248*(0.45-0.0075)*620 = %.6f\n', arg20)
fprintf('        atan(arg) = %.6f rad\n', atan(arg20))
fprintf('        hand-substituted: 620*1.4877*atan(arg) = %.4f N\n', F20_hand)
fprintf('        maxBPAforce(0.45,20,620)            = %.4f N\n', F20)
fprintf('        handoff expects ~1314 N -> abs diff %.4f N (%.3f%%)\n\n', ...
    abs(F20 - 1314), 100 * abs(F20 - 1314) / 1314)

% Consistency between hand-substitution and the function call.
fprintf('function-vs-hand-equation max abs diff: %.3e N (must be 0)\n', ...
    max(abs(F10 - F10_hand), abs(F20 - F20_hand)))

% Zero-contraction caveat (handoff): these are zero-contraction capacities.
fprintf(['NOTE: these are ZERO-contraction capacities (handoff line 32); ' ...
    'the sweep force falls with contraction via festo4.\n'])

% Bundle arithmetic reused by MechanCheck 5 (handoff example, line 36):
fprintf('\nHandoff example inputs (row 5):\n')
fprintf('  3 x 20 mm bundle at 450 mm Fmax = %.2f N\n', 3 * F20)
fprintf('  G=2 ideal shared cable  T = 3*F20/2 = %.2f N (handoff ~1971 N)\n', 3 * F20 / 2)
fprintf('  180-deg reversing sheave R = 2*T     = %.2f N (handoff ~3942 N)\n', 2 * (3 * F20 / 2))
fprintf('=== MechanCheck 1 done ===\n')
end
