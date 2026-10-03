%% goal2_regen_r2025a.m — regenerate SNS_Library + the four demo models as
% R2025a-native on THIS machine (easteregg2, MATLAB R2025a).
%
% WHY (README_SNS_Simscape.md "2026-09-20 Rybak restyle"): the committed
% SNS_Library.slx + demos\*.slx were saved by the laptop (R2025b-native), and
% R2025a REFUSES to load R2025b models. The README's documented route for
% R2025a machines is to re-run the committed builders locally:
%   sns_build_library.m -> sns_build_actuators.m  (library)
%   sns_build_demo.m -> sns_build_circuit_view.m -> sns_build_cpg_demo.m
%   -> sns_build_beer_demo.m                       (demos)
% The builders are the source of truth (the .slx are their deterministic
% outputs). The git tree will show these .slx as modified — that is the
% expected side effect of running on R2025a, listed in the goal-2 report.

rep = fileparts(mfilename('fullpath'));
d = rep;
while exist(fullfile(d, 'Code', 'Matlab', 'SNS_Simscape', 'SNS_Library.slx'), 'file') == 0
    dn = fileparts(d);
    if strcmp(dn, d), error('repo root not found above %s', rep); end
    d = dn;
end
sns = fullfile(d, 'Code', 'Matlab', 'SNS_Simscape');
fprintf('R2025a regen | MATLAB %s | %s\n', version, sns);
cd(sns);
addpath(sns); addpath(fullfile(sns, 'demos'));

%% library (documented order)
sns_build_library;
sns_build_actuators;
sns_test_actuators;      % committed validation, expect PASS ~4.9e-10 N

%% demos (build order matters: CPG demo loads KneeReflexDemo for KneeModel)
cd(fullfile(sns, 'demos'));
sns_build_demo;          % -> KneeReflexDemo.slx
sns_build_circuit_view;  % -> KneeReflexCircuit.slx (loads KneeReflexDemo)
sns_build_cpg_demo;      % -> BPACPGLegDemo.slx
sns_build_beer_demo;     % -> BeerCupReflexDemo.slx

fprintf('R2025a REGEN DONE\n');
