function goal2_run_baselines()
% GOAL 2 TASK A — run the committed SNS_Simscape baselines headless and
% record their metrics in one table (saved to goal2_baselines.mat).
%
% Runs, with the exact demo conventions from the committed runners:
%   1. sns_units_test_2n            (committed PASS bound: dev < 1e-2 mV,
%                                     committed value 6.07e-04 mV)
%   2. KneeReflexDemo   5 s         (rise 15.8 deg @ 0.05 s; settle mean
%                                     43.5 deg, band ~41.6-44.9 deg)
%   3. KneeReflexCircuit 5 s        (runnable 1:1 twin of the demo)
%   4. BPACPGLegDemo    10 s        (theta range 9.7..48.0 deg)
%   5. BeerCupReflexDemo ON/OFF     (max sag after 2 s: ON 2.0 / OFF 7.9 deg;
%                                     triceps A 0.40 -> 0.32 over the pour)
%
% Everything is READ-ONLY with respect to Ben's models: load_system + sim,
% no set_param on the demos beyond StopTime (which the runners also set).

rep = fileparts(mfilename('fullpath'));
d = rep;                                          % walk up to the repo root
while exist(fullfile(d, 'Code', 'Matlab', 'SNS_Simscape', 'SNS_Library.slx'), 'file') == 0
    dn = fileparts(d);
    if strcmp(dn, d), error('repo root with Code/Matlab/SNS_Simscape not found above %s', rep); end
    d = dn;
end
sns = fullfile(d, 'Code', 'Matlab', 'SNS_Simscape');
fprintf('MATLAB %s | SNS_Simscape at %s\n', version, sns);
addpath(sns);
addpath(fullfile(sns, 'demos'));
cd(fullfile(sns, 'demos'));

T = {};

%% 1) units test
try
    sns_units_test_2n();
    units = struct('pass', true, 'devA', NaN, 'devB', NaN);
    % re-extract the devs from the saved mat the test compares against
    S = load(fullfile(sns, 'results', 'units_ref_2n.mat'), 't', 'va', 'vb');
    units.note = 'PASS line printed above (tol 1e-2 mV)';
    T(end+1, :) = {'sns_units_test_2n', 'PASS (see printed line)', ''}; %#ok<SAGROW>
catch ME
    T(end+1, :) = {'sns_units_test_2n', 'FAIL', ME.message}; %#ok<SAGROW>
    rethrow(ME);
end

%% 2) KneeReflexDemo
k = run_knee('KneeReflexDemo', 5);
T(end+1, :) = {'KneeReflexDemo', sprintf( ...
    'theta(0.05s) %.2f deg | settle mean %.1f deg (range %.1f-%.1f, last 2 s) | final %.2f deg | A_ext %.2f A_flex %.2f (win means)', ...
    k.rise, k.settleMean, k.settleMin, k.settleMax, k.final, k.AeWin, k.AfWin), ''};

%% 3) KneeReflexCircuit (runnable twin)
c = run_knee('KneeReflexCircuit', 5);
T(end+1, :) = {'KneeReflexCircuit', sprintf( ...
    'theta(0.05s) %.2f deg | settle mean %.1f deg (range %.1f-%.1f, last 2 s) | final %.2f deg | A_ext %.2f A_flex %.2f (win means)', ...
    c.rise, c.settleMean, c.settleMin, c.settleMax, c.final, c.AeWin, c.AfWin), ...
    'twin of KneeReflexDemo'};

%% 4) BPACPGLegDemo
load_system('SNS_Library');
load_system('BPACPGLegDemo');
set_param('BPACPGLegDemo', 'StopTime', '10');
out = sim('BPACPGLegDemo');
th = out.log_th; Vex = out.log_V_RG_ext; Vfl = out.log_V_RG_flex;
t = th.Time;
dV = Vex.Data - Vfl.Data;
hyst = 0.1*max(dV) - 0.1*min(dV);
upLvl = max(dV) - hyst; dnLvl = min(dV) + hyst;
state = 1*(dV(1) >= 0); tsw = [];
for i = 2:numel(dV)
    if state == 0 && dV(i) > upLvl, state = 1; tsw(end+1) = t(i); %#ok<SAGROW>
    elseif state == 1 && dV(i) < dnLvl, state = 0; tsw(end+1) = t(i); %#ok<SAGROW>
    end
end
if numel(tsw) >= 2, period = 2*mean(diff(tsw)); else, period = NaN; end
cpg = struct('thetaMin', min(th.Data)*180/pi, 'thetaMax', max(th.Data)*180/pi, ...
    'nSw', numel(tsw), 'period', period);
T(end+1, :) = {'BPACPGLegDemo', sprintf( ...
    'theta range %.1f..%.1f deg | %d hyst switches, period ~%.3f s', ...
    cpg.thetaMin, cpg.thetaMax, cpg.nSw, cpg.period), ''};

%% 5) BeerCupReflexDemo ON/OFF
load_system('BeerCupReflexDemo');
runs = struct('label', {'reflex ON', 'reflex OFF'}, 'k', {1, 0});
for r = 1:2
    assignin('base', 'kReflex', runs(r).k);
    out = sim('BeerCupReflexDemo');
    runs(r).th = out.log_th;   runs(r).mCup = out.log_mCup;
    runs(r).A = out.log_A_bi;  runs(r).At = out.log_A_tri;
    runs(r).F = out.log_F_bi;  runs(r).Ft = out.log_F_tri;
    w = runs(r).th.Time > 2;
    runs(r).sagAfter2 = max(abs(runs(r).th.Data(w)))*180/pi;
    runs(r).final = runs(r).th.Data(end)*180/pi;
end
assignin('base', 'kReflex', 1);
T(end+1, :) = {'BeerCupReflexDemo', sprintf( ...
    ['max sag after 2 s: ON %.2f deg / OFF %.2f deg | final ON %+.2f / OFF %+.2f deg | ' ...
     'A_bi %.3f->%.3f, A_tri %.3f->%.3f (ON)'], ...
    runs(1).sagAfter2, runs(2).sagAfter2, runs(1).final, runs(2).final, ...
    runs(1).A.Data(1), runs(1).A.Data(end), runs(1).At.Data(1), runs(1).At.Data(end)), ''};

%% table + save
fprintf('\n=== GOAL 2 BASELINE TABLE ===\n');
for i = 1:size(T, 1)
    fprintf('%-22s | %s %s\n', T{i, 1}, T{i, 2}, T{i, 3});
end
save(fullfile(rep, 'goal2_baselines.mat'), 'T', 'k', 'c', 'cpg', 'runs');
fprintf('saved %s\n', fullfile(rep, 'goal2_baselines.mat'));
end

function r = run_knee(mdl, tf)
load_system('SNS_Library');
load_system(mdl);
set_param(mdl, 'StopTime', num2str(tf));
out = sim(mdl);
th = out.log_th; Ae = out.log_A_ext; Af = out.log_A_flex;
w = th.Time > tf - 2;                     % last-2-s settle window
r.rise = interp1(th.Time, th.Data, 0.05, 'pchip')*180/pi;
r.settleMean = mean(th.Data(w))*180/pi;
r.settleMin = min(th.Data(w))*180/pi;
r.settleMax = max(th.Data(w))*180/pi;
r.final = th.Data(end)*180/pi;
r.AeWin = mean(Ae.Data(Ae.Time > tf - 2));
r.AfWin = mean(Af.Data(Af.Time > tf - 2));
end
