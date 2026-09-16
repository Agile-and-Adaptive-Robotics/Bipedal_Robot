function run_bridge_tests()
% Bridge prove-it tests a/b/c on the sensor-patched simbridge model.
%   a: fire one muscle   (constant ctrl, vas_med_r = 1) -> force + knee respond
%   b: two clock rates   (0.001 s ctrl source -> Rate Transition -> 0.005 s plant)
%   c: sensor readback   (a's trajectory vs Python mujoco 2.3.7 ground truth)
% Models are written to matlab/models, logs to mujoco_bridge/logs.

here = fileparts(mfilename('fullpath'));
repo = fullfile(here, '..', '..', '..', '..', '..');
xml  = strrep(fullfile(repo, 'Solid_Models', 'OpenSim', 'Gait2392_Robotbody', ...
                'mjc', 'gait2392_simbody', 'gait2392_simbody_cvt3_simbridge.xml'), ...
              '\', '/');
modelsDir = fullfile(here, 'models');
if ~isfolder(modelsDir), mkdir(modelsDir); end

NU = 92; VAS = 28;          % vas_med_r actuator index (0-based -> MATLAB 29)

%% ---- test a: fire one muscle --------------------------------------------
uConst = zeros(NU, 1); uConst(VAS + 1) = 1;
a = run_case(fullfile(modelsDir, 'bridge_test_a.slx'), xml, 'const', uConst, [], 0.4);

frcA = a.vals.vas_med_r_frc; kneeA = a.vals.knee_r_pos;
pass_a = abs(frcA(2)) < 1e-6 ...                       % before first step: 0
         && max(abs(frcA)) > 500 ...                   % muscle force responds
         && max(abs(frcA([5:end]))) > 500 ...
         && abs(rad2deg(kneeA(end) - kneeA(1))) > 5;   % knee actually moves
fprintf(['TEST a (fire one muscle): %s | frc 0->%.1f N (max |frc| %.1f N), ' ...
         'knee %.2f->%.2f deg\n'], tf(pass_a), frcA(end), ...
        max(abs(frcA)), rad2deg(kneeA(1)), rad2deg(kneeA(end)));

%% ---- test b: two clock rates --------------------------------------------
n = 401; tv = (0:n-1)' * 0.001;
uTs = zeros(n, NU); uTs(:, VAS + 1) = max(0, sin(2*pi*tv));  % half-rect 1 Hz
b = run_case(fullfile(modelsDir, 'bridge_test_b.slx'), xml, 'ts', [], uTs, 0.4);

dtB = mean(diff(b.time));
kneeB = b.vals.knee_r_pos;
% plumbing: clean 5 ms sensor grid; dynamics: muscle force clearly responds
% to the half-rectified sine (peak 1.0 at t=0.25 s) and the knee moves
pass_b = all(isfinite(b.vals.vas_med_r_frc)) && abs(dtB - 0.005) < 1e-9 ...
         && numel(b.time) == 81 ...
         && max(abs(b.vals.vas_med_r_frc)) > 50 ...
         && rad2deg(max(kneeB) - min(kneeB)) > 5;
fprintf(['TEST b (two clock rates): %s | sensor dt %.4f s (%d samples), ' ...
         'max |frc| %.1f N, knee range %.1f deg\n'], tf(pass_b), dtB, ...
        numel(b.time), max(abs(b.vals.vas_med_r_frc)), ...
        rad2deg(max(kneeB) - min(kneeB)));

%% ---- test c: sensor readback vs Python ground truth ---------------------
gt = load(fullfile(here, '..', 'logs', 'gt_const.mat'), 'sens');
gts = gt.sens;                                     % [80 x 5], step k = row k
fields = fieldnames(a.vals); %#ok<NASGU>
S = cell2mat(struct2cell(a.vals).');               % [81 x 5] Simulink rows
% Simulink row k+1 (time 0.005k) should equal Python step k
dev0 = sum(max(abs(S(2:end, :) - gts), [], 1));    % aligned as derived
dev1 = sum(max(abs(S(1:end-1, :) - gts), [], 1));  % one step later
[~, align] = min([dev0 dev1]);
if align == 1, C = S(2:end, :); else, C = S(1:end-1, :); end
dmax = max(abs(C - gts), [], 1);
pass_c = dmax(1) < deg2rad(0.5) && dmax(2) < 0.05 && dmax(3) < 1e-3 ...
         && dmax(4) < 0.1 && dmax(5) / 2005 < 0.01;
fprintf(['TEST c (sensor readback): %s | alignment shift %d, per-sensor max ' ...
         'abs dev: knee %.2e rad, knee_vel %.2e, len %.2e m, vel %.2e m/s, ' ...
         'frc %.2e N\n'], tf(pass_c), align - 1, dmax);

save(fullfile(here, '..', 'logs', 'bridge_tests_abc.mat'), 'a', 'b', ...
     'dmax', 'pass_a', 'pass_b', 'pass_c');
end

function r = run_case(slxPath, xml, srcType, uConst, uTs, stopTime)
% Build + run one bridge test model; return parsed sensor bus.
mdl = 'bridge_test_tmp';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);

add_block('mjLib/MuJoCo Plant', [mdl '/Plant'], ...
    'xmlFileRel', char(xml), ...
    'renderingType', 'None', 'rgbOutOption', 'off', 'depthOutOption', 'off', ...
    'Position', [260 100 420 240]);

if strcmp(srcType, 'const')
    add_block('simulink/Sources/Constant', [mdl '/u'], 'Value', 'uConst', ...
              'Position', [80 130 140 170]);
    add_line(mdl, 'u/1', 'Plant/1', 'autorouting', 'on');
else
    add_block('simulink/Sources/From Workspace', [mdl '/u'], ...
              'VariableName', 'uTs', 'SampleTime', '0.001', ...
              'Position', [60 130 140 170]);
    add_block('simulink/Signal Attributes/Rate Transition', [mdl '/RT'], ...
              'Position', [190 130 220 170]);
    add_line(mdl, 'u/1', 'RT/1', 'autorouting', 'on');
    add_line(mdl, 'RT/1', 'Plant/1', 'autorouting', 'on');
end

add_block('simulink/Sinks/Terminator', [mdl '/sensT'], 'Position', [520 130 550 170]);
add_line(mdl, 'Plant/1', 'sensT/1', 'autorouting', 'on');
pph = get_param([mdl '/Plant'], 'PortHandles');
set_param(pph.Outport(1), 'DataLogging', 'on', 'DataLoggingNameMode', 'Custom', ...
          'DataLoggingName', 'sb');

set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', ...
          'FixedStep', 'auto', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');
assignin('base', 'uConst', uConst);
assignin('base', 'uTs', timeseries(uTs, (0:size(uTs, 1)-1)' * 0.001));

out = sim(mdl, 'StopTime', num2str(stopTime), 'ReturnWorkspaceOutputs', 'on');
r = parse_bus(out.sigs);
save_system(mdl, slxPath);
close_system(mdl, 0);
end

function r = parse_bus(ds)
% Line logging of a bus -> Dataset element whose Values is a struct of
% timeseries (one per sensor field).
el = ds{1};
V = el.Values;
fn = fieldnames(V);
r.time = V.(fn{1}).Time(:);
r.vals = struct();
for k = 1:numel(fn)
    r.vals.(fn{k}) = V.(fn{k}).Data(:);
end
end

function s = tf(p)
if p, s = 'PASS'; else, s = 'FAIL'; end
end
