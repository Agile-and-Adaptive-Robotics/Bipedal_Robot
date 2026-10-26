% Li twin smoke run: simulate SNS_Li_CPG 30 s under the scripted
% antiphase heel contacts, then report per-side stance-CPG activity and
% the inter-leg phase relationship.
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root); addpath(fullfile(root, 'demos'));
mdl = 'SNS_Li_CPG';
load_system(fullfile(root, 'demos', [mdl '.slx']));
set_param(mdl, 'SignalLogging', 'on', 'SignalLoggingName', 'sigs', ...
          'StopTime', '30');

% log every stance-CPG cell V and every MN V
cpgs = find_system(mdl, 'SearchDepth', 1, 'LookUnderMasks', 'all', ...
                   'Regexp', 'on', 'Name', 'cpg_.*');
fprintf('stance CPGs found: %d\n', numel(cpgs));
for k = 1:numel(cpgs)
    ph = get_param(cpgs{k}, 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', ['cpg' sprintf('%02d', k)]);
end
mns = find_system(mdl, 'SearchDepth', 1, 'LookUnderMasks', 'all', ...
                  'Regexp', 'on', 'Name', '(L|R).*MN|MN.*');
for k = 1:numel(mns)
    ph = get_param(mns{k}, 'PortHandles');
    set_param(ph.Outport(1), 'DataLogging', 'on', ...
              'DataLoggingNameMode', 'Custom', ...
              'DataLoggingName', ['mn' sprintf('%02d', k)]);
end

out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
sigs = out.sigs;
names = {};
try
    names = sigs.SignalInfo.BlockPath(1).OutputName; %#ok<SAGROW>
catch
end
% gather logged cpg signals generically
vs = {}; ts = {};
for k = 1:numel(cpgs)
    nm = ['cpg' sprintf('%02d', k)];
    try
        tsig = sigs.get(nm).Values;
        ts{end + 1} = tsig.Time(:); vs{end + 1} = tsig.Data(:); %#ok<AGROW>
        fprintf('%s (%s): V range [%.1f, %.1f] mV, std %.2f\n', nm, ...
            strrep(cpgs{k}, [mdl '/'], ''), min(vs{end}), max(vs{end}), std(vs{end}));
    catch e
        disp([nm ': ' e.message]);
    end
end
% inter-leg phase: correlate the first L-side vs first R-side stance CPG
if numel(vs) >= 2
    n = min(numel(vs{1}), numel(vs{2}));
    a = vs{1}(1:n) - mean(vs{1}(1:n));
    b = vs{2}(1:n) - mean(vs{2}(1:n));
    r = sum(a .* b) / sqrt(sum(a .^ 2) * sum(b .^ 2));
    fprintf('inter-leg stance correlation: %.3f (negative = antiphase)\n', r);
    % burst counts via prominence peaks
    for j = 1:2
        ipk = find(islocalmax(vs{j}, 'MinProminence', 5));
        fprintf('side %d: %d bursts, mean period %.3f s\n', j, numel(ipk), ...
            mean(diff(ts{j}(ipk))));
    end
end
close_system(mdl, 0);
disp('li smoke done');
