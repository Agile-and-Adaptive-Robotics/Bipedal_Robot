% Debug the flat-V Deng RG, take 3: handle-based line swap + logging on.
root = 'D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape';
cd(root); addpath(root); addpath(fullfile(root, 'demos'));
mdl = 'SNS_Deng_RG';
load_system(mdl);
set_param(mdl, 'StopTime', '2', 'SignalLogging', 'on', ...
          'SignalLoggingName', 'sigs');
% swap the I_stim->HC_ext line for a constant 10 nA source (handle-based)
lh = get_param([mdl '/I_stim'], 'LineHandles');
if lh.Outport(1) ~= -1
    delete_line(lh.Outport(1));
end
add_block('simulink/Sources/Constant', [mdl '/dbg10nA'], 'Value', '10', ...
          'Position', [400 400 430 430]);
add_line(mdl, 'dbg10nA/1', 'pair1_HC_ext/1', 'autorouting', 'on');

tgts = {[mdl '/pair1_HC_ext/Vint'], 'V'; ...
        [mdl '/pair1_HC_ext/GNa_g'], 'INa'; ...
        [mdl '/pair1_HC_ext/hint'], 'h'; ...
        [mdl '/pair1_HC_ext/isyn'], 'Isyn'; ...
        [mdl '/pair1_HC_ext/mInv'], 'm'; ...
        [mdl '/pair1_HC_ext/leak'], 'Ileak'};
for k = 1:size(tgts, 1)
    try
        ph = get_param(tgts{k, 1}, 'PortHandles');
        set_param(ph.Outport(1), 'DataLogging', 'on', ...
                  'DataLoggingNameMode', 'Custom', 'DataLoggingName', tgts{k, 2});
    catch e
        disp(['log fail ' tgts{k, 1} ': ' e.message]);
    end
end
out = sim(mdl, 'ReturnWorkspaceOutputs', 'on');
sigs = out.sigs;
for nm = {'V', 'INa', 'h', 'Isyn', 'm', 'Ileak'}
    try
        ts = sigs.get(nm{1}).Values;
        vv = ts.Data(:); tt = ts.Time(:);
        n = numel(tt);
        i1 = max(1, round(0.05 * n)); i2 = round(0.5 * n);
        fprintf('%-6s: t=%6.3f v=%12.5f | t=%6.3f v=%12.5f | t=%6.3f v=%12.5f\n', ...
            nm{1}, tt(i1), vv(i1), tt(i2), vv(i2), tt(end), vv(end));
    catch e
        disp([nm{1} ': ' e.message]);
    end
end
close_system(mdl, 0);
disp('debug3 done');
