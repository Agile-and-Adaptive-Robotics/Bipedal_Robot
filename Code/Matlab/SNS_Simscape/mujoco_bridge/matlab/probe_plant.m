function probe_plant()
% Mask parameters of the MuJoCo Plant library block.
load_system('mjLib');
blk = 'mjLib/MuJoCo Plant';
dp = get_param(blk, 'DialogParameters');
fn = fieldnames(dp);
fprintf('mask/dialog params:\n');
for k = 1:numel(fn)
    try
        v = get_param(blk, fn{k});
        fprintf('  %-22s = %s\n', fn{k}, char(string(v)));
    catch
        fprintf('  %-22s (no value)\n', fn{k});
    end
end
fprintf('\nS-Function name: %s\n', get_param([blk '/S-Function'], 'FunctionName'));
sf = get_param([blk '/S-Function'], 'DialogParameters');
sfn = fieldnames(sf);
for k = 1:numel(sfn)
    try
        v = get_param([blk '/S-Function'], sfn{k});
        if ischar(v) || isstring(v)
            fprintf('  sfparam %-18s = %s\n', sfn{k}, char(string(v)));
        end
    catch
    end
end
bdclose('mjLib');
end
