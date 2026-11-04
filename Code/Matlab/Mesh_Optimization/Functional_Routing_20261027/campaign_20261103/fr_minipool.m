function ok = fr_minipool()
%FR_MINIPOOL Parfor-plumbing test for the campaign driver paths.
% Pool(10) + broadcast (CASES with anonymous handles, kin with
% containers.Map fields) + struct-array slicing + a real compliant seed
% eval per first-layer group inside workers. NO checkpoints touched.
setenv('OMP_NUM_THREADS', '1');
setenv('MKL_NUM_THREADS', '1');
C = fr_setup();
kin = fr_kinematics(C.osimStock);
CASES = fr_cases27(C, kin);
tc = jsondecode(fileread(fullfile(C.targetsDir, 'target_contract.json')));
[layerIds, regionName] = fr_layerIds(tc, CASES);
if isempty(gcp('nocreate'))
    pool = parpool('Processes', C.parallelMaxWorkers);
end
fprintf('pool: %d workers\n', pool.NumWorkers);
ids = layerIds{1};
fprintf('layer 1 (%s): %s\n', regionName{1}, mat2str(ids));
byId = containers.Map('KeyType', 'double', 'ValueType', 'double');
for c = 1:numel(CASES)
    byId(CASES(c).groupMapId) = c;
end
templ = struct('name', '', 'f', NaN, 'wall', NaN, 'viol', NaN);
outs = repmat(templ, 1, numel(ids));
parfor ii = 1:numel(ids)
    cs = CASES(byId(ids(ii)));
    [fR, repR] = fr_objective(C, cs, kin, cs.xSeed, 'compliant', ...
        struct('quiet', true));
    o = templ;
    o.name = cs.name;
    o.f = fR;
    o.wall = repR.wallSeconds;
    o.viol = repR.violation;
    outs(ii) = o;
end
ok = true;
for ii = 1:numel(ids)
    fprintf('  %-9s compliant seed f %+.3f viol %.2e (%.2f s)\n', ...
        outs(ii).name, outs(ii).f, outs(ii).viol, outs(ii).wall);
    ok = ok && isfinite(outs(ii).f);
end
if ok
    fprintf('MINIPOOL PASSED\n');
else
    fprintf('MINIPOOL FAILED\n');
end
end
