% Fresh-eyes review check (read-only): verbatim 4-route check + all-27 build.
camp = 'D:\GitHub\Bipedal_Robot\tmp\paper_bridge\capability_table_20261022';
addpath(camp);
addpath('D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization');

% 1) verbatim check (prints its own PASS lines; errors on mismatch)
check_verbatim_20261023();

% 2) all-27 build
specs = biPulleySpecsFromMasters_20261023();
fprintf('nSpecs=%d\n', numel(specs));
assert(numel(specs) == 27, 'expected 27');
nBad = 0;
for k = 1:numel(specs)
    sp = specs(k);
    okid = sp.actuatorId == k;
    okpts = all(isfinite(sp.points(:)));
    okfmax = isfinite(sp.fmax) && sp.fmax > 0;
    okprim = strcmp(sp.kin.dofs{sp.kin.primarySlot}, sp.kin.primaryDof);
    if ~(okid && okpts && okfmax && okprim); nBad = nBad + 1; end
    fprintf('%2d %-12s prim=%-17s slot=%d hinges=%s+%s cross=[%d %d] %s\n', ...
        sp.actuatorId, sp.name, sp.kin.primaryDof, sp.kin.primarySlot, ...
        sp.kin.dofs{1}, sp.kin.dofs{2}, sp.cross(1), sp.cross(2), ...
        ternary(okid&&okpts&&okfmax&&okprim, 'OK', 'BAD'));
end
assert(nBad == 0, '%d bad specs', nBad);
fprintf('BUILD27 CHECK PASS\n');

function s = ternary(c, a, b)
if c; s = a; else; s = b; end
end
