function fr_exportRectusSizingRoute(outputDir)
% Export the saved right rectus femoris candidate in native body frames.
C = fr_setup();
kin = fr_kinematics(C.osimStock);
branch = 'compliant';
S = load(fullfile(C.camp,'audit_20261105','reoptimization','rect_fem', ...
    ['coordinated_' branch '.mat']),'chosen');
c = S.chosen;
[pts,bodies,ids] = fr_locForPass(C,c.caseDef,c.design,kin,2);
assert(isequal(bodies,c.caseDef.rowOwning(ids)), ...
    'Expected the saved hip pass to use native body frames');
if ~isfolder(outputDir); mkdir(outputDir); end
record = struct('points_m',pts,'bodies',{bodies},'stock',C.osimStock,'robot',C.osimRobot, ...
    'branch',branch,'endFitting_m',c.caseDef.fittingLength,'pressure_kPa',C.pressure);
fid = fopen(fullfile(outputDir,'rectus_route.json'),'w');
assert(fid>=0);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(record));
end
