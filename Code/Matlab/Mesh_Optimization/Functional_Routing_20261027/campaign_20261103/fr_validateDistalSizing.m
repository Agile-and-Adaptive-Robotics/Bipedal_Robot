function fr_validateDistalSizing(outputDir,budget)
% Validate the stronger EHL seed and the new distal-insertion design variables.
if nargin<2; budget=100; end
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
stage=fullfile(C.camp,'run_20261105_distal_stages');
branch='compliant';
name='ext_hal_compliant';
folder=fullfile(stage,'ankle_mounts_600');
S=load(fullfile(folder,[name '.mat']),'result');
seed=S.result.chosen;
if ~isfolder(outputDir); mkdir(outputDir); end
file=fullfile(outputDir,[name '.mat']);
if isfile(file); return; end
[original,baseDetail]=fr_evalCoupledAnkle(C,kin,seed.caseDef,seed.design,branch,21);
candidate=seed;
candidate.caseDef.diaNum=20;
candidate.caseDef.diaStr='20';
toe=fr_toeShiftPatches(C,kin,candidate.caseDef);
[moved,shift]=fr_shiftToeRoute(candidate.caseDef,toe,8);
assert(all(shift>0) && isequal(moved.rowOwning,candidate.caseDef.rowOwning));
[unchanged,unchangedShift]=fr_shiftToeRoute(candidate.caseDef,toe,0);
assert(all(unchangedShift==0) && isequal(unchanged.seedLoc,candidate.caseDef.seedLoc));
G=fr_ankleCoupledSetup(C,kin,moved,11);
d=fr_refitMountLengths(C,kin,moved,candidate.design,branch,G,.135);
probe=fr_evalCoupledAnkle(C,kin,moved,d,branch,11,G);
fprintf('Distal probe: moves %s mm, BPA rest %.1f mm, tendon %.1f mm, no force %.1f%%\n', ...
    mat2str(1000*shift',3),1000*d.L0,1000*d.tendon,100*mean(probe.noForce));
[chosen,info]=fr_surrogateDistalCase(C,kin,candidate,branch,budget,true,[],true);
[E,D]=fr_evalCoupledAnkle(C,kin,chosen.caseDef,chosen.design,branch,21);
S=load(fullfile(stage,'ankle_guide_arrow','ext_dig_compliant.mat'),'result');
edl=S.result.chosen;
[other,outer]=fr_evalCoupledAnkle(C,kin,edl.caseDef,edl.design,branch,21);
assert(isequal(other.poses,E.poses));
layer=fr_compareAnkleLayers(C,kin,E.poses,outer.sleeves,D.sleeves, ...
    edl.caseDef.fittingLength,chosen.caseDef.fittingLength);
layer.layerPass=layer.overlapPairs>0 && min([layer.orderingMargin_m,layer.rigidMargin_m])>=-1e-4;
record=jsondecode(fileread(fullfile(folder,[name '.json'])));
record=rmfield(record,intersect(fieldnames(record),{'surrogateEvaluations','initialScore','finalScore'}));
record.legacySweep=false;
record.diameter_mm=20;
record.baselineCoupled=original;
record.coupled=E;
record.layerWithEDL=layer;
record.coupled.localPass=E.localPass && other.localPass && layer.layerPass;
record.localConstraintsPass=record.coupled.localPass;
record.foundLengthFit=~any(E.noForce);
record.rest_mm=1000*chosen.design.L0;
record.tendon_mm=1000*chosen.design.tendon;
record.toeShift_mm=1000*(chosen.caseDef.seedLoc(toe.rows,:)-toe.points)*toe.distal';
record.trialEvaluations=info.evaluations;
result=struct('chosen',chosen);
save(file,'result','record','info','D','baseDetail','probe','shift','-v7');
fid=fopen(fullfile(outputDir,[name '.json']),'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
fr_refreshCoupledRoutePlots(outputDir);
fprintf('EHL D20: %.1f mm rest, %.1f mm tendon, toe shifts %s mm, no force %.1f%%, bone %+.2f mm, EDL fitting %+.2f mm\n', ...
    record.rest_mm,record.tendon_mm,mat2str(record.toeShift_mm',3),100*mean(E.noForce), ...
    1000*E.boneMargin_m,1000*layer.rigidMargin_m);
end
