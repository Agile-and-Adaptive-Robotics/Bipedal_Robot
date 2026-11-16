function fr_checkPairedPackagingGrid
% Replay all rectangular combinations as a paired list, including their
% common maximum-span reference used by the pulley closure.
C=fr_setup(); kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
f=fullfile(C.camp,'run_20261106_ankle_assembled','right_side','assembled_validation.mat');
S=load(f,'routes');
r=S.routes{find(cellfun(@(v)strcmp(v.group,'ext_dig'),S.routes),1)};
G=fr_resolveAnkleNative(C,kin,r.caseDef,r.p,3);
original=fr_mechBuild(C,G.cs,r.design,kin,'compliant',G.qg);
qg=rmfield(G.qg,'qList2');
qg.qList=G.cs.passes{1}.qStructs;
qg.pairedPoseGrid=true;
paired=fr_mechBuild(C,G.cs,r.design,kin,'compliant',qg);
assert(max(abs(original.cell.MuscleLength(:)-paired.cell.MuscleLength(:)))<1e-10);
assert(max(abs(original.cell.Ftendon(:)-paired.cell.Ftendon(:)))<1e-8);
assert(isequal(original.cell.PulleySlack(:),paired.cell.PulleySlack(:)));
assert(isequal(original.cell.PulleyInfeasible(:),paired.cell.PulleyInfeasible(:)));
fprintf('Paired-pose grid matches rectangular BiPulley length, tendon force, and closure flags.\n');
end
