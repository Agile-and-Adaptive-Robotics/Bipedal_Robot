function fr_repair_gate_20261105()
C=fr_setup();kin=fr_kinematics(C.osimStock);outdir=fullfile(C.camp,'audit_20261105');
% Analytic check: an X-axis coordinate must use the X wrench component,
% and a translating origin must contribute force dot translation velocity.
cs.passes={struct('sweptCoord','q','qStructs',{{struct('q',.3)}},'T1',@testTransform,'T2',@(~,~) eye(4))};
m.cell.Torque_ins=reshape([7,11,13],1,3,1,1);m.cell.F_ins=reshape([17,19,23],1,3,1,1);
tau=fr_coordinateTorque(cs,kin,m,1,1);
expected=7+dot([17;19;23],fr_rotAxis([1;0;0],.3)'*[.2;0;0]);assert(abs(tau-expected)<1e-7);
fprintf('PASS generalized wrench projection %.9g (expected %.9g)\n',tau,expected);
files=dir(fullfile(C.camp,'run_1','checkpoints','*_compliant.mat'));ov=containers.Map('KeyType','char','ValueType','any');
for i=1:numel(files)
 S=load(fullfile(files(i).folder,files(i).name));ss=strsplit(S.cfgSig,'|');
 ov(S.csName)=struct('diaStr',ss{1},'nBPA',str2double(ss{2}),'family',ss{3},'gain',str2double(ss{4}));
end
cases=fr_cases27(C,kin,ov);rows={};
for i=1:numel(cases)
 cs=cases(i);
 for branch={'rigid','compliant'}
  f=fullfile(C.camp,'run_1','checkpoints',[cs.name '_' branch{1} '.mat']);if ~isfile(f);continue;end
  S=load(f);
  try
   opts=struct('quiet',true,'envelopes',{{}});
   [~,rep]=fr_objective(C,cs,kin,fr_packX(cs,S.design),branch{1},opts);
   qg=struct('passIdx',1,'qList',{cs.passes{1}.qStructs});[m,meta]=fr_mechBuild(C,cs,S.design,kin,branch{1},qg);
   prim=fr_sleevePrimitives(C,cs,S.design,kin,m,meta);
   count=numel(prim{1});nsec=1+double(strcmp(cs.spanKind,'knee')&&S.design.seriesSplit<.999);
   assert(count==cs.seed.nBPA*nsec,'Missing parallel vessels');
   lengths=arrayfun(@(x)norm(x.p2-x.p1),prim{1});
   assert(all(lengths>1e-6),'Collapsed sleeve');
   [tau,rot,trans]=fr_coordinateTorque(cs,kin,m,1,cs.primaryCrossing);
   row=struct('group',cs.name,'branch',branch{1},'oldFractionMet',S.rep.fractionPosesMet, ...
    'correctedFractionMet',rep.fractionPosesMet,'oldViolation',S.violation,'correctedViolation',rep.violation, ...
    'topViolator',rep.topViolator,'nBPA',cs.seed.nBPA,'physicalSections',count, ...
    'oldStrainMin',min(m.legacyRouteStrain),'oldStrainMax',max(m.legacyRouteStrain), ...
    'physicalStrainMin',min(m.strain_p),'physicalStrainMax',max(m.strain_p), ...
    'legacyTau',S.rep.tauPrimary,'correctedTau',tau,'rotationContribution',rot,'translationContribution',trans, ...
    'constraints',{rep.cinfo});rows{end+1}=row;
   fprintf('AUDIT %s %s old meet %.3g corrected %.3g violation %.4g (%s)\n',cs.name,branch{1},S.rep.fractionPosesMet,rep.fractionPosesMet,rep.violation,rep.topViolator);
  catch ME
   fprintf(2,'AUDIT ERROR %s %s %s\n',cs.name,branch{1},getReport(ME,'extended'));
  end
 end
end
fid=fopen(fullfile(outdir,'corrected_checkpoint_audit.json'),'w');fwrite(fid,jsonencode(rows));fclose(fid);
assert(numel(rows)==54,'fr_repair_gate:Incomplete','Only %d of 54 branches checked',numel(rows));
fprintf('REPAIR GATE PASSED: 54 branches reevaluated, all physical sleeve counts/lengths checked.\n');
end
function T=testTransform(~,q)
T=eye(4);T(1:3,1:3)=fr_rotAxis([1;0;0],q.q);T(1,4)=.2*q.q;
end
