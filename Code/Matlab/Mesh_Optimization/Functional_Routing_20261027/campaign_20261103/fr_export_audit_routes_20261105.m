function fr_export_audit_routes_20261105()
% Read existing checkpoints without changing their recorded results.
C=fr_setup(); kin=fr_kinematics(C.osimStock);
outdir=fullfile(C.camp,'audit_20261105');if ~isfolder(outdir);mkdir(outdir);end
files=dir(fullfile(C.camp,'run_1','checkpoints','*.mat'));rows={};
for i=1:numel(files)
 if contains(files(i).name,'refine');continue;end
 S=load(fullfile(files(i).folder,files(i).name));
 if ~isfield(S,'design') || ~isfield(S,'rep');continue;end
 sig=strsplit(S.cfgSig,'|');ov=containers.Map('KeyType','char','ValueType','any');
 ov(S.csName)=struct('diaStr',sig{1},'nBPA',str2double(sig{2}),'family',sig{3},'gain',str2double(sig{4}));
 cases=fr_cases27(C,kin,ov);cs=cases(strcmp({cases.name},S.csName));
 ids=find(cs.rowUse{1});loc=S.rep.loc;
 assert(size(loc,1)==numel(ids),'Row count differs for %s',S.csName);
 assigned=cs.rowFrame{1}(ids);own=cs.rowOwning(ids);ownloc=loc;
 for r=1:numel(ids)
  ex=cs.rowExpr{1}{ids(r)};
  if ~isempty(ex);ownloc(r,:)=(ex(1:3,1:3)'*(loc(r,:)'-ex(1:3,4)))';end
 end
 rows{end+1}=struct('group',S.csName,'branch',S.branch,'checkpoint',fullfile(files(i).folder,files(i).name), ...
 'status',S.status,'cfgSig',S.cfgSig,'nominal_diameter_mm',str2double(sig{1}),'nBPA',str2double(sig{2}), ...
 'family',sig{3},'gain',str2double(sig{4}),'design',S.design,'recordedLoc',loc,'assignedBodies',{assigned}, ...
 'owningBodies',{own},'owningLoc',ownloc,'seedLoc',cs.seedLoc, ...
 'rowIndices',ids,'cross',cs.crossFor{1},'fractionMet',S.rep.fractionPosesMet,'recordedViolation',S.violation, ...
 'constraints',{S.rep.cinfo},'terms',S.rep.terms,'tauRobot',S.rep.tauPrimary,'tauTarget',S.rep.tauTarget, ...
 'primaryDof',cs.passes{1}.sweptCoord,'qDefault',kin.qDefault, ...
 'seedSignatureMatches',S.nX==cs.nX && abs(S.xSeedSum-sum(cs.xSeed))<1e-9,'envelope',S.envelope);
 fprintf('Export %s %s n%d G%g\n',S.csName,S.branch,str2double(sig{2}),str2double(sig{4}));
end
fid=fopen(fullfile(outdir,'recorded_routes.json'),'w');assert(fid>0);fwrite(fid,jsonencode(rows));fclose(fid);
fprintf('Exported %d checkpoint routes to %s\n',numel(rows),outdir);
end
