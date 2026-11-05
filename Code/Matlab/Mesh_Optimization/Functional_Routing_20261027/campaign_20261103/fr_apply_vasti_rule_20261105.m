function fr_apply_vasti_rule_20261105()
% Combined vasti remains a parallel 20 mm design, as requested by Ben.
C=fr_setup();C.exactSegmentDistance=true;C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);root=fullfile(C.camp,'audit_20261105','reoptimization');
if isfile(fullfile(root,'finalization_complete.json'))
 done=jsondecode(fileread(fullfile(root,'finalization_complete.json')));
 if isfield(done,'vas_int_20mm_only')&&done.vas_int_20mm_only;return;end
 input=fullfile(root,'final_results.mat');
else
 input=fullfile(root,'coordinated_progress.mat');
 P=load(input,'branch','start');assert(strcmp(P.branch,'compliant')&&P.start==25,'Coordinated round incomplete');
end
S=load(input,'results');results=S.results;
save(fullfile(root,'before_vasti_diameter_rule.mat'),'results','-v7');
idx=find(cellfun(@(r)strcmp(r.compliant.caseDef.name,'vas_int'),results));assert(isscalar(idx));
for br={'rigid','compliant'}
 branch=br{1};chosen=results{idx}.(branch);base=chosen.caseDef;base.diaStr='20';base.diaNum=20;
 envs=cellfun(@(r)r.(branch).envelope,results,'UniformOutput',false);envs(idx)=[];
 opts=struct('quiet',true,'envelopes',{envs});x=fr_packX(base,chosen.design);
 bestF=inf;bestCS=base;timer=tic;
 counts=unique([base.seed.nBPA,3,4,8,12,16],'stable');
 for nb=counts
  for gain=unique([max(base.txConfig.gain),1,2,3],'stable')
   families={'none'};if gain>1;families={'S','M','L'};end
   for family=families
    if toc(timer)>30;break;end
    cs=base;cs.seed.nBPA=nb;cs.hasTackle=gain>1;cs.txConfig.nPulleyBPA=[nb,1];
    cs.txConfig.gain=[gain,1];cs.txConfig.tackleLineParts=[gain,1];cs.txConfig.tackleFamily=family{1};
    f=score(x,C,cs,kin,branch,opts);
    if f<bestF;bestF=f;bestCS=cs;end
   end
   if toc(timer)>30;break;end
  end
  if toc(timer)>30;break;end
 end
 cs=bestCS;free=find(cs.ub-cs.lb>1e-10);fixed=x;
 fun=@(z)score(expand(z,fixed,free),C,cs,kin,branch,opts);
 po=optimoptions('patternsearch','Display','off','UseParallel',false,'MaxFunctionEvaluations',300,'MaxTime',40,'MeshTolerance',1e-4,'StepTolerance',1e-4);
 [z,f]=patternsearch(fun,x(free),[],[],[],[],cs.lb(free),cs.ub(free),[],po);
 if f<bestF;x=expand(z,fixed,free);bestF=f;end
 [~,rep]=fr_objective(C,cs,kin,x,branch,opts);
 chosen.caseDef=cs;chosen.design=fr_unpackX(cs,x);chosen.xBest=x;chosen.fBest=bestF;chosen.rep=rep;
 chosen.envelope=fr_envelope(C,cs,chosen.design,kin,cs.passes{1}.qStructs,branch);
 chosen.diameterPolicy='Combined vasti: 20 mm only. Cam remains unresolved.';
 results{idx}.(branch)=chosen;
 fprintf('VASTI 20MM %s n%d G%g meet %.3f violation %.4g\n',branch,cs.seed.nBPA,max(cs.txConfig.gain),rep.fractionPosesMet,rep.violation);
end
fr_validate_repaired_layouts_20261105(results);
end
function x=expand(z,x,idx)
x(idx)=z;
end
function f=score(x,C,cs,kin,branch,opts)
try;f=fr_objective(C,cs,kin,x,branch,opts);if ~isreal(f)||~isfinite(f);f=1e12;end
catch;f=1e12;end
end
