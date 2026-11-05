function fr_finalize_repaired_20261105()
% One coordinated refinement round and independent, branch-consistent checks.
C=fr_setup();kin=fr_kinematics(C.osimStock);root=fullfile(C.camp,'audit_20261105','reoptimization');
if isfile(fullfile(root,'finalization_complete.json'));return;end
C.exactSegmentDistance=true;C.fullFittingEnvelopes=true;
% Preliminary search retains its constructed C. Restart workers so final
% refinement and checks uniformly use the verified distance/fitting model.
pool=gcp('nocreate');if ~isempty(pool);delete(pool);end
clear fr_objective fr_obstacles fr_retainedEnvelope fr_envelope
pc=parcluster('local');pc.JobStorageLocation=fullfile(root,'parallel_jobs');parpool(pc,3);
fprintf('FINAL GEOMETRY: exact segment minima; full fitting extents; fresh three-worker pool.\n');
S=load(fullfile(root,'final_results.mat'),'results');results=S.results;
save(fullfile(root,'first_pass_results.mat'),'results','-v7');
for br={'rigid','compliant'}
 branch=br{1};envs=cellfun(@(r)r.(branch).envelope,results,'UniformOutput',false);
 for start=1:3:numel(results)
  ids=start:min(start+2,numel(results));batch=cell(size(ids));prior=envs;
  parfor k=1:numel(ids)
   idx=ids(k);chosen=results{idx}.(branch);cs=chosen.caseDef;
   other=prior;other(idx)=[];opts=struct('quiet',true,'envelopes',{other});
   x=fr_packX(cs,chosen.design);free=find(cs.ub-cs.lb>1e-10);
   [f0,r0]=fr_objective(C,cs,kin,x,branch,opts);
   % Re-screen nearby discrete choices under the FINAL collision model;
   % preliminary ranking must not freeze a now-invalid sizing choice.
   base=cs;bestCS=cs;screenTimer=tic;
   counts=unique([base.seed.nBPA,max(1,base.seed.nBPA-1),base.seed.nBPA+1],'stable');
   dias=unique({base.diaStr,'10','20'},'stable');
   gains=unique([max(base.txConfig.gain),1,2,3],'stable');
   for nb=counts
    for dd=1:numel(dias)
     for gain=gains
      families={'none'};if gain>1;families=unique({base.txConfig.tackleFamily,'S','M','L'},'stable');families(strcmp(families,'none'))=[];end
      for familyIndex=1:numel(families)
       if toc(screenTimer)>30;break;end
       trial=base;trial.seed.nBPA=nb;trial.diaStr=dias{dd};trial.diaNum=str2double(trial.diaStr);
       trial.hasTackle=gain>1;trial.txConfig.nPulleyBPA=[nb,1];trial.txConfig.gain=[gain,1];
       trial.txConfig.tackleLineParts=[gain,1];trial.txConfig.tackleFamily=families{familyIndex};
       try
        [fv,rv]=fr_objective(C,trial,kin,x,branch,opts);
        if isreal(fv)&&isfinite(fv)&&fv<f0;f0=fv;r0=rv;bestCS=trial;end
       catch
        % Rejected candidate retains no force/feasibility credit.
       end
      end
      if toc(screenTimer)>30;break;end
     end
     if toc(screenTimer)>30;break;end
    end
    if toc(screenTimer)>30;break;end
   end
   cs=bestCS;chosen.caseDef=cs;chosen.fBest=f0;chosen.xBest=x;
   fprintf('FINAL SCREEN %s %s dia%s n%d G%g %.1fs\n',cs.name,branch,cs.diaStr,cs.seed.nBPA,max(cs.txConfig.gain),toc(screenTimer));
   fun=@(z)score(z,x,free,C,cs,kin,branch,opts);
   po=optimoptions('patternsearch','Display','off','UseParallel',false, ...
    'MaxFunctionEvaluations',200,'MaxTime',20,'MeshTolerance',1e-4,'StepTolerance',1e-4);
   [z,f]=patternsearch(fun,x(free),[],[],[],[],cs.lb(free),cs.ub(free),[],po);
   if isfinite(f)&&f<f0
    x(free)=z;[~,r0]=fr_objective(C,cs,kin,x,branch,opts);
    chosen.design=fr_unpackX(cs,x);chosen.xBest=x;chosen.fBest=f;
   end
   chosen.rep=r0;chosen.globalRep=r0;
   chosen.envelope=fr_envelope(C,cs,chosen.design,kin,cs.passes{1}.qStructs,branch);
   chosen.refinement='Final-model local discrete screen <=30s, then all-other-group refinement <=200 eval /20s';
   batch{k}=chosen;
  end
  for k=1:numel(ids)
   idx=ids(k);results{idx}.(branch)=batch{k};envs{idx}=batch{k}.envelope;
   fprintf('COORDINATED %s %s meet %.3f violation %.4g\n',batch{k}.caseDef.name,branch,batch{k}.rep.fractionPosesMet,batch{k}.rep.violation);
  end
  save(fullfile(root,'coordinated_progress.mat'),'results','branch','start','-v7');
 end
end
fr_apply_vasti_rule_20261105();
end
function f=score(z,x,free,C,cs,kin,branch,opts)
x(free)=z;
try;f=fr_objective(C,cs,kin,x,branch,opts);if ~isfinite(f)||~isreal(f);f=1e12;end
catch;f=1e12;end
end
