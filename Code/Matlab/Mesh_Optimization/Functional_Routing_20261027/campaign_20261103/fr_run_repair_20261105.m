function fr_run_repair_20261105()
% Separate, checkpointed bounded reoptimization; old runs are never replaced.
C=fr_setup();kin=fr_kinematics(C.osimStock);root=fullfile(C.camp,'audit_20261105','reoptimization');
if ~isfolder(root);mkdir(root);end
if isfile(fullfile(root,'finalization_complete.json'))&&isfile(fullfile(root,'final_results.mat'))
 fprintf('Completed repaired campaign retained. Export/inspect final_results.mat; use a new output directory for a new campaign.\n');return
end
diary(fullfile(root,'master.log'));clean=onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('REPAIR RUN START %s %s\n',char(datetime('now')),getenv('COMPUTERNAME'));
files=dir(fullfile(C.camp,'run_1','checkpoints','*_compliant.mat'));ov=containers.Map('KeyType','char','ValueType','any');
for i=1:numel(files)
 S=load(fullfile(files(i).folder,files(i).name));sig=strsplit(S.cfgSig,'|');
 ov(S.csName)=struct('diaStr',sig{1},'nBPA',str2double(sig{2}),'family',sig{3},'gain',str2double(sig{4}));
end
cases=fr_cases27(C,kin,ov);tc=jsondecode(fileread(fullfile(C.targetsDir,'target_contract.json')));
[layers,~]=fr_layerIds(tc,cases);map=containers.Map({cases.name},num2cell(1:numel(cases)));
idmap=containers.Map([cases.groupMapId],num2cell(1:numel(cases)));
% Three independent regional fronts, one group per region per wave. Each
% region's anatomical order is preserved; other current designs are checked
% at simultaneous poses during the final all-group validation.
pool=gcp('nocreate');
if isempty(pool)
 pc=parcluster('local');pc.JobStorageLocation=fullfile(root,'parallel_jobs');if ~isfolder(pc.JobStorageLocation);mkdir(pc.JobStorageLocation);end
 pool=parpool(pc,min(3,C.parallelMaxWorkers));
end
fprintf('Workers %d. Per branch: discrete seed screen, surrogate <=60 eval/12s, pattern <=300 eval/40s.\n',pool.NumWorkers);
envelopes={};results=cell(1,numel(cases));
for wave=1:max(cellfun(@numel,layers))
 indices=[];
 for region=1:numel(layers)
  if wave<=numel(layers{region});indices(end+1)=idmap(layers{region}(wave));end %#ok<AGROW>
 end
 batch=cell(1,numel(indices));prior=envelopes;
 parfor k=1:numel(indices)
  idx=indices(k);cs=cases(idx);folder=fullfile(root,cs.name);if ~isfolder(folder);mkdir(folder);end
  result=repairOne(C,kin,cs,folder,prior);
  batch{k}=result;
 end
 for k=1:numel(indices)
  idx=indices(k);results{idx}=batch{k};chosen=batch{k}.compliant;
  envelopes{end+1}=chosen.envelope; %#ok<AGROW>
  fprintf('WAVE %d %s: rigid meet %.3f viol %.4g; compliant meet %.3f viol %.4g\n',wave,cases(idx).name, ...
   batch{k}.rigid.rep.fractionPosesMet,batch{k}.rigid.rep.violation,chosen.rep.fractionPosesMet,chosen.rep.violation);
 end
 save(fullfile(root,'progress.mat'),'results','wave','-v7');
end
% Final global status belongs to fr_finalize_repaired_20261105, which uses
% matching layouts for each branch and the final verified collision model.
save(fullfile(root,'final_results.mat'),'results','-v7');
fprintf('INITIAL SWEEP COMPLETE %s: %d branches checkpointed; global validation pending.\n',char(datetime('now')),2*numel(results));
end

function result=repairOne(C,kin,base,folder,prior)
for br={'rigid','compliant'}
 branch=br{1};destination=fullfile(folder,[branch '.mat']);
 if isfile(destination);old=load(destination,'chosen');result.(branch)=old.chosen;continue;end
 diary(fullfile(folder,[branch '.log']));fprintf('START %s %s %s\n',base.name,branch,char(datetime('now')));
 rng(C.rngSeed+base.groupMapId+double(strcmp(branch,'compliant'))*1000);
 S=load(fullfile(C.camp,'run_1','checkpoints',[base.name '_' branch '.mat']));
n0=base.seed.nBPA;counts=unique([n0,min(8,n0+1),min(8,max(4,2*n0))]);
if any(strcmp(base.name,{'vas_int','glut_max3','add_mag3'}));counts=unique([counts,8,12,16]);end
dias={'10','20'};gains=unique([max(base.txConfig.gain),1,2,3]);bestF=inf;bestX=[];csBest=[];
 opts=struct('quiet',true,'envelopes',{prior});
 for dd=1:numel(dias)
  for nb=counts
   for gain=gains
    families={'none'};if gain>1;families={'S','M','L'};end
    for family=families
    cs=base;cs.diaStr=dias{dd};cs.diaNum=str2double(cs.diaStr);cs.seed.nBPA=nb;
    cs.hasTackle=gain>1;cs.txConfig.nPulleyBPA=[nb,1];cs.txConfig.gain=[gain,1];cs.txConfig.tackleLineParts=[gain,1];
    cs.txConfig.tackleFamily=family{1};
    seeds={cs.xSeed,min(cs.ub,max(cs.lb,fr_packX(cs,S.design)))};
    for seed=seeds
     f=safeObjective(C,cs,kin,seed{1},branch,opts);
     if f<bestF;bestF=f;bestX=seed{1};csBest=cs;end
    end
    end
   end
  end
 end
 assert(~isempty(bestX),'No evaluable configuration for %s',base.name);
 cs=csBest;free=find(cs.ub-cs.lb>1e-10);fixed=bestX;
 fun=@(z)safeObjective(C,cs,kin,expand(z,fixed,free),branch,opts);
 fprintf('SCREEN %s %s dia%s n%d G%g objective %.5g\n',cs.name,branch,cs.diaStr,cs.seed.nBPA,max(cs.txConfig.gain),bestF);
 if exist('surrogateopt','file')==2 && numel(free)>1
  try
   so=optimoptions('surrogateopt','Display','off','UseParallel',false,'MaxFunctionEvaluations',60,'MaxTime',12,'InitialPoints',bestX(free)');
   [z,fs]=surrogateopt(fun,cs.lb(free),cs.ub(free),so);
   if fs<bestF;bestX=expand(z,fixed,free);bestF=fs;end
  catch ME;fprintf('SURROGATE %s\n',ME.message);end
 end
 po=optimoptions('patternsearch','Display','off','UseParallel',false,'MaxFunctionEvaluations',300,'MaxTime',40,'MeshTolerance',1e-4,'StepTolerance',1e-4);
 fixed=bestX;fun=@(z)safeObjective(C,cs,kin,expand(z,fixed,free),branch,opts);
 [z,fp,exitflag,output]=patternsearch(fun,bestX(free),[],[],[],[],cs.lb(free),cs.ub(free),[],po);
 if fp<bestF;bestX=expand(z,fixed,free);bestF=fp;end
 [~,rep]=fr_objective(C,cs,kin,bestX,branch,opts);design=fr_unpackX(cs,bestX);
 env=fr_envelope(C,cs,design,kin,cs.passes{1}.qStructs,branch);
 chosen=struct('caseDef',cs,'design',design,'rep',rep,'envelope',env,'exitflag',exitflag,'solverOutput',output,'xBest',bestX,'fBest',bestF);
 save(destination,'chosen','-v7');result.(branch)=chosen;
 fprintf('DONE %s %s meet %.3f viol %.4g %s\n',cs.name,branch,rep.fractionPosesMet,rep.violation,char(datetime('now')));diary('off');
end
end
function x=expand(z,x,idx)
x(idx)=z;
end
function f=safeObjective(C,cs,kin,x,branch,opts)
try
 f=fr_objective(C,cs,kin,x,branch,opts);if ~isscalar(f)||~isreal(f)||~isfinite(f);f=1e12;end
catch
 f=1e12;
end
end
