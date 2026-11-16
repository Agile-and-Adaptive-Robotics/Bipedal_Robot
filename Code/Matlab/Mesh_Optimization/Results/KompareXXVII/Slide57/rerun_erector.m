function rerun_erector(out)
% Fresh OpenSim erector-spinae search for the Slide 57 schematic.
here=fileparts(mfilename('fullpath'));
root=fileparts(fileparts(fileparts(here)));
camp=fullfile(root,'Functional_Routing_20261027','campaign_20261103');
addpath(camp);
C=fr_setup(); kin=fr_kinematics(C.osimStock);
cases=fr_cases27(C,kin,[],struct());
seed=cases(strcmp({cases.name},'ercspn'));
assert(isscalar(seed));
seed.diaStr='20'; seed.diaNum=20; seed.xi3Mode='perPose';
seed.bundlePacking='rows';
if ~isfolder(out); mkdir(out); end
pool=parcluster('local');
pool.JobStorageLocation=fullfile(out,'workers');
if ~isfolder(pool.JobStorageLocation); mkdir(pool.JobStorageLocation); end
parpool(pool,min(8,pool.NumWorkers));
[n,g,b]=ndgrid(1:4,1:3,1:2);
branches={'rigid','compliant'};
trials=cell(numel(n),1);
parfor j=1:numel(n)
    trials{j}=one(C,kin,seed,n(j),g(j),branches{b(j)},out);
end
for k=1:2
    candidates=trials(b(:)==k);
    rank=cell2mat(cellfun(@(x)x.rank,candidates,'UniformOutput',false));
    [~,i]=sortrows(rank);
    record=candidates{i(1)}.record;
    record.openSimSeed_m=seed.seedLoc;
    record.seedModel=C.osimStock;
    record.torso_to_pelvis=fr_chainT(kin,kin.qDefault,'torso','pelvis');
    record.femur_to_pelvis=fr_chainT(kin,kin.qDefault,'femur_r','pelvis');
    write(fullfile(out,['erector_' branches{k} '.json']),record);
    fprintf('Selected %s: %d BPAs, gain %d, %d topology sheaves, %.1f%% target, %.1f%% no force.\n', ...
        record.branch,record.count,record.gain,record.sheaves,100*record.fractionMet,100*record.noForce);
end
end

function r=one(C,kin,cs,n,g,branch,out)
file=fullfile(out,sprintf('erector_n%d_g%d_%s.mat',n,g,branch));
if isfile(file); S=load(file,'r'); r=S.r; return; end
rng(n+10*g+sum(double(branch)),'twister');
reeving=fr_travelReeving(g);
cs.seed.nBPA=n;
cs.txConfig.nPulleyBPA=[1 1];
cs.txConfig.gain=[g 1];
cs.txConfig.tackleLineParts=[g 1];
cs.hasTackle=g>1;
if cs.hasTackle; cs.txConfig.tackleFamily='L'; else; cs.txConfig.tackleFamily='none'; end
cs.bundleSpread=1.9;
opts=struct('quiet',true,'envelopes',{{}});
ids=find(cs.ub>cs.lb+1e-10);
fun=@(v)fr_objective(C,cs,kin,replace(cs.xSeed,ids,v),branch,opts);
so=optimoptions('surrogateopt','Display','off','UseParallel',false, ...
    'MaxFunctionEvaluations',100,'InitialPoints',cs.xSeed(ids)');
[v,f]=surrogateopt(fun,cs.lb(ids),cs.ub(ids),so);
po=optimoptions('patternsearch','Display','off','MaxFunctionEvaluations',90,'MeshTolerance',1e-4);
[vp,fp]=patternsearch(fun,v,[],[],[],[],cs.lb(ids),cs.ub(ids),[],po);
if fp<f; v=vp; end
x=replace(cs.xSeed,ids,v); d=fr_unpackX(cs,x);
dense=fr_densePilotCase(C,kin,cs);
[~,rep]=fr_objective(C,dense,kin,x,branch,opts);
[tau,active,m]=fr_primaryRoutingSweep(C,kin,dense,d,branch);
h=dense.target.tau(:);
hit=mean(active & tau.*sign(h)>=abs(h)-1e-6);
shortfall=sum(max(0,abs(h)-tau.*sign(h)))/sum(abs(h));
[p,frames]=fr_locForPass(C,cs,d,kin,1);
record=struct('branch',branch,'count',n,'gain',g,'sheaves',n*reeving.sheaveCount, ...
    'sheavesPerLine',reeving.sheaveCount,'points_m',p,'frames',{frames}, ...
    'localPass',rep.feasible,'constraints',rep.cinfo,'fractionMet',hit,'noForce',mean(~active), ...
    'signedCapacity',sum(min(abs(h),max(0,tau.*sign(h))))/sum(abs(h)), ...
    'rest_mm',1000*d.L0,'tendon_mm',1000*d.tendon,'robot_Nm',tau,'human_Nm',h, ...
    'angles_deg',rad2deg(dense.passes{1}.th1),'peakTendon_N',max(m.Ftendon,[],'omitnan'));
r=struct('record',record,'rank',[~rep.feasible,~all(active),-hit,shortfall,n,n*reeving.sheaveCount]);
temp=[tempname(out) '.mat']; save(temp,'r','cs','d','rep'); movefile(temp,file,'f');
fprintf('%s %d BPAs G%d: %.1f%% target, local pass %d.\n',branch,n,g,100*hit,rep.feasible);
end

function x=replace(x,i,v)
x(i)=v;
end

function write(file,x)
temp=[tempname(fileparts(file)) '.json'];
f=fopen(temp,'w'); assert(f>=0);
fprintf(f,'%s\n',jsonencode(x)); fclose(f); movefile(temp,file,'f');
end
