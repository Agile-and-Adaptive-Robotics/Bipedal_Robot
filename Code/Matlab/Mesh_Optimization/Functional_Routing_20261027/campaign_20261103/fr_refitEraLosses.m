function fr_refitEraLosses(out,source)
% Give the current routes the same fixed-attachment length search as 2019.
C=fr_setup(); kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
file=fullfile(source,'era_comparison_by_group.csv');
import=detectImportOptions(file,'TextType','string');
import=setvartype(import,{'group','era','branch'},'string');
era=readtable(file,import);
A=sortrows(era(era.era=="2019",:),{'branch','group'});
B=sortrows(era(era.era=="current",:),{'branch','group'});
ids=find(B.torque_margin_met_pct<A.torque_margin_met_pct-1e-8);
upper=fullfile(C.camp,'run_20261106_remaining20');
lower=fullfile(C.camp,'run_20261106_ankle20_independent');
cluster=parcluster('local');
cluster.JobStorageLocation=fullfile(out,'era_refit_workers');
if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
parpool(cluster,min(5,cluster.NumWorkers));
parts=cell(numel(ids),1);
parfor k=1:numel(ids)
    j=ids(k);
    parts{k}=one(C,kin,B(j,:),A(j,:),upper,lower,out);
end
writetable(vertcat(parts{:}),fullfile(out,'matched_length_refit.csv'));
end

function R=one(C,kin,current,old,upper,lower,out)
name=sprintf('%s_%s',current.group,current.branch);
checkpoint=fullfile(out,[name '_matched_refit.mat']);
if isfile(checkpoint); S=load(checkpoint,'R'); R=S.R; return; end
file=fullfile(upper,'paired_refit',[name '.mat']);
if ~isfile(file); file=fullfile(upper,[name '.mat']); end
if isfile(file)
    S=load(file,'chosen'); chosen=S.chosen;
else
    S=load(fullfile(lower,[name '.mat']),'result'); chosen=S.result.chosen;
end
cs=fr_densePilotCase(C,kin,chosen.caseDef);
d=chosen.design;
h=cs.target.tau(:);
branch=char(current.branch);
opts=optimoptions('patternsearch','Display','off','MaxFunctionEvaluations',65, ...
    'InitialMeshSize',.03,'MeshTolerance',1e-4);
lb=[max(.012,.30*d.L0),max(.004,.30*d.tendon)];
ub=[max(d.L0*1.55,lb(1)+.02),max(d.tendon*2.5,.12)];
x0=[d.L0,d.tendon];
best=score(x0,C,kin,cs,d,branch,h); selected=x0;
for k=[.75 1 1.2]
    x=min(ub,max(lb,[k*d.L0,d.tendon+(1-k)*d.L0]));
    [x,f]=patternsearch(@(x)score(x,C,kin,cs,d,branch,h),x,[],[],[],[],lb,ub,[],opts);
    if f<best; best=f; selected=x; end
end
d.L0=selected(1); d.tendon=selected(2);
[t,a]=fr_primaryRoutingSweep(C,kin,cs,d,branch);
R=table(current.group,current.branch,old.torque_margin_met_pct,current.torque_margin_met_pct, ...
    100*mean(t.*sign(h)>=abs(h)-1e-6 & a),current.active_rom_pct,100*mean(a), ...
    100*d.L0,100*d.tendon, ...
    'VariableNames',{'group','branch','old_met_pct','current_met_pct','refit_met_pct', ...
    'current_active_pct','refit_active_pct','rest_cm','tendon_cm'});
save(checkpoint,'R','t','a','h','d');
fprintf('%s old %.1f%% current %.1f%% refit %.1f%%\n',name,R.old_met_pct,R.current_met_pct,R.refit_met_pct);
end

function f=score(x,C,kin,cs,d,branch,h)
d.L0=x(1); d.tendon=x(2);
[t,a]=fr_primaryRoutingSweep(C,kin,cs,d,branch);
deficit=max(0,abs(h)-sign(h).*t);
f=sum(deficit)/max(sum(abs(h)),eps)+.4*mean(deficit>1e-6)+2*mean(~a);
end
