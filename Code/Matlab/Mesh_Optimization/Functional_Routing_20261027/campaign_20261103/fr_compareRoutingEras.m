function fr_compareRoutingEras(out)
% Replay both eras on the same primary coordinates and signed human curves.
C=fr_setup();
kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
upper=fullfile(C.camp,'run_20261106_remaining20');
lower=fullfile(C.camp,'run_20261106_ankle20_independent');
old=fr_cases27(C,kin,[],struct('routeFromRobotbody',true));
if ~isfolder(out); mkdir(out); end
pool=gcp('nocreate');
if isempty(pool)
    cluster=parcluster('local');
    cluster.JobStorageLocation=fullfile(out,'compare_workers');
    if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
    parpool(cluster,min(8,cluster.NumWorkers));
end
parts=cell(numel(old),1);
parfor j=1:numel(old)
    parts{j}=compareGroup(C,kin,upper,lower,out,old(j));
end
T=vertcat(parts{:});
writetable(T,fullfile(out,'era_comparison_by_group.csv'));
summary=groupsummary(T,{'era','branch'},'mean',{'torque_margin_met_pct','active_rom_pct'});
writetable(summary,fullfile(out,'era_comparison.csv'));
disp(summary);
end

function T=compareGroup(C,kin,upper,lower,out,ref)
file=fullfile(out,['era_' ref.name '.mat']);
if isfile(file); s=load(file,'T'); T=s.T; return; end
T=table;
for b=["rigid","compliant"]
    key=sprintf('%s_%s.mat',ref.name,b);
    f=fullfile(upper,'paired_refit',key);
    if ~isfile(f); f=fullfile(upper,key); end
    if isfile(f)
        S=load(f,'chosen'); cur=S.chosen;
    else
        S=load(fullfile(lower,key),'result'); cur=S.result.chosen;
    end
    cs=fr_densePilotCase(C,kin,cur.caseDef);
    [tau,active]=fr_primaryRoutingSweep(C,kin,cs,cur.design,char(b));
    target=cs.target.tau(:);
    T=[T; row(ref.name,b,"current",tau,active,target)]; %#ok<AGROW>
    baseline=ref;
    baseline.diaNum=cs.diaNum; baseline.diaStr=cs.diaStr;
    baseline.seed.nBPA=cs.seed.nBPA;
    baseline.txConfig=cs.txConfig;
    baseline.hasTackle=cs.hasTackle;
    baseline.xi3Mode='perPose';
    baseline.movable=baseline.movable([]);
    baseline=fr_densePilotCase(C,kin,baseline);
    assert(strcmp(baseline.passes{1}.sweptCoord,cs.passes{1}.sweptCoord));
    assert(max(abs(baseline.passes{1}.th1-cs.passes{1}.th1))<1e-10);
    d=struct('L0',baseline.L0Seed,'tendon',baseline.tendonSeed, ...
        'seriesSplit',cur.design.seriesSplit,'uv',zeros(0,3));
    % Settle actuator/tendon lengths at the fixed reference attachments.
    % Use multiple rest-length starts to give the reference a fair fit.
    opts=optimoptions('patternsearch','Display','off','MaxFunctionEvaluations',65, ...
        'InitialMeshSize',.03,'MeshTolerance',1e-4);
    lb=[max(.012,.30*d.L0),max(.004,.30*d.tendon)];
    ub=[max(d.L0*1.55,lb(1)+.02),max(d.tendon*2.5,.12)];
    best=inf;
    for k=[.75 1 1.2]
        x=min(ub,max(lb,[k*d.L0,d.tendon+(1-k)*d.L0]));
        [x,score]=patternsearch(@(x)cost(x,C,kin,baseline,d,char(b),target), ...
            x,[],[],[],[],lb,ub,[],opts);
        if score<best; best=score; selected=x; end
    end
    d.L0=selected(1); d.tendon=selected(2);
    [tau,active]=fr_primaryRoutingSweep(C,kin,baseline,d,char(b));
    T=[T; row(ref.name,b,"2019",tau,active,target)]; %#ok<AGROW>
end
temp=[tempname(out) '.mat']; save(temp,'T'); movefile(temp,file,'f');
fprintf('Era comparison checkpoint %s\n',ref.name);
end

function s=cost(x,C,kin,cs,d,branch,h)
d.L0=x(1); d.tendon=x(2);
[t,a]=fr_primaryRoutingSweep(C,kin,cs,d,branch);
deficit=max(0,abs(h)-sign(h).*t);
s=sum(deficit)/max(sum(abs(h)),eps)+.4*mean(deficit>1e-6)+2*mean(~a);
end

function R=row(group,branch,era,t,a,h)
signed=t.*sign(h);
R=table(string(group),branch,era,100*mean(a & signed>=abs(h)-1e-6), ...
    100*mean(a),100*sum(min(abs(h),max(0,signed)))/max(sum(abs(h)),eps), ...
    'VariableNames',{'group','branch','era','torque_margin_met_pct','active_rom_pct','signed_capacity_pct'});
end
