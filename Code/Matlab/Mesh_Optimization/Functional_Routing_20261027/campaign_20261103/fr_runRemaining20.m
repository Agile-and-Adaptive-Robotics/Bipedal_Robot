function fr_runRemaining20(out,budget,workers,only)
% New checkpointed 20 mm force/torque search for the 17 non-ankle groups.
if nargin<2; budget=600; end
if nargin<3; workers=10; end
C=fr_setup();
kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
groups={'med_gas','lat_gas','bifemsh','vas_int','ercspn','intobl','extobl', ...
    'psoas','iliacus','add_mag3','glut_max3','bifemlh','semimem','grac','tfl','sar','rect_fem'};
if nargin>=4; groups=intersect(groups,only,'stable'); end
branches={'rigid','compliant'};
if ~isfolder(out); mkdir(out); end
[B,G]=ndgrid(1:2,1:numel(groups));
jobs=[G(:),B(:)];
pool=[];
if workers>1
    try
        cluster=parcluster('local');
        cluster.JobStorageLocation=fullfile(out,'parallel_jobs');
        if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
        pool=parpool(cluster,min([workers,cluster.NumWorkers,size(jobs,1)]));
    catch ME
        fprintf('Worker startup failed. Running serially. %s\n',ME.message);
    end
end
if isempty(pool)
    for j=1:size(jobs,1)
        runJob(C,kin,out,groups{jobs(j,1)},branches{jobs(j,2)},budget);
    end
else
    parfor j=1:size(jobs,1)
        runJob(C,kin,out,groups{jobs(j,1)},branches{jobs(j,2)},budget);
    end
end
end

function runJob(C,kin,out,group,branch,budget)
name=[group '_' branch];
file=fullfile(out,[name '.mat']);
reportFile=fullfile(out,[name '.json']);
if isfile(file) && isfile(reportFile); return; end
try
    if isfile(file)
        S=load(file,'record'); record=S.record;
    else
        [chosen,record,trials]=fitJob(C,kin,out,group,branch,budget);
        temp=[tempname(out) '.mat'];
        save(temp,'chosen','record','trials','-v7');
        movefile(temp,file,'f');
    end
    temp=[tempname(out) '.json'];
    fid=fopen(temp,'w'); assert(fid>=0);
    fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
    movefile(temp,reportFile,'f');
    fprintf('DONE %s %d x 20 mm, G%d, target %.1f%%, no force %.1f%%, local pass %d\n', ...
        name,record.count,record.gain,100*record.fractionMet,100*record.noForce,record.localPass);
catch ME
    fid=fopen(fullfile(out,[name '_error.txt']),'w');
    fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
    fprintf('FAILED %s: %s\n',name,ME.message);
end
end

function [chosen,record,trials]=fitJob(C,kin,out,group,branch,budget)
name=[group '_' branch];
rng(sum(double(name)),'twister');
audit=fullfile(C.camp,'audit_20261105','reoptimization',group);
S=load(fullfile(audit,['coordinated_' branch '.mat']),'chosen');
warm=S.chosen;
cs=warm.caseDef;
cs.diaStr='20'; cs.diaNum=20;
opts=struct('quiet',true,'envelopes',{{}});
targets=loadTargets(out,group,'search');
counts=unique([cs.seed.nBPA,1,2,3,min(8,cs.seed.nBPA+1)],'stable');
if strcmp(group,'vas_int'); counts=unique([max(2,cs.seed.nBPA),2,3,4,6,8],'stable'); end
gains=unique([1,max(cs.txConfig.gain),2,3],'stable');
screen=struct('caseDef',{},'x',{},'score',{},'gain',{});
for n=counts
    for gain=gains
        c=cs;
        c.seed.nBPA=n;
        c.txConfig.nPulleyBPA=[1 1];
        c.txConfig.gain=[gain,1];
        reeving=fr_travelReeving(gain);
        c.txConfig.tackleLineParts=[reeving.movingFalls,1];
        c.hasTackle=gain>1;
        if c.hasTackle
            if strcmp(c.txConfig.tackleFamily,'none'); c.txConfig.tackleFamily='L'; end
        else
            c.txConfig.tackleFamily='none';
        end
        G=fr_proximalCoupled(C,kin,c,targets);
        x=min(c.ub,max(c.lb,fr_packX(c,warm.design)));
        y=min(c.ub,max(c.lb,c.xSeed));
        values=[objective(C,kin,c,x,branch,opts,G),objective(C,kin,c,y,branch,opts,G)];
        [f,k]=min(values);
        if k==2; x=y; end
        screen(end+1)=struct('caseDef',c,'x',x,'score',f,'gain',gain); %#ok<AGROW>
    end
end
[~,order]=sort([screen.score]);
pick=order(1:min(3,numel(order)));
per=max(3,floor(budget/numel(pick)));
trials=cell(size(pick));
fprintf('%s: %d sizing screens, %d continuous evaluations scheduled.\n',name,numel(screen),per*numel(pick));
for j=1:numel(pick)
    f=fullfile(out,sprintf('%s_sizing_%d.mat',name,j));
    if isfile(f)
        S=load(f,'trial'); trials{j}=S.trial; continue
    end
    s=screen(pick(j)); c=s.caseDef; x=s.x; score=s.score;
    G=fr_proximalCoupled(C,kin,c,targets);
    free=find(c.ub-c.lb>1e-10);
    fun=@(v)objective(C,kin,c,replace(x,free,v),branch,opts,G);
    so=optimoptions('surrogateopt','Display','off','UseParallel',false, ...
        'MaxFunctionEvaluations',per,'InitialPoints',x(free)');
    [v,fx,~,info]=surrogateopt(fun,c.lb(free),c.ub(free),so);
    if fx<score; x(free)=v; score=fx; end
    trial=struct('caseDef',c,'x',x,'score',score,'evaluations',info.funccount);
    temp=[tempname(out) '.mat']; save(temp,'trial','-v7'); movefile(temp,f,'f');
    trials{j}=trial;
end
[~,i]=min(cellfun(@(t)t.score,trials));
best=trials{i}; cs=best.caseDef; x=best.x;
free=find(cs.ub-cs.lb>1e-10);
G=fr_proximalCoupled(C,kin,cs,targets);
fun=@(v)objective(C,kin,cs,replace(x,free,v),branch,opts,G);
po=optimoptions('patternsearch','Display','off','UseParallel',false, ...
    'MaxFunctionEvaluations',125,'MeshTolerance',1e-4,'StepTolerance',1e-4);
[v,f]=patternsearch(fun,x(free),[],[],[],[],cs.lb(free),cs.ub(free),[],po);
if f<best.score; x(free)=v; best.score=f; end
[guide,~]=fr_screenViaGuides(C,cs,kin,x,branch,opts,fr_objective(C,cs,kin,x,branch,opts));
if numel(guide.rows)>numel(cs.rows)
    test=objective(C,kin,guide,x,branch,opts,fr_proximalCoupled(C,kin,guide,targets));
    if test<best.score; cs=guide; best.score=test; end
end
d=fr_unpackX(cs,x);
check=fr_densePilotCase(C,kin,cs);
[~,rep]=fr_objective(C,check,kin,fr_packX(check,d),branch,opts);
noForce=0; loadN=0;
for j=1:numel(check.passes)
    qg=struct('passIdx',j,'qList',{check.passes{j}.qStructs});
    m=fr_mechBuild(C,check,d,kin,branch,qg);
    noForce=max(noForce,mean(~isfinite(m.Ftendon)|m.Ftendon<=1e-6|m.PulleySlack|m.PulleyInfeasible));
    loadN=max(loadN,max([0; m.Ftendon(isfinite(m.Ftendon))]));
end
[P,B,ids]=fr_locForPass(C,cs,d,kin,1);
for j=1:numel(ids)
    T=fr_chainT(kin,kin.qDefault,B{j},cs.rowOwning{ids(j)});
    P(j,:)=(T(1:3,1:3)*P(j,:)'+T(1:3,4))';
end
gain=max(cs.txConfig.gain); topology=fr_travelReeving(gain);
record=struct('group',group,'branch',branch,'diameter_mm',20,'count',cs.seed.nBPA, ...
    'gain',gain,'sheaves',cs.seed.nBPA*topology.sheaveCount,'localPass',rep.feasible, ...
    'fractionMet',rep.fractionPosesMet,'noForce',noForce,'peakTendon_N',loadN, ...
    'constraints',rep.cinfo,'rest_mm',1000*d.L0,'tendon_mm',1000*d.tendon, ...
    'coordinate',check.passes{1}.sweptCoord,'angles_deg',rad2deg(check.passes{1}.th1), ...
    'robot_Nm',rep.tauPrimary,'human_Nm',rep.tauTarget, ...
    'evaluations',sum(cellfun(@(t)t.evaluations,trials)), ...
    'guideAdded',numel(cs.rows)>numel(warm.caseDef.rows),'points_m',P, ...
    'bodies',{cs.rowOwning(ids)},'model',C.osimRobot);
G=fr_proximalCoupled(C,kin,cs,loadTargets(out,group,'dense'));
if ~isempty(G)
    record.coupled=coupledTorque(C,kin,G,d,branch);
    record.noForce=max(record.noForce,record.coupled.noForceFraction);
    record.localPass=record.localPass && record.coupled.noForceFraction==0 && ...
        record.coupled.boneMargin_m>=-1e-4 && record.coupled.tendonMargin_m>=-1e-6;
end
chosen=struct('caseDef',cs,'design',d,'rep',rep,'xBest',x,'fBest',best.score);
end

function f=objective(C,kin,cs,x,branch,opts,G)
try
    f=fr_objective(C,cs,kin,x,branch,opts);
    if ~isempty(G)
        v=coupledTorque(C,kin,G,fr_unpackX(cs,x),branch);
        f=f+100*v.deficit+40*v.romDeficit+1000*v.noForceFraction+ ...
            1e5*max(0,-min(v.boneMargin_m,v.tendonMargin_m));
    end
    if ~isfinite(f) || ~isreal(f); f=1e12; end
catch
    f=1e12;
end
end

function r=coupledTorque(C,kin,G,d,branch)
cs=G.cs;
[m,meta]=fr_mechBuild(C,cs,d,kin,branch,G.qg);
a=fr_coordinateTorque(cs,kin,m,1,1,G.indices);
cs.passes{1}.sweptCoord=G.target.coordinate2;
b=fr_coordinateTorque(cs,kin,m,1,2,G.indices);
F=m.cell.Ftendon(:);
bad=~isfinite(a)|~isfinite(b)|~isfinite(F)|F<=1e-6| ...
    m.cell.PulleySlack(:)|m.cell.PulleyInfeasible(:)| ...
    m.cell.sContraction(:)>m.rest-m.kmax+1e-9;
a(bad)=0; b(bad)=0;
ha=G.target.target1_Nm(:); hb=G.target.target2_Nm(:);
da=max(0,abs(ha)-a.*sign(ha)); db=max(0,abs(hb)-b.*sign(hb));
meta.qg.qList=cs.passes{1}.qStructs;
m.nG=numel(meta.qg.qList);
m.strain_p=m.cell.sContraction(:)/m.rest;
m.sleeveLength_m=m.rest-m.cell.sContraction(:);
caps=fr_sleevePrimitives(C,cs,d,kin,m,meta);
ob=fr_obstacles(C,cs,meta.loc,meta.rowBodies,kin,meta.qg.qList,{},caps);
span=fr_spans(C,kin,meta.loc,meta.rowBodies,meta.qg.qList,'knee');
r=struct('axes',G.target,'robot1_Nm',a,'robot2_Nm',b, ...
    'fraction1',mean(da<=1e-6 & ~bad),'fraction2',mean(db<=1e-6 & ~bad), ...
    'deficit',.5*(sum(da)/max(sum(abs(ha)),eps)+sum(db)/max(sum(abs(hb)),eps)), ...
    'romDeficit',.5*(mean(da>1e-6)+mean(db>1e-6)), ...
    'noForceFraction',mean(bad),'boneMargin_m',ob.boneMargin, ...
    'tendonMargin_m',d.tendon-span.requiredSpan_m);
end

function T=loadTargets(out,group,label)
file=fullfile(out,'targets',[group '_' label '.json']);
T=[];
if isfile(file); T=jsondecode(fileread(file)); end
end

function x=replace(x,free,z)
x(free)=z;
end
