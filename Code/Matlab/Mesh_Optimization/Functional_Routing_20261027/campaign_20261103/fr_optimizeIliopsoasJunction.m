function fr_optimizeIliopsoasJunction(out)
% Join the distal load transfer on the femur, retaining independent tendons.
C=fr_setup();
kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
for br={'rigid','compliant'}
    branch=br{1};
    file=fullfile(out,['iliopsoas_junction_' branch '.mat']);
    if isfile(file); continue; end
    a=load(fullfile(out,['psoas_' branch '.mat']),'chosen');
    b=load(fullfile(out,['iliacus_' branch '.mat']),'chosen');
    a=a.chosen; b=b.chosen;
    endA=nativeEnd(C,kin,a); endB=nativeEnd(C,kin,b);
    patch=fr_patch(C,'Femur','femur_r',mean([endA;endB]),25);
    xa=lengthBounds(a); xb=lengthBounds(b);
    low=[-1 -1 C.standoff_mm(1),xa(1,:),xb(1,:)];
    high=[1 1 C.standoff_mm(2),xa(2,:),xb(2,:)];
    x0=[0 0 12,a.design.L0,a.design.tendon,b.design.L0,b.design.tendon];
    x0=max(low,min(high,x0));
    rng(sum(double(branch))+170,'twister');
    so=optimoptions('surrogateopt','Display','off','UseParallel',false, ...
        'MaxFunctionEvaluations',140,'InitialPoints',x0);
    [x,score]=surrogateopt(@(x)trial(C,kin,a,b,patch,x,branch),low,high,so);
    [~,A,B,junction]=trial(C,kin,a,b,patch,x,branch);
    rows={A,B};
    records=cell(1,2);
    for k=1:2
        c=rows{k};
        d=fr_densePilotCase(C,kin,c.caseDef);
        [~,r]=fr_objective(C,d,kin,fr_packX(d,c.design),branch,struct('quiet',true,'envelopes',{{}}));
        m=fr_mechBuild(C,d,c.design,kin,branch,struct('passIdx',1,'qList',{d.passes{1}.qStructs}));
        records{k}=struct('group',c.caseDef.name,'fractionMet',r.fractionPosesMet, ...
            'localPass',r.feasible,'constraints',r.cinfo, ...
            'robot_Nm',r.tauPrimary,'human_Nm',r.tauTarget, ...
            'force_N',m.Ftendon,'angles_deg',rad2deg(d.passes{1}.th1), ...
            'noForceFraction',mean(~isfinite(m.Ftendon)|m.Ftendon<=1e-6|m.PulleyInfeasible|m.PulleySlack));
    end
    assert(max(abs(records{1}.angles_deg-records{2}.angles_deg))<1e-8);
    total=records{1}.robot_Nm+records{2}.robot_Nm;
    target=records{1}.human_Nm+records{2}.human_Nm;
    pass=all(cellfun(@(r)r.localPass,records));
    result=struct('branch',branch,'junction_femur_m',junction,'individual',{records}, ...
        'combinedRobot_Nm',total,'combinedHuman_Nm',target,'score',score, ...
        'combinedFractionMet',mean(total.*sign(target)>=abs(target)-1e-6), ...
        'localPass',pass);
    save(file,'result','rows','a','b','-v7');
    fid=fopen(fullfile(out,['iliopsoas_junction_' branch '_result.json']),'w');
    fprintf(fid,'%s\n',jsonencode(result)); fclose(fid);
    fprintf('Iliopsoas %s common femoral load transfer: combined target %.1f%%, local checks %d.\n', ...
        branch,100*result.combinedFractionMet,pass);
end
end

function [f,A,B,junction]=trial(C,kin,a,b,patch,x,branch)
junction=patch.project(x(1),x(2),x(3));
A=atJunction(a,junction,x(4:5));
B=atJunction(b,junction,x(6:7));
opts=struct('quiet',true,'envelopes',{{}});
try
    [fa,ra]=fr_objective(C,A.caseDef,kin,fr_packX(A.caseDef,A.design),branch,opts);
    [fb,rb]=fr_objective(C,B.caseDef,kin,fr_packX(B.caseDef,B.design),branch,opts);
    assert(numel(ra.tauPrimary)==numel(rb.tauPrimary));
    h=ra.tauTarget+rb.tauTarget;
    t=ra.tauPrimary+rb.tauPrimary;
    f=fa+fb+100*sum(max(0,abs(h)-t.*sign(h)))/max(sum(abs(h)),eps);
    if ~isfinite(f) || ~isreal(f); f=1e12; end
catch
    f=1e12;
end
end

function r=atJunction(r,p,lengths)
cs=r.caseDef;
row=size(cs.seedLoc,1);
assert(strcmp(cs.rowOwning{row},'femur_r'));
cs.seedLoc(row,:)=p;
k=find([cs.movable.row]==row,1);
if ~isempty(k); cs.movable(k).patch=[]; end
r.caseDef=cs;
r.design.L0=lengths(1);
r.design.tendon=lengths(2);
end

function p=nativeEnd(C,kin,r)
[P,B]=fr_locForPass(C,r.caseDef,r.design,kin,1);
T=fr_chainT(kin,kin.qDefault,B{end},'femur_r');
p=(T(1:3,1:3)*P(end,:)'+T(1:3,4))';
end

function x=lengthBounds(r)
i=3*numel(r.caseDef.movable)+(1:2);
x=[r.caseDef.lb(i)';r.caseDef.ub(i)']/10;
end
