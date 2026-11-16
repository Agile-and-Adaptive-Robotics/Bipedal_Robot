function fr_splitVastiGrac(out,pack,refine,polishFrom)
% Three independent vasti cam routes and a two-guide gracilis bypass.
if nargin<3; refine=false; end
if nargin<4; polishFrom=''; end
if ~isfolder(out); mkdir(out); end
C=fr_setup(); kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
upper=fullfile(C.camp,'run_20261106_remaining20');
T=jsondecode(fileread(fullfile(upper,'targets','vasti_individual.json')));
P=jsondecode(fileread(fullfile(pack,'dense','poses.json')));
stress=cell(1,size(P.angles_deg,1));
for i=1:numel(stress)
    q=kin.qDefault;
    for j=1:numel(P.keys); q.(P.keys{j})=deg2rad(P.angles_deg(i,j)); end
    stress{i}=q;
end
old=read(fullfile(pack,'dense','grac_compliant.json')); old=old(end);
fit=.025;
cam=loadCam(C,kin);
groups={'vas_med_r','vas_int_r','vas_lat_r'};
branches={'rigid','compliant'};
for b=1:2
    branch=branches{b};
    packV=cell(1,3);
    for j=1:3
        file=fullfile(out,[groups{j} '_' branch '.mat']);
        if ~isempty(polishFrom)
            S=load(fullfile(polishFrom,'flat_rows_refined',[groups{j} '_' branch '.mat']),'winner');
            winner=S.winner; save(file,'winner','-v7');
        elseif isfile(file)
            S=load(file,'winner'); winner=S.winner;
        else
            winner=fitVastus(C,kin,cam,T.(groups{j}),j,branch,fit,refine);
            save(file,'winner','-v7');
        end
        packV{j}=winner;
        fprintf('%s %s: %g BPAs, %.1f%% target, %.1f%% no force, %.1f mm femoral clearance\n', ...
            groups{j},branch,winner.count,100*winner.coverage,100*winner.noForce,1000*winner.bone);
    end
    file=fullfile(out,['grac_' branch '.mat']);
    if ~isempty(polishFrom)
        dense=read(fullfile(upper,'targets','grac_dense.json'));
        a=load(fullfile(polishFrom,['grac_' branch '.mat']),'winner');
        b0=load(fullfile(polishFrom,'flat_rows_refined',['grac_' branch '.mat']),'winner');
        winner=polish(C,kin,upper,{a.winner,b0.winner},packV,dense,stress,branch,fit);
        save(file,'winner','-v7');
    elseif isfile(file)
        S=load(file,'winner'); winner=S.winner;
    else
        dense=read(fullfile(upper,'targets','grac_dense.json'));
        search=read(fullfile(upper,'targets','grac_search.json'));
        winner=fitGrac(C,kin,old.lines.points_m,packV,search,dense,stress,branch,fit,refine);
        winner=checkBiPulley(C,kin,upper,winner,packV,dense,stress,branch,fit);
        save(file,'winner','-v7');
    end
    [pairGap,insideGap]=vastiGaps(packV,fit);
    result=struct('branch',branch,'vasti',{packV},'grac',winner, ...
        'vasti_pair_gap_mm',1000*pairGap,'within_pair_gap_mm',1000*insideGap, ...
        'stressPoseCount',numel(stress),'model',C.osimRobot);
    write(fullfile(out,['split_' branch '.json']),result);
    fprintf('Gracilis %s: %s, hip %.1f%%, knee %.1f%%, no force %.1f%%, vasti gap %+.1f mm\n', ...
        branch,winner.topology,100*winner.coverage(1),100*winner.coverage(2), ...
        100*winner.noForce,1000*winner.vastiGap);
end
end

function best=polish(C,kin,upper,seeds,V,dense,stress,branch,fit)
score=inf;
q=kin.qDefault;
for key=fieldnames(q)'; q.(key{1})=0; end
home=fr_chainT(kin,q,'femur_r','pelvis');
for r=seeds
    seed=r{1};
    for change=[0,0,0;-.025,0,0;-.045,0,0;.025,0,0; ...
            0,-.020,0;-.025,-.020,0;0,0,-.012;-.025,0,-.012;0,0,.012]'
        t=seed;
        t.points_m(3,:)=t.points_m(3,:)+change';
        t.points_m(2,3)=t.points_m(2,3)+change(3);
        atHome=apply(home,t.points_m(2:3,:));
        if min(atHome(:,3))<.001 || norm(diff(t.points_m(2:3,:)))<.020; continue; end
        t=checkBiPulley(C,kin,upper,t,V,dense,stress,branch,fit,true);
        fail=[t.bone,t.vastiGap,t.viaBone,t.bay];
        f=1e4*sum(max(0,-fail))+1200*t.noForce+170*t.commonNoForce+ ...
            45*(1-t.coverage(2))+45*(1-t.signedCapacity(2))+30*(1-t.coverage(1));
        if f<score
            score=f; best=t;
            fprintf('Polish %s: knee %.1f%%, gap %+.1f mm, common no force %.1f%%\n', ...
                branch,100*t.coverage(2),1000*t.vastiGap,100*t.commonNoForce);
        end
    end
end
end

function cam=loadCam(C,kin)
path=fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Results','Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S=load(path,'ctxCamBest','predBest');
ctx=S.ctxCamBest; pred=S.predBest;
ids=find(ctx.phiD>=-120-.02 & ctx.phiD<=.02);
ang=ctx.phiD(ids); [q,I]=sort(ang(:)); ids=ids(I);
parts=cell(numel(ids),1); tib=zeros(numel(ids),3);
for k=1:numel(ids)
    i=ids(k);
    route=pred.Location(:,:,i);
    active=logical(pred.routeInfo.active(:,i)); active(1:2)=false; active(end)=true;
    r=find(active); d=r>=ctx.CrossPoint;
    nodes=route(r,:);
    nodes(d,:)=apply(ctx.T_Pam(:,:,i),nodes(d,:));
    parts{k}=nodes;
    tib(k,:)=apply(inv(fr_robotKneeT(kin.robotKnee,q(k))),nodes(end,:));
end
[~,h]=min(abs(q));
cam=struct('q',q,'parts',{parts},'entry',pred.routeInfo.seed20(2,:), ...
    'tibialEnd',tib(h,:),'home',h);
end

function best=fitVastus(C,kin,cam,h,slot,branch,fit,refine)
target=interp1(h.knee_deg,h.target_Nm,cam.q,'pchip','extrap');
% Compact anterior array. Two channels per group, all parallel to the shaft.
centers=[.076,-.025;.121,-.025;.098,.056];
shift=[0,.015,.028]; pitch=.043;
if refine
    centers=[.078,-.025;.128,-.025;.101,.068];
    shift=[0,.002]; pitch=.046;
end
base=centers(slot,:);
best=struct('score',inf);
q0=kin.qDefault; q0.knee_angle_r=0;
home=fr_chainT(kin,q0,'femur_r','pelvis');
for y0=[-.025,-.055,-.085,-.115]
    for y1=[-.330,-.370,-.405]
        for dx=shift
            x=base(1)+dx; z=base(2);
            origin=[x,y0,z]; guide=[x,y1,z];
            channel=fr_parallelRows(2,pitch);
            p=[origin;guide];
            for ap=[0,.030,.060]
                terminal=cam.tibialEnd+[ap,0,z];
                [lengths,paths,arm,bend,lead]=camLength(C,kin,cam,p,terminal);
                bay=abs(y1-y0)-2*fit;
                clearance=inf;
                for n=1:2
                    A=origin+[channel(n,2),0,channel(n,1)];
                    B=guide+[channel(n,2),0,channel(n,1)];
                    clearance=min(clearance,segmentBone(C,'femur_r',A,B)-.019-.001);
                    w=apply(home,[A;B]); 
                    if min(w(:,3))<.001; clearance=min(clearance,min(w(:,3))-.001); end
                end
                for gain=1:3
                    [d,e]=sizeLine(C,lengths,bay,arm,target,2,gain,fit,branch,'ext');
                    score=e.score+3e3*max(0,-clearance)+3e3*max(0,-lead)+ ...
                        .3*norm(origin-h.origin_m)/.025;
                    if score<best.score
                        best=struct('score',score,'points_m',p,'terminal_m',terminal, ...
                            'channel_m',channel,'count',2,'gain',gain,'rest_m',d(1), ...
                            'tendon_m',d(2),'radius_m',e.radius,'sleeve_m',e.sleeve, ...
                            'knee_deg',cam.q,'robot_Nm',e.tau,'human_Nm',target,'tension_N',e.tension, ...
                            'coverage',e.coverage,'noForce',mean(e.bad),'bay',e.bay, ...
                            'bone',clearance,'cableBone',lead,'camBendRadius_m',bend, ...
                            'home_cam_femur_m',paths{cam.home},'contactPath_m',{paths}, ...
                            'origin_movement_mm',1000*norm(origin-h.origin_m));
                    end
                end
            end
        end
    end
end
end

function [L,P,arm,bend,clearance]=camLength(C,kin,cam,p,terminal)
L=zeros(size(cam.q)); P=cell(size(L)); bend=inf; clearance=inf;
for i=1:numel(L)
    T=fr_robotKneeT(kin.robotKnee,cam.q(i));
    path=cam.parts{i}(1:end-1,:)+[0,0,p(2,3)];
    dest=apply(T,terminal);
    path=[p(2,:);path;dest];
    path=path([true;vecnorm(diff(path),2,2)>1e-6],:);
    P{i}=path;
    L(i)=norm(diff(p))+sum(vecnorm(diff(path),2,2));
    % Radii of the recorded contact polyline. Keep free guide corners separate.
    for k=3:size(path,1)-2
        a=path(k,:)-path(k-1,:); b=path(k+1,:)-path(k,:);
        bend=min(bend,norm(a)*norm(b)*norm(a+b)/(2*max(norm(cross(a,b)),eps)));
    end
    if mod(i,4)==1 && size(path,1)>1
        a=path(end-1,:); b=path(end,:);
        a=a+.04*(b-a); b=b-.08*(b-a);
        clearance=min(clearance,min(segmentBone(C,'femur_r',a,b), ...
            segmentBone(C,'tibia_r',apply(inv(T),a),apply(inv(T),b)))-.002);
    end
end
arm=-gradient(L)./gradient(deg2rad(cam.q));
end

function best=fitGrac(C,kin,old,vasti,search,dense,stress,branch,fit,refine)
coarse=gracGrid(kin,search);
fine=gracGrid(kin,dense);
common=prep(kin,stress);
best=struct('score',inf);
candidates=struct('score',{},'p',{},'gain',{},'d',{});
q0=kin.qDefault;
for key=fieldnames(q0)'; q0.(key{1})=0; end
home=fr_chainT(kin,q0,'femur_r','pelvis');
leads=[-.030,0,.030]; stopY=-.418; jogZ=0;
if refine
    leads=[-.085,-.060,-.035,0,.025];
    stopY=[-.405,-.448]; jogZ=[-.012,0,.012];
end
for x=[-.095,-.070,-.045,-.020,.005]
    for z=[-.078,-.060,-.042,-.024]
        for y=[-.285,-.325,-.365,-.400]
            for dx=leads
                for y2=stopY
                for dz=jogZ
                p=[old(1,:);x,y,z;x+dx,y2,z+dz;old(2:3,:)];
                if norm(p(3,:)-p(2,:))<.024; continue; end
                atHome=apply(home,p(2:3,:));
                if min(atHome(:,3))<.001; continue; end
                guideClear=min(fr_surfaceDistance(C,'femur_r',p(2:3,:)))-.008;
                if guideClear<0; continue; end
                [L,bay,arm]=gracPath(coarse,p,fit);
                for gain=1:3
                    [d,e]=sizeLine(C,L,bay,arm,coarse.target,1,gain,fit,branch,'flx');
                    score=e.score+0.025*norm(p(2,:)-p(3,:))*1000;
                    candidates(end+1)=struct('score',score,'p',p,'gain',gain,'d',d); %#ok<AGROW>
                end
                end
                end
            end
        end
    end
end
[~,I]=sort([candidates.score]);
I=unique([I(1:min(85,numel(I))), diverse(candidates)],'stable');
for j=I
    c=candidates(j);
    [L,bay,arm]=gracPath(fine,c.p,fit);
    [d,e]=sizeLine(C,L,bay,arm,fine.target,1,c.gain,fit,branch,'flx');
    [S,bayS]=gracPath(common,c.p,fit);
    M=closure(C,S,bayS,zeros(numel(S),1),zeros(numel(S),1),1,c.gain,d,fit,branch,'flx');
    [bone,gap,via,hit]=gracClearance(C,common,c.p,M,vasti,fit);
    score=e.score+100*mean(M.bad)+4e3*sum(max(0,-[bone,gap,via,M.bay]));
    if score<best.score
        best=struct('score',score,'topology','two femoral cable guides','points_m',c.p, ...
            'bodies',{{'pelvis','femur_r','femur_r','tibia_r','tibia_r'}}, ...
            'count',1,'gain',c.gain,'rest_m',d(1),'tendon_m',d(2),'radius_m',M.radius, ...
            'sleeve_m',M.sleeve,'coverage',e.coverage,'noForce',mean(e.bad), ...
            'commonNoForce',mean(M.bad),'bay',min(e.bay,M.bay),'bone',bone,'vastiGap',gap, ...
            'viaBone',via,'governingPose',hit,'robot_Nm',e.tau,'human_Nm',fine.target, ...
            'signedCapacity',e.capacity,'axes',dense);
    end
end
assert(isfield(best,'topology'),'No gracilis guide candidate');
end

function I=diverse(c)
I=[];
for x=unique(arrayfun(@(v)v.p(2,1),c))
    for z=unique(arrayfun(@(v)v.p(2,3),c))
        ii=find(arrayfun(@(v)v.p(2,1)==x && v.p(2,3)==z,c));
        if isempty(ii); continue; end
        [~,k]=min([c(ii).score]); I(end+1)=ii(k); %#ok<AGROW>
    end
end
end

function G=gracGrid(kin,T)
[a,b]=ndgrid(T.axis1_deg,T.axis2_deg);
q=cell(1,numel(a));
for i=1:numel(a)
    p=kin.qDefault; p.hip_adduction_r=deg2rad(a(i)); p.knee_angle_r=deg2rad(b(i)); q{i}=p;
end
G=prep(kin,q); G.shape=size(a); G.a=deg2rad(T.axis1_deg); G.b=deg2rad(T.axis2_deg);
G.target=[T.target1_Nm(:),T.target2_Nm(:)];
end

function G=prep(kin,q)
N=numel(q); G.P=zeros(4,4,N); G.F=G.P; G.T=G.P; G.knee=zeros(N,1);
for i=1:N
    G.P(:,:,i)=fr_bodyT(kin,q{i},'pelvis'); G.F(:,:,i)=fr_bodyT(kin,q{i},'femur_r');
    G.T(:,:,i)=fr_bodyT(kin,q{i},'tibia_r'); G.knee(i)=rad2deg(q{i}.knee_angle_r);
end
end

function [L,bay,arm]=gracPath(G,p,fit)
O=transform(G.P,p(1,:)); A=transform(G.F,p(2,:)); B=transform(G.F,p(3,:));
C=transform(G.T,p(4,:)); D=transform(G.T,p(5,:));
l=vecnorm(O-A,2,2); bay=l-2*fit;
L=l+vecnorm(A-B,2,2)+vecnorm(B-C,2,2)+vecnorm(C-D,2,2);
arm=[];
if isfield(G,'shape')
    z=reshape(L,G.shape);
    [b,a]=gradient(z,G.b,G.a);
    arm=-[a(:),b(:)];
end
end

function [design,best]=sizeLine(C,L,bay,arm,h,n,g,fit,branch,type)
best=struct('score',inf);
lo=max(.105,range(L)/(.24*g)*.98); hi=max(lo*1.03,min(bay)*1.03);
for rest=unique([linspace(lo,hi,4),max(.105,min(bay)+.004)])
    tendon=max(.015,max(L)-rest-2*fit);
    for shift=[-.010,-.005,.001,.007]
        d=[rest,max(.015,tendon+shift)];
        e=closure(C,L,bay,arm,h,n,g,d,fit,branch,type);
        if e.score<best.score; design=d; best=e; end
    end
end
end

function E=closure(C,L,bay,arm,h,n,g,d,fit,branch,type)
rest=d(1); cap=.24*rest; F=maxBPAforce(rest,'20',C.pressure)/g;
xi=0; compliance=0;
if strcmp(branch,'compliant')
    xi=C.Xi.(['Xi0_' type]);
    % One load path, two directional brackets and one tendon.
    compliance=2/(.75*C.Xi.Xi1+.25*C.Xi.Xi2)+d(2)/1.8e4;
end
demand=rest+d(2)+2*fit+xi-L;
% Fast monotone closure screen, using the existing Festo force law.
grid=linspace(0,1,41)';
force=F*festo4(20,grid,C.pressure);
residual=g*cap*grid-compliance*force;
s=interp1(residual,grid,demand,'linear','extrap');
bad=s<0|s>1|~isfinite(s);
s=max(0,min(1,s));
t=F*festo4(20,s,C.pressure); t(bad)=0;
torque=n*t.*arm;
deficit=max(0,abs(h)-torque.*sign(h));
capacity=1-sum(deficit,1)./max(sum(abs(h),1),eps);
axial=min(bay-(rest-cap*s));
E=struct('score',mean(sum(deficit,1)./max(sum(abs(h),1),eps))*100+ ...
    35*mean(deficit>1e-6,'all')+800*mean(bad)+3500*max(0,-axial), ...
    'tau',torque,'coverage',mean(deficit<=1e-6 & ~bad,1),'bad',bad, ...
    'radius',fr_diameter(C,'20',.24*s)/2000,'sleeve',rest-cap*s, ...
    'tension',t,'bay',axial,'capacity',capacity);
end

function [bone,gap,via,witness]=gracClearance(C,G,p,E,vasti,fit)
O=transform(G.P,p(1,:)); A=transform(G.F,p(2,:)); B=transform(G.F,p(3,:));
D=transform(G.T,p(4,:)); last=transform(G.T,p(5,:));
bone=inf; gap=inf; via=inf; witness=0;
for i=1:size(O,1)
    u=A(i,:)-O(i,:); u=u/norm(u);
    start=O(i,:)+fit*u; finish=start+E.sleeve(i)*u;
    F=G.F(:,:,i); invF=inv(F); invT=inv(G.T(:,:,i));
    prox=apply(invF,[start;finish]);
    bone=min(bone,segmentBone(C,'femur_r',prox(1,:),prox(2,:))-E.radius(i)-.001);
    wire=[finish;A(i,:);B(i,:);D(i,:);last(i,:)];
    via=min(via,min(segmentBone(C,'femur_r',p(2,:),p(3,:))-.002, ...
        segmentBone(C,'tibia_r',apply(invT,B(i,:)),D(i,:)*invT(1:3,1:3)'+invT(1:3,4)')-.002));
    for v=1:3
        solid=vastiParts(vasti{v},G.knee(i),fit);
        for j=1:size(solid,1)
            t=solid(j,:); a=apply(F,t(1:3)); b=apply(F,t(4:6));
            g=fr_segmentDistance(start-fit*u,finish+fit*u,a,b)-E.radius(i)-t(7);
            g=min(g,fr_segmentDistance(A(i,:),A(i,:),a,b)-.008-t(7));
            g=min(g,fr_segmentDistance(B(i,:),B(i,:),a,b)-.008-t(7));
            for k=1:size(wire,1)-1
                g=min(g,fr_segmentDistance(wire(k,:),wire(k+1,:),a,b)-.001-t(7));
            end
            if g<gap; gap=g; witness=i; end
        end
    end
end
end

function s=vastiParts(v,q,fit)
[~,i]=min(abs(v.knee_deg-q));
p=v.points_m; u=diff(p); u=u/norm(u);
s=zeros(v.count,7);
for j=1:v.count
    d=[v.channel_m(j,2),0,v.channel_m(j,1)];
    a=p(1,:)+d; b=a+(v.sleeve_m(i)+2*fit)*u;
    s(j,:)=[a,b,v.radius_m(i)];
end
end

function [gap,inside]=vastiGaps(V,fit)
gap=inf; inside=inf;
for q=linspace(-120,0,61)
    P=cellfun(@(v)vastiParts(v,q,fit),V,'UniformOutput',false);
    for j=1:3
        A=P{j}; inside=min(inside,partGap(A(1,:),A(2,:)));
        for k=j+1:3
            B=P{k};
            for u=1:2
                for v=1:2; gap=min(gap,partGap(A(u,:),B(v,:))); end
            end
        end
    end
end
end

function r=checkBiPulley(C,kin,upper,r,V,target,stress,branch,fit,fixedRest)
if nargin<10; fixedRest=false; end
file=fullfile(upper,'paired_refit',['grac_' branch '.mat']);
S=load(file,'chosen');
cs=S.chosen.caseDef;
cs.rowOwning=r.bodies;
cs.seed.nBPA=1;
cs.hasTackle=r.gain>1;
cs.fittingLength=fit;
cs.txConfig.gain=[r.gain,1];
cs.txConfig.tackleLineParts=[r.gain,1];
if isfield(cs.txConfig,'pulleyExitIndex'); cs.txConfig=rmfield(cs.txConfig,'pulleyExitIndex'); end
G=fr_proximalCoupled(C,kin,cs,target);
G.qg.resolvedRoute=struct('points',r.points_m,'frames',{r.bodies});
fine=gracGrid(kin,target);
[~,bay,arm]=gracPath(fine,r.points_m,fit);
common=prep(kin,stress);
[~,bayS]=gracPath(common,r.points_m,fit);
csS=G.cs;
csS.passes{1}.qStructs=stress;
qg=rmfield(G.qg,'qList2'); qg.qList=stress; qg.pairedPoseGrid=true;
best=inf;
rests=[-.012,0,.012]; tendons=[-.018,-.009,0,.009,.018];
if fixedRest; rests=0; tendons=[-.009,0,.009]; end
for dL=rests
    for dT=tendons
        d=struct('L0',r.rest_m+dL,'tendon',r.tendon_m+dT,'seriesSplit',1,'uv',zeros(0,3));
        if min(d.L0,d.tendon)<=0; continue; end
        m=fr_mechBuild(C,G.cs,d,kin,branch,G.qg);
        a=active(m,d);
        force=m.cell.Ftendon(:); force(~a)=0;
        torque=force.*arm;
        deficit=max(0,abs(fine.target)-torque.*sign(fine.target));
        s=fr_mechBuild(C,csS,d,kin,branch,qg);
        aS=active(s,d);
        L=d.L0-s.cell.sContraction(:);
        gapBay=min([bay-(d.L0-m.cell.sContraction(:));bayS-L]);
        val=800*mean(~a)+180*mean(~aS)+3500*max(0,-gapBay)+ ...
            100*mean(sum(deficit)./sum(abs(fine.target)))+35*mean(deficit>1e-6,'all');
        if val<best
            best=val;
            selected=struct('d',d,'coverage',mean(deficit<=1e-6 & a,1),'noForce',mean(~a), ...
                'commonNoForce',mean(~aS),'robot_Nm',torque,'bay',gapBay,'sleeve',L, ...
                'radius',fr_diameter(C,'20',s.cell.sContraction(:)/d.L0)/2000);
        end
    end
end
[bone,gap,via,hit]=gracClearance(C,common,r.points_m,selected,V,fit);
r.reducedCoverage=r.coverage;
r.rest_m=selected.d.L0; r.tendon_m=selected.d.tendon;
for k={'coverage','noForce','commonNoForce','robot_Nm','bay'}
    r.(k{1})=selected.(k{1});
end
r.sleeve_m=selected.sleeve; r.radius_m=selected.radius;
r.bone=bone; r.vastiGap=gap; r.viaBone=via; r.governingPose=hit;
r.signedCapacity=sum(min(abs(r.human_Nm),max(0,r.robot_Nm.*sign(r.human_Nm)))) ./ ...
    sum(abs(r.human_Nm));
r.equilibrium='BiPulley with signed torque from finite differences of the full guided path';
end

function a=active(m,d)
F=m.cell.Ftendon(:);
a=isfinite(F)&F>1e-6&~m.cell.PulleySlack(:)&~m.cell.PulleyInfeasible(:)& ...
    m.cell.sContraction(:)<=.24*d.L0+1e-9;
end

function d=partGap(a,b)
d=fr_segmentDistance(a(1:3),a(4:6),b(1:3),b(4:6))-a(7)-b(7);
end

function d=segmentBone(C,body,a,b)
n=max(2,ceil(norm(b-a)/.006)+1);
p=a+linspace(0,1,n)'.*(b-a);
d=min(fr_surfaceDistance(C,body,p))-norm(b-a)/(2*(n-1));
end

function p=transform(T,x)
n=size(T,3); p=zeros(n,3);
for i=1:n; p(i,:)=apply(T(:,:,i),x); end
end

function p=apply(T,p)
p=p*T(1:3,1:3)'+T(1:3,4)';
end

function x=read(p)
x=jsondecode(fileread(p));
end

function write(p,x)
fid=fopen(p,'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(x)); fclose(fid);
end
