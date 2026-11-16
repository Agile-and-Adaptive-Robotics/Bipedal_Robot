function fr_screenGastrocPinnedEnds(out)
% Knee Bowden with separate active bays on the femur and the shank.
C=fr_setup(); kin=fr_robotKneeContext(C,fr_kinematics(C.osimStock));
for group={'med_gas','lat_gas'}
    g=group{1};
    source=fullfile(C.camp,'audit_20261105','reoptimization',g,'coordinated_compliant.mat');
    S=load(source,'chosen');
    seed=S.chosen.caseDef;
    origin=seed.seedLoc(1,:);
    foot=seed.seedLoc(end,:);
    target=jsondecode(fileread(fullfile(out,'targets',[g '_search.json'])));
    dense=jsondecode(fileread(fullfile(out,'targets',[g '_dense.json'])));
    medial=strcmp(g,'med_gas');
    z=.026*(1-2*medial);
    fit=seed.fittingLength;
    cs=seed;
    cs.diaNum=20; cs.diaStr='20';
    cs.rowOwning={'femur_r','femur_r','tibia_r','tibia_r','calcn_r'};
    cs.txConfig.routingMode={'bowden','moving_via'};
    cs.txConfig.nPulleyBPA=[1 1];
    cs.seed.nBPA=1;
    cs.fittingLength=2*fit;
    for br={'rigid','compliant'}
        branch=br{1};
        file=fullfile(out,[g '_pinned_ends_' branch '.mat']);
        if isfile(file); continue; end
        best=struct('score',inf);
        for m=[.055,.050;.075,.060;.095,.067;.075,.083;.100,.080]'
            p=[origin(1)-m(2),origin(2)+m(1),origin(3); ...
                -.065,-.390,z; ...
                -.065,-.067,z; ...
                -.060,-.325,z; ...
                foot];
            bay=vecnorm(p([2,4],:)-p([1,3],:),2,2)-2*fit;
            if any(bay<.025); continue; end
            G=fr_proximalCoupled(C,kin,cs,target);
            G.qg.resolvedRoute=struct('points',p,'frames',{cs.rowOwning});
            K=geometry(C,kin,G,p,bay,fit);
            for gain=1:3
                G.cs.txConfig.gain=[gain 1];
                G.cs.txConfig.tackleLineParts=[gain 1];
                % Existence of an axial sizing interval precedes force sizing.
                lower=max(.060,range(K.length)/(gain*.24));
                upper=sum(bay)*1.08;
                if lower>=upper; continue; end
                closes=false;
                for L=linspace(lower,upper,6)
                    base=max(K.tendonMin,max(K.length)-4*fit-L);
                    for shift=[-.015,-.005,0,.006,.015,.028]
                        tendon=max(K.tendonMin,base+shift);
                        d=struct('L0',L,'tendon',tendon,'uv',zeros(0,3), ...
                            'seriesSplit',bay(1)/sum(bay));
                        e=run(C,kin,G,K,d,p,fit,branch);
                        closes=closes || ~any(e.bad);
                        for n=1:4
                            [score,~]=scoreTrial(e,target,n);
                            if score<best.score
                                best=struct('score',score,'p',p,'bay',bay,'d',d,'gain',gain,'count',n);
                            end
                        end
                    end
                end
                if closes; break; end
            end
        end
        if ~isfield(best,'p')
            result=struct('group',g,'branch',branch,'failure','No axial length window');
        else
            G=fr_proximalCoupled(C,kin,cs,dense);
            G.cs.txConfig.gain=[best.gain 1];
            G.cs.txConfig.tackleLineParts=[best.gain 1];
            G.qg.resolvedRoute=struct('points',best.p,'frames',{cs.rowOwning});
            K=geometry(C,kin,G,best.p,best.bay,fit);
            e=run(C,kin,G,K,best.d,best.p,fit,branch);
            [~,hit]=scoreTrial(e,dense,best.count);
            result=struct('group',g,'branch',branch,'count',best.count,'gain',best.gain, ...
                'rest_mm',1000*best.d.L0,'sectionRest_mm',1000*best.d.L0* ...
                    [best.d.seriesSplit,1-best.d.seriesSplit], ...
                'tendon_mm',1000*best.d.tendon,'noForceFraction',mean(e.bad), ...
                'kneeMet',hit(1),'ankleMet',hit(2),'boneMargin_mm',1000*e.boneMargin, ...
                'bayMargin_mm',1000*e.bayMargin,'tendonMargin_mm',1000*e.tendonMargin, ...
                'points_m',best.p,'bodies',{cs.rowOwning},'axes',dense, ...
                'knee_Nm',best.count*e.tau(:,1),'ankle_Nm',best.count*e.tau(:,2));
        end
        save(file,'best','result','-v7');
        fid=fopen(fullfile(out,[g '_pinned_ends_' branch '_result.json']),'w');
        fprintf(fid,'%s\n',jsonencode(result)); fclose(fid);
        if isfield(result,'noForceFraction')
            fprintf('%s %s pinned ends: no force %.1f%%, knee %.1f%%, ankle %.1f%%, bone %+.1f mm.\n', ...
                g,branch,100*result.noForceFraction,100*result.kneeMet,100*result.ankleMet,result.boneMargin_mm);
        end
    end
end
end

function K=geometry(C,kin,G,p,bay,fit)
poses=G.cs.passes{1}.qStructs;
K.length=zeros(numel(poses),1); K.arm=zeros(numel(poses),2);
tendon=K.length;
for i=1:numel(poses)
    q=poses{i};
    [K.length(i),tendon(i)]=lengthAt(kin,q,G.cs.rowOwning,p);
    for j=1:2
        coord=G.target.(sprintf('coordinate%d',j)); h=1e-5;
        a=q; b=q; a.(coord)=q.(coord)+h; b.(coord)=q.(coord)-h;
        K.arm(i,j)=-(lengthAt(kin,a,G.cs.rowOwning,p)-lengthAt(kin,b,G.cs.rowOwning,p))/(2*h);
    end
end
K.tendonMin=max(tendon)+2*(2*C.bowden.ferrule_mm+C.bowden.stopAllowance_mm)/1000;
K.bay=bay;
K.boneDistance=inf;
for j=1:2
    r=2*j-1;
    u=p(r+1,:)-p(r,:);
    n=norm(u);
    unit=u/n;
    points=p(r,:)+linspace(fit,n-fit,45)'.*unit;
    K.boneDistance=min(K.boneDistance,min(fr_surfaceDistance(C,G.cs.rowOwning{r},points)));
end
end

function [L,T]=lengthAt(kin,q,b,p)
for j=1:size(p,1)
    X=fr_bodyT(kin,q,b{j});
    p(j,:)=(X(1:3,1:3)*p(j,:)'+X(1:3,4))';
end
s=vecnorm(diff(p),2,2);
L=sum(s); T=sum(s([2 4]));
end

function e=run(C,kin,G,K,d,p,fit,branch)
G.qg.resolvedRoute=struct('points',p,'frames',{G.cs.rowOwning});
m=fr_mechBuild(C,G.cs,d,kin,branch,G.qg);
force=m.cell.Ftendon(:);
s=m.cell.sContraction(:);
bad=~isfinite(force)|force<=1e-6|m.cell.PulleySlack(:)|m.cell.PulleyInfeasible(:)| ...
    s/d.L0>.24+1e-9;
force(bad)=0;
% A series line cannot exceed the weaker section's force capacity.
split=[d.seriesSplit,1-d.seriesSplit];
capacity=maxBPAforce(d.L0*min(split),'20',C.pressure)/maxBPAforce(d.L0,'20',C.pressure);
force=force*min(1,capacity);
tau=K.arm.*force;
radius=max(fr_diameter(C,'20',s/d.L0))/2000;
bay=min(K.bay-max(d.L0-s)*split');
e=struct('tau',tau,'bad',bad,'boneMargin',K.boneDistance-radius-C.clearanceBone_m, ...
    'bayMargin',bay,'tendonMargin',d.tendon-K.tendonMin);
end

function [f,hit]=scoreTrial(e,T,n)
h=[T.target1_Nm(:),T.target2_Nm(:)];
d=max(0,abs(h)-n*e.tau.*sign(h));
hit=mean(d<=1e-6 & ~e.bad,1);
deficit=mean(sum(d)./max(sum(abs(h)),eps));
f=2e3*mean(e.bad)+1e4*sum(max(0,-[e.boneMargin,e.bayMargin,e.tendonMargin]))+ ...
    100*deficit+40*mean(1-hit)+.05*n;
end
