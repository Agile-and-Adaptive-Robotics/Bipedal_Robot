function fr_tryParallelAnkleExtensors(outputDir,exactLengths)
% Screen two direct-pull BPAs at the saved compliant extensor mount locations.
if nargin<2; exactLengths=false; end
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
if ~isfolder(outputDir); mkdir(outputDir); end
stage=fullfile(C.camp,'run_20261105_distal_stages');
names={'ext_dig_compliant','ext_hal_compliant'};
folders={'ankle_guide_arrow','ankle_mounts_600'};
branch='compliant';
n=4;
base=cell(1,2); fitted=base; coarse=base;
for g=1:2
    S=load(fullfile(stage,folders{g},[names{g} '.mat']),'result');
    base{g}=S.result.chosen;
    [e,d]=fr_evalCoupledAnkle(C,kin,base{g}.caseDef,base{g}.design,branch,n);
    coarse{g}=struct('e',e,'d',d);
    c=base{g};
    assert(c.caseDef.seed.nBPA==1 && all(c.caseDef.txConfig.gain==1));
    c.caseDef.seed.nBPA=2;
    c.caseDef.txConfig.nPulleyBPA=[2 1];
    lengthFile=fullfile(outputDir,[names{g} '_lengths.mat']);
    if isequal(exactLengths,'same')
        % Compare a literal duplicate before allowing either length to change.
    elseif isfile(lengthFile)
        S=load(lengthFile,'design'); c.design=S.design;
    else
        G=fr_ankleCoupledSetup(C,kin,c.caseDef,21);
        estimate=fr_refitMountLengths(C,kin,c.caseDef,c.design,branch,G);
        if exactLengths
            [design,lengthTrials]=solveLengths(C,kin,c,estimate,G,branch);
            save(lengthFile,'design','lengthTrials','-v7');
            c.design=design;
        else
            c.design=estimate;
        end
    end
    fitted{g}=c;
end

[angle,spread]=ndgrid(0:30:150,[1.8,2.3,2.8,3.2]);
layout=[angle(:),spread(:)];
screen=fullfile(outputDir,'layout_screen.mat');
if isfile(screen)
    S=load(screen,'trials'); trials=S.trials;
else
    pool=gcp('nocreate');
    if isempty(pool)
        try
            cluster=parcluster('local');
            cluster.JobStorageLocation=fullfile(outputDir,'parallel_jobs');
            if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
            pool=parpool(cluster,min([C.parallelMaxWorkers,cluster.NumWorkers,6]));
        catch ME
            fprintf('Worker startup failed. Screening serially. %s\n',ME.message);
        end
    end
    jobs=cell(size(layout,1),2);
    if isempty(pool)
        for j=1:numel(jobs)
            jobs{j}=evaluateLayout(C,kin,fitted,coarse,layout,n,j);
        end
    else
        parfor j=1:numel(jobs)
            jobs{j}=evaluateLayout(C,kin,fitted,coarse,layout,n,j);
        end
    end
    trials=jobs;
    save(screen,'trials','layout','-v7');
end

% Test the strongest individual layouts together on shared poses.
rank=cell(1,2);
for g=1:2
    s=cellfun(@(t) penalty(t.e,t.layer,t.c),trials(:,g));
    [~,ix]=sort(s);
    rank{g}=ix(1:min(6,numel(ix)));
end
best=inf; pair=[];
for i=rank{1}'
    for j=rank{2}'
        a=trials{i,1}; b=trials{j,2};
        L=compare(C,kin,a,b);
        score=penalty(a.e,L,a.c)+penalty(b.e,L,b.c);
        if score<best
            best=score; pair=[i,j];
        end
    end
end
final=cell(1,2); retained=final;
for g=1:2
    c=trials{pair(g),g}.c;
    [e,d]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,branch,21);
    final{g}=struct('c',c,'e',e,'d',d);
    c0=base{g};
    [e0,d0]=fr_evalCoupledAnkle(C,kin,c0.caseDef,c0.design,branch,21);
    retained{g}=struct('c',c0,'e',e0,'d',d0);
end
layer=compare(C,kin,final{1},final{2});
independent={compare(C,kin,final{1},retained{2}),compare(C,kin,retained{1},final{2})};
save(fullfile(outputDir,'parallel_extensors_dense.mat'),'final','retained','layer','independent','pair','-v7');
for g=1:2
    name=names{g}; t=final{g}; c=t.c;
    record=jsondecode(fileread(fullfile(stage,folders{g},[name '.json'])));
    drop=intersect(fieldnames(record),{'surrogateEvaluations','initialScore','finalScore', ...
        'guideShift_mm','previous','baselineFolder','innerRoute','layer'});
    record=rmfield(record,drop);
    record.count=2;
    record.rest_mm=1000*c.design.L0;
    record.tendon_mm=1000*c.design.tendon;
    record.bundleSpread=c.caseDef.bundleSpread;
    record.bundleClock_deg=c.caseDef.bundleClock_deg;
    record.coupled=t.e;
    record.localConstraintsPass=t.e.localPass;
    record.foundLengthFit=~any(t.e.noForce);
    record.legacySweep=false;
    record.baselineCoupled=retained{g}.e;
    record.sameGroupLocalPass=t.e.localPass;
    record.pairedLocalPass=t.e.localPass && final{3-g}.e.localPass && layerPass(layer);
    record.coupled.localPass=record.pairedLocalPass;
    record.layerWithSingleOther=independent{g};
    if g==1; record.layer=layer; end
    if g==1
        record.innerRoute=jsondecode(fileread(fullfile(stage,folders{2},[names{2} '.json'])));
        record.innerRoute=struct('points_m',record.innerRoute.points_m,'bodies',{record.innerRoute.bodies});
    end
    fid=fopen(fullfile(outputDir,[name '.json']),'w'); assert(fid>=0);
    fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
    fprintf('%s: %d x %g mm, spread %.1f, orientation %.0f deg, bone %+.2f mm, sleeve gap %+.2f mm, fitting gap %+.2f mm, no force %.1f%%\n', ...
        name,record.count,record.diameter_mm,record.bundleSpread,record.bundleClock_deg, ...
        1000*t.e.boneMargin_m,1000*t.e.bundleMargin_m,1000*t.e.bundleHardwareMargin_m,100*mean(t.e.noForce));
end
fprintf('Two bundles: ordering %+.2f mm, fitting margin %+.2f mm, soft sleeve overlap %.2f mm\n', ...
    1000*layer.orderingMargin_m,1000*layer.rigidMargin_m,1000*layer.sleeveOverlap_m);
end

function t=evaluateLayout(C,kin,fitted,coarse,layout,n,j)
[i,g]=ind2sub([size(layout,1),2],j);
c=fitted{g};
c.caseDef.bundleClock_deg=layout(i,1);
c.caseDef.bundleSpread=layout(i,2);
[e,d]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',n);
t=struct('c',c,'e',e,'d',d);
other=coarse{3-g}; other.c=fitted{3-g};
if g==1
    t.layer=compare(C,kin,t,other);
else
    t.layer=compare(C,kin,other,t);
end
end

function L=compare(C,kin,a,b)
assert(isequal(a.e.poses,b.e.poses));
L=fr_compareAnkleLayers(C,kin,a.e.poses,a.d.sleeves,b.d.sleeves, ...
    a.c.caseDef.fittingLength,b.c.caseDef.fittingLength);
L.layerPass=layerPass(L);
end

function yes=layerPass(L)
yes=L.overlapPairs>0 && min([L.orderingMargin_m,L.rigidMargin_m])>=-1e-4;
end

function f=penalty(e,L,c)
limits=[e.boneMargin_m,e.leadMargin_m,e.shaftCableMargin_m,e.footMargin_m, ...
    e.distalGuideMargin_m,e.wrapMargin_m,e.bundleHardwareMargin_m, ...
    L.orderingMargin_m,L.rigidMargin_m]+1e-4;
limits=[limits,e.bundleMargin_m+1e-9,c.design.tendon-e.tendonRequired_m+1e-6];
f=1000*mean(e.noForce)+1000*sum(max(0,-limits))+ ...
    1000*double(L.overlapPairs==0)+100*L.sleeveOverlap_m;
end

function [d,trial]=solveLengths(C,kin,c,estimate,dense,branch)
cs=c.caseDef;
G=fr_ankleCoupledSetup(C,kin,cs,7);
[P,B]=fr_locForPass(C,dense.cs,c.design,kin,1);
resolved=struct('points',P,'frames',{B});
G.qg.resolvedRoute=resolved;
dense.qg.resolvedRoute=resolved;
bay=sum(vecnorm(diff(P(1:cs.distalSleeveGuideRow,:)),2,2))-2*cs.fittingLength;
n=3*numel(cs.movable);
span=fr_spans(C,kin,P,B,dense.cs.passes{1}.qStructs,'foot');
lower=max(cs.lb(n+2)/10,span.requiredSpan_m);
upper=cs.ub(n+2)/10;
assert(lower<=upper,'Tendon span exceeds the permitted length');
rest=unique([c.design.L0,estimate.L0,linspace(cs.lb(n+1)/10,min(cs.ub(n+1)/10,bay),9)]);
trial=zeros(0,5);
for L=rest
    lo=lower; hi=upper;
    a=assess(L,lo,G); b=assess(L,hi,G);
    if a(3)==0 && b(3)>0
        while hi-lo>0.00025
            mid=(lo+hi)/2;
            state=assess(L,mid,G);
            if state(3)==0; lo=mid; else; hi=mid; end
        end
    end
    for tendon=unique(min(upper,max(lower,[lo-(0:4)*.0005, ...
            lower,c.design.tendon,estimate.tendon])))
        trial(end+1,:)=assess(L,tendon,G); %#ok<AGROW>
    end
end
[~,idx]=sortrows([sum(trial(:,3:4),2),-trial(:,5)]);
fine=zeros(0,5);
for i=idx(1:min(12,numel(idx)))'
    fine(end+1,:)=assess(trial(i,1),trial(i,2),dense); %#ok<AGROW>
end
[~,i]=min(1e4*sum(fine(:,3:4),2)-fine(:,5));
d=c.design; d.L0=fine(i,1); d.tendon=fine(i,2);
fprintf('%s exact length screen: rest %.2f mm, tendon %.2f mm, slack %.1f%%, stroke %.1f%%\n', ...
    cs.name,1000*d.L0,1000*d.tendon,100*fine(i,3),100*fine(i,4));

    function row=assess(L,t,grid)
        probe=c.design; probe.L0=L; probe.tendon=t;
        m=fr_mechBuild(C,grid.cs,probe,kin,branch,grid.qg);
        slack=m.cell.PulleySlack(:);
        f=m.cell.Ftendon(:);
        stroke=m.cell.PulleyInfeasible(:) | ~isfinite(f) | (f<=1e-6 & ~slack);
        f(slack|stroke)=0;
        row=[L,t,mean(slack),mean(stroke),mean(f)];
    end
end
