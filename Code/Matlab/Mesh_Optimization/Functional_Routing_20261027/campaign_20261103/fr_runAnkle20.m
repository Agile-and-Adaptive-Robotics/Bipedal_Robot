function fr_runAnkle20(outputDir,budget,workers,groups,branches)
% Refit the tibia-down campaign with separate 20 mm actuator load paths.
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
if nargin<2; budget=600; end
if nargin<3; workers=C.parallelMaxWorkers; end
if nargin<4
    groups={'soleus','tib_post','tib_ant','per_brev','per_long','per_tert', ...
        'flex_dig','flex_hal','ext_hal','ext_dig'};
end
if nargin<5; branches={'rigid','compliant'}; end
if ~isfolder(outputDir); mkdir(outputDir); end
[G,B]=ndgrid(1:numel(groups),1:numel(branches));
jobs=[G(:),B(:)];
pool=[];
if workers>1
    try
        pool=gcp('nocreate');
        if isempty(pool)
            cluster=parcluster('local');
            cluster.JobStorageLocation=fullfile(outputDir,'parallel_jobs');
            if ~isfolder(cluster.JobStorageLocation); mkdir(cluster.JobStorageLocation); end
            pool=parpool(cluster,min([workers,cluster.NumWorkers,size(jobs,1)]));
        end
    catch ME
        fprintf('Worker startup failed: %s. Continuing serially.\n',ME.message);
    end
end
if isempty(pool)
    for j=1:size(jobs,1)
        solveGroup(C,kin,outputDir,groups{jobs(j,1)},branches{jobs(j,2)},budget);
    end
else
    parfor j=1:size(jobs,1)
        solveGroup(C,kin,outputDir,groups{jobs(j,1)},branches{jobs(j,2)},budget);
    end
end
end

function solveGroup(C,kin,folder,group,branch,budget)
name=[group '_' branch];
checkpoint=fullfile(folder,[name '.mat']);
reportFile=fullfile(folder,[name '.json']);
if isfile(checkpoint) && isfile(reportFile); return; end
if isfile(checkpoint)
    S=load(checkpoint,'record'); record=S.record;
else
    rng(sum(double(name)),'twister');
    old=fullfile(C.camp,'run_20261105_distal_stages');
    source=fullfile(old,'ankle_mounts_600',[name '.mat']);
    S=load(source,'result');
    seed=S.result.chosen;
    seed.caseDef.diaNum=20;
    seed.caseDef.diaStr='20';
    mono=strcmp(group,'soleus');
    opts=struct('quiet',true,'envelopes',{{}});
    retained=[];
    if strcmp(group,'ext_dig')
        S=load(fullfile(old,'ankle_mounts_600',['ext_hal_' branch '.mat']),'result');
        retained=S.result.chosen;
        retained.caseDef.diaNum=20; retained.caseDef.diaStr='20';
        G=fr_ankleCoupledSetup(C,kin,retained.caseDef,7);
        retained.design=fr_refitMountLengths(C,kin,retained.caseDef,retained.design,branch,G);
        pilot=fullfile(old,'ankle_distal_sizing_pilot',['ext_hal_' branch '.mat']);
        if isfile(pilot)
            S=load(pilot,'result'); retained=S.result.chosen;
        end
    end
    counts=unique([seed.caseDef.seed.nBPA,1,2,3],'stable');
    each=floor(budget/numel(counts));
    trials=cell(size(counts));
    rank=inf(numel(counts),4);
    for k=1:numel(counts)
        count=counts(k);
        file=fullfile(folder,sprintf('%s_%dBPA.mat',name,count));
        if isfile(file)
            S=load(file,'trial'); trial=S.trial;
        else
            c=seed;
            c.caseDef.seed.nBPA=count;
            c.caseDef.txConfig.nPulleyBPA=[1 1];
            c.caseDef.bundleSpread=max(1.9,c.caseDef.bundleSpread);
            fprintf('%s: %d x 20 mm, starting %d evaluations\n',name,count,each);
            [c,info]=fr_surrogateDistalCase(C,kin,c,branch,each,~mono,retained,strcmp(group,'ext_hal'));
            trial=struct('chosen',c,'info',info);
            if mono
                dense=fr_densePilotCase(C,kin,c.caseDef);
                [~,rep]=fr_objective(C,dense,kin,fr_packX(dense,c.design),branch,opts);
                m=fr_mechBuild(C,dense,c.design,kin,branch,struct('passIdx',1,'qList',{dense.passes{1}.qStructs}));
                trial.dense=rep;
                trial.noForce=mean(m.PulleySlack|m.PulleyInfeasible|~isfinite(m.Ftendon)|m.Ftendon<=1e-6);
                trial.rank=[~rep.feasible,trial.noForce,0,info.score];
            else
                [e,d]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,branch,21);
                trial.coupled=e;
                trial.detail=d;
                if ~isempty(retained)
                    [inner,I]=fr_evalCoupledAnkle(C,kin,retained.caseDef,retained.design,branch,21);
                    assert(isequal(e.poses,inner.poses));
                    L=fr_compareAnkleLayers(C,kin,e.poses,d.sleeves,I.sleeves, ...
                        c.caseDef.fittingLength,retained.caseDef.fittingLength);
                    L.layerPass=L.overlapPairs>0 && min([L.rigidMargin_m,L.orderingMargin_m])>=-1e-4;
                    trial.layer=L;
                    trial.coupled.localPass=e.localPass && L.layerPass;
                end
                margins=[e.boneMargin_m,e.leadMargin_m,e.footMargin_m,e.distalGuideMargin_m, ...
                    e.wrapMargin_m,e.shaftCableMargin_m,e.bundleMargin_m,e.bundleHardwareMargin_m];
                if isfield(trial,'layer')
                    margins=[margins,trial.layer.rigidMargin_m,trial.layer.orderingMargin_m];
                end
                trial.rank=[~trial.coupled.localPass,mean(e.noForce),sum(max(0,-margins)),info.score];
            end
            saveAtomic(file,'trial',trial);
            fprintf('%s: %d x 20 complete, pass %d, no force %.1f%%\n', ...
                name,count,~trial.rank(1),100*trial.rank(2));
        end
        trials{k}=trial;
        rank(k,:)=trial.rank;
    end
    [~,order]=sortrows(rank);
    best=trials{order(1)};
    c=best.chosen; cs=c.caseDef;
    record=jsondecode(fileread(fullfile(old,'ankle_mounts_600',[name '.json'])));
    record=rmfield(record,intersect(fieldnames(record),{'surrogateEvaluations','initialScore','finalScore'}));
    record.count=cs.seed.nBPA;
    record.diameter_mm=20;
    record.gain=max(cs.txConfig.gain);
    record.topologySheaves=record.count*fr_travelReeving(record.gain).sheaveCount;
    record.physicalReevingVerified=record.gain==1;
    record.independentLoadPaths=true;
    record.rest_mm=1000*c.design.L0;
    record.tendon_mm=1000*c.design.tendon;
    record.bundleSpread=cs.bundleSpread;
    record.trialEvaluations=sum(cellfun(@(t)t.info.evaluations,trials));
    if mono
        rep=best.dense;
        record.fractionMet=rep.fractionPosesMet;
        record.constraints=rep.cinfo;
        record.robot_Nm=rep.tauPrimary;
        record.human_Nm=rep.tauTarget;
        dense=fr_densePilotCase(C,kin,cs);
        record.angles_deg=rad2deg(dense.passes{1}.th1);
        record.noForceFraction=best.noForce;
        [P,frames,ids]=fr_locForPass(C,cs,c.design,kin,1);
        record.localConstraintsPass=rep.feasible;
    else
        record.legacySweep=false;
        record.coupled=best.coupled;
        record.localConstraintsPass=best.coupled.localPass;
        record.mountRows=[1 cs.distalSleeveGuideRow cs.shankStopRow];
        if isfield(best,'layer'); record.layer=best.layer; end
        G=fr_ankleCoupledSetup(C,kin,cs,2);
        [P,frames,ids]=fr_locForPass(C,G.cs,c.design,kin,1);
    end
    record.foundLengthFit=best.rank(2)==0;
    bodies=cs.rowOwning(ids);
    for r=1:numel(ids)
        T=fr_chainT(kin,kin.qDefault,frames{r},bodies{r});
        P(r,:)=(T(1:3,1:3)*P(r,:)'+T(1:3,4))';
    end
    record.points_m=P;
    record.bodies=bodies;
    tmp=[tempname(folder) '.mat'];
    result=struct('chosen',c,'trials',{trials}); %#ok<NASGU>
    save(tmp,'result','record','-v7');
    movefile(tmp,checkpoint);
end
tmp=[tempname(folder) '.json'];
fid=fopen(tmp,'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
movefile(tmp,reportFile);
end

function saveAtomic(file,field,value)
S=struct(); S.(field)=value;
tmp=[tempname(fileparts(file)) '.mat'];
save(tmp,'-struct','S','-v7');
movefile(tmp,file);
end
