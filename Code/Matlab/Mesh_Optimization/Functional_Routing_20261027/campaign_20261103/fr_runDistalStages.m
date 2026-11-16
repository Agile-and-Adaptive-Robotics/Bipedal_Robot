function fr_runDistalStages(outputDir, phase)
% Bounded, resumable group/branch solve for the distal routing stages.
C = fr_setup();
C.exactSegmentDistance = true;
C.fullFittingEnvelopes = true;
kin = fr_kinematics(C.osimStock);
switch phase
    case {'ankle','ankle_coupled','ankle_surrogate','ankle_mounts'}
        groups = {'soleus','tib_post','tib_ant','per_brev','per_long', ...
            'per_tert','flex_dig','flex_hal','ext_dig','ext_hal'};
    case 'knee'
        groups = {'med_gas','lat_gas','bifemsh','vas_int'};
    case 'ankle_extensor_layers'
        groups = {'ext_dig'};
    otherwise
        error('Unknown distal phase %s',phase);
end
layerSearch = strcmp(phase,'ankle_extensor_layers');
coupledPhase = ismember(phase,{'ankle_coupled','ankle_surrogate','ankle_mounts'}) || layerSearch;
optimizePhase = ismember(phase,{'ankle_surrogate','ankle_mounts'}) || layerSearch;
mountSearch = strcmp(phase,'ankle_mounts') || layerSearch;
if ~isfolder(outputDir); mkdir(outputDir); end
jobs = cell(2*numel(groups),2);
for i=1:numel(groups)
    jobs(2*i-1,:) = {groups{i},'rigid'};
    jobs(2*i,:) = {groups{i},'compliant'};
end
pool = gcp('nocreate');
try
    if isempty(pool)
        cluster = parcluster('local');
        jobDir = fullfile(outputDir,'parallel_jobs');
        if ~isfolder(jobDir); mkdir(jobDir); end
        cluster.JobStorageLocation = jobDir;
        pool = parpool(cluster,min([C.parallelMaxWorkers,cluster.NumWorkers,size(jobs,1)]));
    end
catch ME
    fprintf('Worker pool unavailable (%s). Resuming serially.\n',ME.message);
end
if isempty(pool)
    for k = 1:size(jobs,1)
        solveOne(C,kin,outputDir,jobs{k,1},jobs{k,2},coupledPhase,optimizePhase,mountSearch,layerSearch);
    end
else
    parfor k = 1:size(jobs,1)
        solveOne(C,kin,outputDir,jobs{k,1},jobs{k,2},coupledPhase,optimizePhase,mountSearch,layerSearch);
    end
end
end

function solveOne(C,kin,outputDir,group,branch,coupledPhase,optimizePhase,mountSearch,layerSearch)
name = [group '_' branch];
checkpoint = fullfile(outputDir,[name '.mat']);
reportFile = fullfile(outputDir,[name '.json']);
if isfile(checkpoint) && isfile(reportFile); return; end
if isfile(checkpoint)
    saved=load(checkpoint,'result','record'); result=saved.result; record=saved.record;
else
    opts=struct('quiet',true,'envelopes',{{}});
    if coupledPhase && ~strcmp(group,'soleus'); opts.coupledTravelN=5; end
    if optimizePhase
        priorFolder='ankle_coupled';
        if mountSearch; priorFolder='ankle_surrogate_600'; end
        retained=[];
        if layerSearch
            priorFolder='ankle_mounts_600';
            inner=load(fullfile(C.camp,'run_20261105_distal_stages',priorFolder, ...
                ['ext_hal_' branch '.mat']),'result','record');
            retained=inner.result.chosen;
        end
        prior=load(fullfile(C.camp,'run_20261105_distal_stages',priorFolder, ...
            [name '.mat']),'result');
        result=prior.result;
        fitCheckpoint=fullfile(outputDir,[name '_surrogate.mat']);
        if isfile(fitCheckpoint)
            priorFit=load(fitCheckpoint,'chosen','info');
            chosen=priorFit.chosen; info=priorFit.info;
        else
            fprintf('%s: starting 600 evaluations\n',name);
            [chosen,info]=fr_surrogateDistalCase(C,kin,result.chosen,branch,600,mountSearch,retained);
            tmp=[tempname(outputDir) '.mat'];
            save(tmp,'chosen','info','-v7'); movefile(tmp,fitCheckpoint);
        end
        result.surrogate=info;
        found=false;
        fprintf('%s: %d evaluations, objective %.4g -> %.4g\n', ...
            name,info.evaluations,info.initialScore,info.score);
    else
    S = load(fullfile(C.camp,'audit_20261105','reoptimization',group, ...
        ['coordinated_' branch '.mat']),'chosen');
    seed = S.chosen;
    if coupledPhase && strcmp(group,'tib_post')
        newer=fullfile(C.camp,'run_20261105_length_ladder_pilot', ...
            ['tib_post_' branch '_sizing.mat']);
        if isfile(newer)
            R=load(newer,'best');
            seed=R.best;
        end
    end
    original = seed.caseDef;
    original.allowFootLeadCorrection = any(ismember(original.rowOwning, ...
        {'talus_r','calcn_r','toes_r'}));
    direct = original;
    direct.txConfig.gain = [1 1]; direct.txConfig.tackleLineParts = [1 1];
    direct.hasTackle = false;
    result = fr_fitTravelLengths(C,direct,kin,seed.design,branch,opts);
    found = ~isempty(result.best);
    chosen = result.best;
    if ~found
        % Preserve and report the prior candidate when no sampled refit closes.
        chosen = struct('caseDef',original,'design',seed.design);
        [chosen.f,chosen.rep] = fr_objective(C,original,kin, ...
            fr_packX(original,seed.design),branch,opts);
        chosen.reeving = fr_travelReeving(max(original.txConfig.gain));
    elseif ~isfield(opts,'coupledTravelN')
        cs = chosen.caseDef;
        x = fr_packX(cs,chosen.design);
        [via,fVia] = fr_screenViaGuides(C,cs,kin,x,branch,opts,chosen.f);
        if fVia<chosen.f
            % Refit the added guide at the current transmission first.
            [fVia,rep] = fr_objective(C,via,kin,x,branch,opts);
            tx = ismember({rep.cinfo.name},{'tx_slack_cells','tx_infeasible_cells', ...
                'physical_contraction_ceiling','strain_floor'});
            if ~any([rep.cinfo(tx).violated])
                chosen.caseDef=via; chosen.rep=rep; chosen.f=fVia;
            end
        end
    end
    end
    cs=chosen.caseDef;
    dense=fr_densePilotCase(C,kin,cs);
    [~,rep]=fr_objective(C,dense,kin,fr_packX(dense,chosen.design),branch,opts);
    scorePass=1;
    if isfield(dense,'scorePass'); scorePass=dense.scorePass; end
    sweep=dense.passes{scorePass};
    [P,frames,ids] = fr_locForPass(C,cs,chosen.design,kin,scorePass);
    q=kin.qDefault;
    owning=cs.rowOwning(ids);
    for i=1:size(P,1)
        T=fr_chainT(kin,q,frames{i},owning{i});
        P(i,:)=(T(1:3,1:3)*P(i,:)'+T(1:3,4))';
    end
    gain=max(cs.txConfig.gain);
    record=struct('group',group,'branch',branch,'foundLengthFit',found, ...
        'gain',gain,'count',cs.seed.nBPA,'diameter_mm',cs.diaNum, ...
        'topologySheaves',chosen.reeving.sheaveCount,'physicalReevingVerified',gain==1, ...
        'rest_mm',chosen.design.L0*1000,'tendon_mm',chosen.design.tendon*1000, ...
        'localConstraintsPass',rep.feasible,'fractionMet',rep.fractionPosesMet, ...
        'constraints',rep.cinfo,'coordinate',sweep.sweptCoord, ...
        'angles_deg',rad2deg(sweep.th1),'robot_Nm',rep.tauPrimary,'human_Nm',rep.tauTarget, ...
        'points_m',P,'bodies',{owning},'model',C.osimRobot);
    result.chosen=chosen;
    result.dense=rep;
    if isfield(opts,'coupledTravelN')
        if layerSearch
            [record.coupled,outer]=fr_evalCoupledAnkle(C,kin,cs,chosen.design,branch,11);
            [ehl,innerDetail]=fr_evalCoupledAnkle(C,kin,retained.caseDef,retained.design,branch,11);
            assert(isequal(ehl.ankle_deg,record.coupled.ankle_deg) && ...
                isequal(ehl.second_deg,record.coupled.second_deg));
            record.layer=fr_compareAnkleLayers(C,kin,ehl.poses,outer.sleeves, ...
                innerDetail.sleeves,cs.fittingLength,retained.caseDef.fittingLength);
            record.layer.ehlLocalPass=ehl.localPass;
            record.layer.layerPass=record.layer.overlapPairs>0 && ...
                min([record.layer.orderingMargin_m,record.layer.rigidMargin_m])>=-1e-4;
            record.coupled.localPass=record.coupled.localPass && record.layer.layerPass;
            innerReport=fullfile(C.camp,'run_20261105_distal_stages','ankle_mounts_600', ...
                ['ext_hal_' branch '.json']);
            innerRoute=jsondecode(fileread(innerReport));
            record.innerRoute=struct('points_m',innerRoute.points_m,'bodies',{innerRoute.bodies});
            fprintf('%s above EHL: ordering %.2f mm, fittings %.2f mm, layer pass %d\n', ...
                name,1000*record.layer.orderingMargin_m,1000*record.layer.rigidMargin_m,record.layer.layerPass);
        else
            record.coupled=fr_evalCoupledAnkle(C,kin,cs,chosen.design,branch,11);
        end
        fprintf('%s coupled: %d slack/stroke/invalid, local pass %d\n', ...
            name,nnz(record.coupled.slack|record.coupled.stroke|record.coupled.invalid), ...
            record.coupled.localPass);
    end
    if optimizePhase
        record.surrogateEvaluations=info.evaluations;
        record.initialScore=info.initialScore;
        record.finalScore=info.score;
        if isfield(record,'coupled')
            e=record.coupled;
            record.foundLengthFit=~any(e.slack|e.stroke|e.invalid);
        else
            tx=ismember({rep.cinfo.name},{'tx_slack_cells','tx_infeasible_cells', ...
                'physical_contraction_ceiling','strain_floor'});
            record.foundLengthFit=~any([rep.cinfo(tx).violated]);
        end
    end
    if mountSearch
        record.bundleSpread=cs.bundleSpread;
        record.mountRows=[1 cs.distalSleeveGuideRow cs.shankStopRow];
    end
    tmp=[tempname(outputDir) '.mat'];
    save(tmp,'result','record','-v7');
    movefile(tmp,checkpoint);
end
tmp=[tempname(outputDir) '.json'];
fid=fopen(tmp,'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(record)); fclose(fid);
movefile(tmp,reportFile);
fprintf('%s: length fit %d, dense local constraints %d, torque %.1f%%, gain %d\n', ...
    name,record.foundLengthFit,record.localConstraintsPass,100*record.fractionMet,record.gain);
end
