function fr_refineTibialStopPilot(outputDir)
% Search stop positions on a local tibial patch; preserve the foot anchors.
C = fr_setup();
C.exactSegmentDistance = true;
C.fullFittingEnvelopes = true;
kin = fr_kinematics(C.osimStock);
opts = struct('quiet',true,'envelopes',{{}});
for br = {'compliant','rigid'}
    branch = br{1};
    S = load(fullfile(outputDir,['tib_post_' branch '.mat']),'result');
    assert(~isempty(S.result.best), 'Stop refinement requires a direct length fit');
    seed = S.result.best;
    cs = seed.caseDef;
    firstFoot = find(ismember(cs.rowOwning,{'calcn_r','talus_r','toes_r'}),1);
    stopRow = firstFoot-1;
    assert(strcmp(cs.rowOwning{stopRow},'tibia_r'));
    assert(~any([cs.movable.row]==stopRow),'Pilot expects a fixed initial tibial stop');
    patch = fr_patch(C,'Tibia','tibia_r',cs.seedLoc(stopRow,:),25);
    best = seed;
    offsets = [0 0; .6 0; -.6 0; 0 .6; 0 -.6];
    candidates = cell(size(offsets,1),1);
    clearance = -inf(size(candidates));
    for j = 1:size(offsets,1)
        candidate = cs;
        point = patch.project(offsets(j,1),offsets(j,2),C.standoff_mm(2));
        candidate.seedLoc(stopRow,:) = point(:)';
        candidate.rows{stopRow}.loc = point(:)';
        candidates{j} = candidate;
        margin = inf;
        for p = 1:numel(candidate.passes)
            [P,F,I] = fr_locForPass(C,candidate,seed.design,kin,p,true);
            margin = min(margin,fr_tibialLeadMargin(C,kin,P,F, ...
                candidate.rowOwning(I),candidate.passes{p}.qStructs,.002));
        end
        clearance(j) = margin;
    end
    [~,order] = sort(clearance,'descend');
    for j = order(1:min(3,numel(order)))'
        if clearance(j)<-1e-4; continue; end
        result = fr_fitTravelLengths(C,candidates{j},kin,seed.design,branch,opts,1);
        if ~isempty(result.best) && result.best.f<best.f
            best = result.best;
        end
    end
    stopMove = norm(best.caseDef.seedLoc(stopRow,:)-cs.seedLoc(stopRow,:));
    fprintf('Stop refinement %s: moved %.2f mm, torque %.1f%%, objective %.5g -> %.5g; %s\n', ...
        branch,stopMove*1000,100*best.rep.fractionPosesMet,seed.f,best.f,best.rep.topViolator);
    save(fullfile(outputDir,['tib_post_' branch '_stop.mat']),'best','clearance','stopMove','-v7');
    fr_previewFootLead(fullfile(outputDir,[branch '_stop']),best,branch);
end
end
