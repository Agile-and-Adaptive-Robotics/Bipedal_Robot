function fr_screenTibialSizingPilot(outputDir)
% Bounded follow-up: compare direct transmissions at the validated stop.
C = fr_setup();
C.exactSegmentDistance = true;
C.fullFittingEnvelopes = true;
kin = fr_kinematics(C.osimStock);
opts = struct('quiet',true,'envelopes',{{}});
for br = {'compliant','rigid'}
    branch = br{1};
    checkpoint = fullfile(outputDir,['tib_post_' branch '_sizing.mat']);
    if isfile(checkpoint); continue; end
    seed = load(fullfile(outputDir,['tib_post_' branch '_stop.mat']),'best');
    best = seed.best;
    variants = [10 2; 10 3; 20 1; 20 2; 20 3];
    trials = cell(size(variants,1),1);
    for i = 1:size(variants,1)
        cs = seed.best.caseDef;
        cs.diaNum = variants(i,1);
        cs.diaStr = num2str(cs.diaNum);
        cs.seed.nBPA = variants(i,2);
        cs.txConfig.nPulleyBPA = [cs.seed.nBPA 1];
        cap = 1-C.LminFraction.(['d' cs.diaStr]);
        cs.bundleSpread = fr_diameter(C,cs.diaStr,cap)/fr_diameter(C,cs.diaStr,0);
        result = fr_fitTravelLengths(C,cs,kin,seed.best.design,branch,opts,1);
        trials{i} = result;
        if isempty(result.best); continue; end
        c = result.best;
        fprintf('Sizing %s D%d n%d: objective %.3f, torque %.1f%%, %s\n', ...
            branch,cs.diaNum,cs.seed.nBPA,c.f,100*c.rep.fractionPosesMet,c.rep.topViolator);
        if c.rep.feasible && c.f < best.f
            dense = fr_densePilotCase(C,kin,c.caseDef);
            [~,r] = fr_objective(C,dense,kin,fr_packX(dense,c.design),branch,opts);
            if r.feasible; best = c; end
        end
    end
    dense = fr_densePilotCase(C,kin,best.caseDef);
    [~,denseRep] = fr_objective(C,dense,kin,fr_packX(dense,best.design),branch,opts);
    save(checkpoint,'best','trials','denseRep','-v7');
    fprintf('Best %s: D%s n%d, signed torque meet %.1f%%, constraints feasible %d\n', ...
        branch,best.caseDef.diaStr,best.caseDef.seed.nBPA, ...
        100*denseRep.fractionPosesMet,denseRep.feasible);
    folder = fullfile(outputDir,[branch '_sizing']);
    fr_previewFootLead(folder,best,branch);
    f = figure('Visible','off','Color','w');
    ax = axes(f);
    x = rad2deg(dense.passes{1}.th1);
    plot(ax,x,denseRep.tauTarget,'--','Color',[.45 .45 .45],'LineWidth',2);
    hold(ax,'on');
    plot(ax,x,denseRep.tauPrimary,'-','Color',[0 .35 .8],'LineWidth',2);
    grid(ax,'on'); xlabel(ax,'Subtalar angle (deg)'); ylabel(ax,'Signed torque (N m)');
    title(ax,sprintf('Tibialis posterior | %s | %d x %s mm BPA', ...
        branch,best.caseDef.seed.nBPA,best.caseDef.diaStr));
    legend(ax,{'Human target (interpolated)','Robot candidate'},'Location','best');
    exportgraphics(ax,fullfile(folder,'tibialis_posterior_torque.png'),'Resolution',180);
    close(f);
end
end
