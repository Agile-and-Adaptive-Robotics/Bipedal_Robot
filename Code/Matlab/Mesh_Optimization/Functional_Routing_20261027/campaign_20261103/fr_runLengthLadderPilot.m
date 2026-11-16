function fr_runLengthLadderPilot(outputDir, groups)
% Separate checkpoint and figures for each group and branch.
C = fr_setup();
kin = fr_kinematics(C.osimStock);
C.exactSegmentDistance = true;
C.fullFittingEnvelopes = true;
if ~isfolder(outputDir); mkdir(outputDir); end
opts = struct('quiet',true,'envelopes',{{}});
for i = 1:numel(groups)
    group = groups{i};
    for b = {'compliant','rigid'}
        branch = b{1};
        key = [group '_' branch];
        checkpoint = fullfile(outputDir,[key '.mat']);
        if isfile(checkpoint)
            previous = load(checkpoint,'result');
            result = previous.result;
        else
        src = fullfile(C.camp,'audit_20261105','reoptimization',group, ...
            ['coordinated_' branch '.mat']);
        S = load(src,'chosen');
        cs = S.chosen.caseDef;
        cs.allowFootLeadCorrection = any(ismember(cs.rowOwning,{'talus_r','calcn_r','toes_r'}));
        cs.txConfig.gain = [1 1];
        cs.txConfig.tackleLineParts = [1 1];
        cs.hasTackle = false;
        result = fr_fitTravelLengths(C,cs,kin,S.chosen.design,branch,opts);
        temp = [tempname(outputDir) '.mat'];
        save(temp,'result','-v7');
        movefile(temp,checkpoint);
        end
        if isempty(result.best); continue; end
        chosen = result.best;
        rep = chosen.rep;
        dense = fr_densePilotCase(C,kin,chosen.caseDef);
        [~,denseRep] = fr_objective(C,dense,kin,fr_packX(dense,chosen.design),branch,opts);
        save(fullfile(outputDir,[key '_dense.mat']),'denseRep','-v7');
        fprintf('%s: BPA %.1f mm, tendon %.1f mm, gain %d, topology sheaves %d, torque %.1f%%, %s\n', ...
            key,chosen.design.L0*1000,chosen.design.tendon*1000,max(chosen.caseDef.txConfig.gain), ...
            chosen.reeving.sheaveCount,rep.fractionPosesMet*100,rep.topViolator);
        fprintf('%s dense: torque %.1f%%, failed constraints %s\n',key, ...
            denseRep.fractionPosesMet*100,strjoin({denseRep.cinfo([denseRep.cinfo.violated]).name},', '));
        if strcmp(group,'tib_post')
            folder = fullfile(outputDir,branch);
            fr_previewFootLead(folder,chosen,branch);
        end
        f = figure('Visible','off','Color','w');
        ax = axes(f);
        scorePass = 1;
        if isfield(chosen.caseDef,'scorePass'); scorePass = chosen.caseDef.scorePass; end
        p = chosen.caseDef.passes{scorePass};
        plot(ax,rad2deg(p.th1),rep.tauTarget,'--','Color',[.45 .45 .45],'LineWidth',1.5);
        hold(ax,'on');
        plot(ax,rad2deg(p.th1),rep.tauPrimary,'-o','Color',[0 .35 .8],'LineWidth',1.5);
        xlabel(ax,[strrep(p.sweptCoord,'_',' ') ' (deg)']);
        ylabel(ax,'Signed torque (N m)');
        title(ax,strrep(key,'_',' ')); grid(ax,'on');
        legend(ax,{'Human target','Robot candidate'},'Location','best');
        exportgraphics(ax,fullfile(outputDir,[key '_torque.png']),'Resolution',180);
        close(f);
    end
end
end

