function fr_explainDuplicateAnkle(outputDir)
% Closure term audit for the unchanged one- and two-BPA ankle routes.
C=fr_setup();
kin=fr_kinematics(C.osimStock);
stage=fullfile(C.camp,'run_20261105_distal_stages');
groups={'ext_dig','ext_hal'}; folders={'ankle_guide_arrow','ankle_mounts_600'};
rows={};
for g=1:2
    S=load(fullfile(stage,folders{g},[groups{g} '_compliant.mat']),'result');
    c=S.result.chosen;
    for count=1:2
        cs=c.caseDef;
        cs.seed.nBPA=count;
        cs.txConfig.nPulleyBPA=[count 1];
        G=fr_ankleCoupledSetup(C,kin,cs,11);
        for mode=1:4
            context=C;
            branch='compliant';
            label={'full compliance','no insertion compliance','origin compliance only','rigid'};
            if mode==2; context.Xi.noInsBracket=true; end
            if mode==3; context.Xi.noCompliance=true; end
            if mode==4; branch='rigid'; end
            m=fr_mechBuild(context,G.cs,c.design,kin,branch,G.qg);
            x=m.cell;
            i=find(abs(G.ankle_deg+1)<1e-9 & G.second_deg==0,1);
            assert(~isempty(i));
            rest=m.rest;
            cap=rest-m.kmax;
            initial=rest-x.MuscleLength(:)+x.Xi0+x.TendonL+2*x.FittingLength;
            initial=initial+x.deltaRoute(:);
            origin=x.cbBracket(:);
            added=x.DeflClosure_cAdd(:);
            factor=added-origin-count/x.kSpr;
            s=linspace(0,cap,101);
            F=x.Fmax*festo4(cs.diaNum,s/cap,C.pressure);
            residual=s-initial+factor*F;
            zeroRoots=sum(residual(:,1)>0 & min(residual,[],2)<=0);
            force=x.Ftendon(:);
            rows{end+1}=struct('group',groups{g},'count',count,'mode',label{mode}, ... %#ok<AGROW>
                'slack_pct',100*mean(x.PulleySlack(:)), ...
                'stroke_pct',100*mean(x.PulleyInfeasible(:)), ...
                'minimum_force_N',min(force),'maximum_force_N',max(force), ...
                'ankle_deg',G.ankle_deg(i),'mtp_deg',G.second_deg(i), ...
                'initialDemand_mm',1000*initial(i),'addedCompliance_mm',1000*added(i)*F(1), ...
                'originCompliance_mm',1000*origin(i)*F(1), ...
                'tendonStretch_mm',1000*count*F(1)/x.kSpr, ...
                'zeroContractionResidual_mm',1000*residual(i,1), ...
                'minimumResidual_mm',1000*min(residual(i,:)), ...
                'positiveStartInteriorRootCount',zeroRoots,'forceAtPose_N',force(i));
        end
    end
end
T=struct2table(vertcat(rows{:}));
writetable(T,fullfile(outputDir,'duplicate_closure_audit.csv'));
disp(T(:,{'group','count','mode','slack_pct','stroke_pct','initialDemand_mm', ...
    'addedCompliance_mm','zeroContractionResidual_mm','minimumResidual_mm', ...
    'positiveStartInteriorRootCount','forceAtPose_N'}));
end
