function fr_refreshCoupledRoutePlots(outputDir)
% The coupled route keeps calcaneal guides omitted by the old MTP pass.
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
files=dir(fullfile(outputDir,'*.json'));
for k=1:numel(files)
    reportFile=fullfile(outputDir,files(k).name);
    r=jsondecode(fileread(reportFile));
    [~,name]=fileparts(reportFile);
    if ~isfield(r,'coupled')
        if ~strcmp(r.group,'soleus') || isfield(r,'noForceFraction'); continue; end
        S=load(fullfile(outputDir,[name '.mat']),'result');
        c=S.result.chosen;
        dense=fr_densePilotCase(C,kin,c.caseDef);
        qg=struct('passIdx',1,'qList',{dense.passes{1}.qStructs});
        m=fr_mechBuild(C,dense,c.design,kin,r.branch,qg);
        r.noForceFraction=mean(m.PulleySlack|m.PulleyInfeasible|~isfinite(m.Ftendon)|m.Ftendon<=1e-6);
        writeRecord(outputDir,reportFile,r);
        continue
    end
    S=load(fullfile(outputDir,[name '.mat']),'result');
    c=S.result.chosen;
    if ~isfield(r.coupled,'noForce')
        r.coupled=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,r.branch,11);
        if isfield(r,'layer')
            r.coupled.localPass=r.coupled.localPass && r.layer.layerPass;
        end
    end
    G=fr_ankleCoupledSetup(C,kin,c.caseDef,11);
    [P,frames,ids]=fr_locForPass(C,G.cs,c.design,kin,1);
    native=G.cs.rowOwning(ids);
    for j=1:numel(ids)
        T=fr_chainT(kin,kin.qDefault,frames{j},native{j});
        P(j,:)=(T(1:3,1:3)*P(j,:)'+T(1:3,4))';
    end
    r.points_m=P;
    r.bodies=native;
    writeRecord(outputDir,reportFile,r);
end
end

function writeRecord(outputDir,file,r)
temp=[tempname(outputDir) '.json'];
fid=fopen(temp,'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(r)); fclose(fid);
movefile(temp,file,'f');
end
