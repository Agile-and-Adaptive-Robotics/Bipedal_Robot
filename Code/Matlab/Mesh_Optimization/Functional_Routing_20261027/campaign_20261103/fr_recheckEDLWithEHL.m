function fr_recheckEDLWithEHL(folder,pairFolder,ehlCount)
% Recheck the separate fixed-length EDL paths against the final 20 mm EHL.
C=fr_setup();
C.exactSegmentDistance=true; C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
S=load(fullfile(pairFolder,'independent_pair_dense.mat'),'chosen','c','record');
pair=S.chosen;
c=S.c;
name='final_ehl_recheck';
if nargin<3
    Q=load(fullfile(folder,'ext_hal_compliant.mat'),'result');
    ehl=Q.result.chosen;
else
    Q=load(fullfile(folder,sprintf('ext_hal_compliant_%dBPA.mat',ehlCount)),'trial');
    ehl=Q.trial.chosen;
    name=sprintf('ehl_%dBPA_recheck',ehlCount);
end
[E,D]=fr_evalCoupledAnkle(C,kin,ehl.caseDef,ehl.design,'compliant',21);
parts={pair.a,pair.b};
clearance=cell(1,2);
for k=1:2
    b=parts{k};
    assert(isequal(E.poses,b.e.poses));
    L=fr_compareAnkleLayers(C,kin,E.poses,b.d.sleeves,D.sleeves, ...
        c.caseDef.fittingLength,ehl.caseDef.fittingLength);
    A=b.d.loc(c.caseDef.shankStopRow,:);
    B=D.loc(ehl.caseDef.shankStopRow,:);
    diameter=hypot(C.bowden.stopAllowance_mm,C.bowden.housingOD_mm)/1000;
    L.stopClearance_m=norm(A-B)-diameter-C.clearanceEnvelope_m;
    L.wireMargin_m=fr_ankleWireClearance(kin,E.poses,b.d,c.caseDef.shankStopRow, ...
        D,ehl.caseDef.shankStopRow);
    L.pass=L.overlapPairs>0 && min([L.rigidMargin_m,L.orderingMargin_m, ...
        L.stopClearance_m,L.wireMargin_m])>=-1e-4;
    clearance{k}=L;
end
result=struct('ehlCount',ehl.caseDef.seed.nBPA,'ehlLocalPass',E.localPass, ...
    'edlLocalPass',S.record.localConstraintsPass,'clearance',{clearance});
result.interfacePass=all(cellfun(@(v)v.pass,clearance));
result.pairedPass=result.ehlLocalPass && result.edlLocalPass && result.interfacePass;
save(fullfile(pairFolder,[name '.mat']),'result','E','D','-v7');
f=fopen(fullfile(pairFolder,[name '.json']),'w'); assert(f>=0);
fprintf(f,'%s\n',jsonencode(result,'PrettyPrint',true)); fclose(f);
fprintf('EDL with final EHL: fittings %+.2f / %+.2f mm, interface pass %d, EHL local pass %d.\n', ...
    1000*clearance{1}.rigidMargin_m,1000*clearance{2}.rigidMargin_m, ...
    result.interfacePass,result.ehlLocalPass);
end
