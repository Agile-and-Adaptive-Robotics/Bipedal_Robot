function fr_checkEDLWrapXi3(outputDir)
% Compare the exact same saved geometry and length fit with vessel Xi3 off.
C=fr_setup(); C.exactSegmentDistance=true; C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
S=load(fullfile(outputDir,'ext_dig_compliant.mat'),'result');
c=S.result.chosen;
e=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',11);
C.Xi.xi3Override=0;
off=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,'compliant',11);
assert(~isfield(off,'bpaBendLossMax_m') || off.bpaBendLossMax_m==0);
fprintf(['EDL bend-only Xi3 comparison: peak ankle torque difference %.3f N m, ' ...
    'peak MTP torque difference %.3f N m, no-force %.1f%% with / %.1f%% without\n'], ...
    max(abs(e.ankle_Nm-off.ankle_Nm)),max(abs(e.second_Nm-off.second_Nm)), ...
    100*mean(e.noForce),100*mean(off.noForce));
save(fullfile(outputDir,'edl_bend_xi3_check.mat'),'e','off','-v7');
end
