function fr_check_metadata_20261105()
% Unresolved sheave metadata must not alter a saved force/torque result.
C=fr_setup();kin=fr_kinematics(C.osimStock);
root=fullfile(C.camp,'audit_20261105');
S=load(fullfile(root,'reoptimization','bifemsh','coordinated_compliant.mat'),'chosen');a=S.chosen;
[~,rep]=fr_objective(C,a.caseDef,kin,fr_packX(a.caseDef,a.design),'compliant',struct('quiet',true));
error=max(abs(rep.tauPrimary-a.globalRep.tauPrimary));assert(error<1e-8,'Metadata changed the saved torque');
[m,~]=fr_mechBuild(C,a.caseDef,a.design,kin,'compliant',struct('passIdx',1,'qList',{a.caseDef.passes{1}.qStructs}));
assert(isnan(m.cfg.sheaveCountPhysical));
fid=fopen(fullfile(root,'metadata_force_verification.json'),'w');fwrite(fid,jsonencode(struct('status','PASS','max_saved_torque_error_Nm',error,'physical_sheave_count','UNRESOLVED')));fclose(fid);
fprintf('PASS: unresolved physical sheave metadata; saved torque unchanged (%.3g Nm).\n',error);
end
