root='D:/Github/Bipedal_Robot';addpath(genpath(fullfile(root,'Code/Matlab')));addpath(fullfile(root,'Testing_Data/2022_02_Festo'));
o=fullfile(root,'Documentation/Reports and Papers/Dissertation/Archive/appendix_c_corrections_20261028');
out=struct;
files={'KneeFlx_10mm_42cm.mat','KneeFlx_20mm_42cm.mat','KneeFlx_20mm_42cm.mat'};vars={'Bifemsh_Pam','Bifemsh_Pam3','Bifemsh_Pam2'};
for i=1:3,s=load(fullfile(root,'Testing_Data/2022_02_Festo',files{i}),vars{i});k=s.(vars{i});out.flx(i)=struct('rest',k.RestingL,'tendon',k.TendonL,'diameter',k.Diameter,'pressure',k.Pressure,'Fmax',k.Fmax);end
s=load(fullfile(root,'Testing_Data/2022_02_Festo/Plot_KneeExt10mm_52cm.mat'),'Vas_Pam');k=s.Vas_Pam;out.ext=struct('rest',k.RestingL,'tendon',k.TendonL,'diameter',k.Diameter,'pressure',k.Pressure,'Fmax',k.Fmax);
s=load(fullfile(root,'Code/Matlab/Mesh_Optimization/Results/Bifemsh_20mm_Result.mat'),'xBest');out.finalFlexor_xBest=s.xBest;
s=load(fullfile(root,'Code/Matlab/Mesh_Optimization/Results/Vas_Pam_20mm_Result_cam_20261014_0416.mat'),'predBest','ctxCamBest');out.finalExtensor=struct('rest',s.predBest.rest,'tendon',s.predBest.tendon);out.geo=s.ctxCamBest.geo;
out.T0=s.ctxCamBest.T(:,:,find(abs(s.ctxCamBest.phiD)==min(abs(s.ctxCamBest.phiD)),1));
f=fopen(fullfile(o,'source_constants.json'),'w');fprintf(f,'%s',jsonencode(out,PrettyPrint=true));fclose(f);

