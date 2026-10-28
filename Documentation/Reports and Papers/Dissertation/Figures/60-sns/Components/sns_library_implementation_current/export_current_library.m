root='D:/Github/Bipedal_Robot';here=fullfile(root,'Code/Matlab/SNS_Simscape');out=fullfile(root,'tmp/diss_history_native');addpath(genpath(here));load_system(fullfile(here,'SNS_Library.slx'));set_param('SNS_Library','Lock','off');
names={'NonSpikingNeuron','NonSpikingSynapse','SynSum','SpikingLIFNeuron','IaMuscleSpindle','IbGolgiTendon','MuscleActivation','BioMuscle','BPAForce','BPA_10mm','BPA_20mm','BPA_40mm'};
for i=1:numel(names)
 x=60+mod(i-1,4)*240;y=60+floor((i-1)/4)*185;set_param(['SNS_Library/' names{i}],'Position',[x y x+170 y+115],'FontSize','14');
end
set_param('SNS_Library','PaperOrientation','landscape','PaperType','A4','PaperPositionMode','auto');print('-sSNS_Library','-dpdf',fullfile(out,'sns_library_grid_current.pdf'));close_system('SNS_Library',0);
