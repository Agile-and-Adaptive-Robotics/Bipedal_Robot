function target = refresh_sns_library_icons(target)
% Export a display-only library copy after fixing the shared icon emitter.
% Dynamics, mask values, and source SNS_Library.slx are preserved.
root=fileparts(mfilename('fullpath'));
while ~endsWith(root,'Bipedal_Robot'), root=fileparts(root); end
src=fullfile(root,'Code','Matlab','SNS_Simscape');
addpath(src);
if nargin<1
    target=fullfile(root,'tmp','sns_library_icons_20261116','SNS_Library.slx');
end
assert(~strcmpi(char(java.io.File(target).getCanonicalPath()), ...
    char(java.io.File(fullfile(src,'SNS_Library.slx')).getCanonicalPath())), ...
    'Choose a separate output copy; the source library must remain unchanged.');
folder=fileparts(target); if ~isfolder(folder),mkdir(folder);end
load_system(fullfile(src,'SNS_Library.slx'));
cleanup=onCleanup(@()close_system('SNS_Library',0));
set_param('SNS_Library','Lock','off');
families={'BPAForce','BPA';'BioMuscle','BIO';'BPA_10mm','BPA 10'; ...
    'BPA_20mm','BPA 20';'BPA_40mm','BPA 40'};
for k=1:size(families,1)
    mask=Simulink.Mask.get(['SNS_Library/' families{k,1}]);
    mask.Display=sns_icon_muscle(families{k,2});
end
set_param('SNS_Library','Lock','on');
save_system('SNS_Library',target);
fprintf('Display-only library copy: %s\n',target);
end
