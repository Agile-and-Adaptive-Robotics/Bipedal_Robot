tblDir = 'D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables';
files = {'PinnedFlexor_Preopt_GoF.xlsx','PinnedFlexor_GoF.xlsx','PinnedExtensor_GoF.xlsx','AdoptedCorrectionTerms.xlsx'};
for i = 1:numel(files)
    f = fullfile(tblDir, files{i});
    fprintf('=== %s\n', files{i});
    try
        sh = sheetnames(f);
        for s = 1:numel(sh)
            fprintf('-- sheet: %s\n', sh{s});
            C = readcell(f, 'Sheet', sh{s});
            disp(C);
        end
    catch me
        fprintf('READ FAILED: %s\n', me.message);
    end
end
