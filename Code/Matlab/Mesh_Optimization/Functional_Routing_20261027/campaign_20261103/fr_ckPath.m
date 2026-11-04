function p = fr_ckPath(name, branch, ckDir)
%FR_CKPATH Checkpoint file path for one (group, branch).
p = fullfile(ckDir, sprintf('%s_%s.mat', name, branch));
end
