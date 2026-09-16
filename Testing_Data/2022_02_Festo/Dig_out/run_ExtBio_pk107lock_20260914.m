% run_ExtBio_pk107lock_20260914.m
% Wrapper: rerun Dig_ExtPin_frontBio.m and save FULL workspace (Ben 2026-09-12).
run('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/Dig_ExtPin_frontBio.m');
save('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo/Dig_out/Dig_ExtPin_frontBio_pk107lock_20260914.mat');  % full workspace
fprintf('SAVED FULL WORKSPACE: Dig_out/Dig_ExtPin_frontBio_pk107lock_20260914.mat\n');
