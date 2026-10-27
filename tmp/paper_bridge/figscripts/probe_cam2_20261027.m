% probe_cam2_20261027.m -- dump remaining cam-mat reference fields (read-only)
R = 'D:\GitHub\Bipedal_Robot';
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
fprintf('xBestCam  = [%s]\n', num2str(C.xBestCam(:).', 8));
fprintf('xFullBest = [%s]\n', num2str(C.xFullBest(:).', 8));
fprintf('xFixed    = [%s]\n', num2str(C.xFixed(:).', 8));
fprintf('xRef8     = [%s]\n', num2str(C.xRef8(:).', 8));
fprintf('cBest     = [%s]\n', num2str(C.cBest(:).', 8));
fprintf('tendonRef = %g   thetaP6BaseD = [%s]\n', C.tendonRef, num2str(C.thetaP6BaseD(:).', 6));
disp('seedBase ='); disp(C.seedBase);
disp('refPts ='); disp(C.refPts);
f = fieldnames(C.routeInfoBase);
fprintf('routeInfoBase fields: %s\n', strjoin(f', ', '));
for k = 1:numel(f)
  v = C.routeInfoBase.(f{k});
  if isnumeric(v) && numel(v) <= 12
    fprintf('  %s = [%s]\n', f{k}, num2str(v(:).', 6));
  end
end
fprintf('routeInfoBase.route =\n');
if isfield(C.routeInfoBase, 'route'); disp(size(C.routeInfoBase.route)); end
fprintf('requiredBest = [%s]\n', num2str(C.requiredBest(:).', 6));
fprintf('CAM2 PROBE DONE\n');
