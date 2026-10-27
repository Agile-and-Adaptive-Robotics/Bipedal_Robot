% probe_cam_20261027.m -- dump cam design-of-record values (read-only)
R = 'D:\GitHub\Bipedal_Robot';
C = load(fullfile(R, 'Code', 'Matlab', 'Mesh_Optimization', 'Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
fprintf('xBestCam = [%s]\n', num2str(C.xBestCam(:).', 8));
fprintf('xFullBest = [%s]\n', num2str(C.xFullBest(:).', 8));
fprintf('xFixed = [%s]\n', num2str(C.xFixed(:).', 8));
fprintf('stamp = %s\n', C.stamp);
if isfield(C, 'fBestCam'); fprintf('fBestCam = %g\n', C.fBestCam); end
if isfield(C, 'marginBest'); fprintf('marginBest = %g\n', C.marginBest); end
if isfield(C, 'marginReqBest'); fprintf('marginReqBest = %g\n', C.marginReqBest); end
if isfield(C, 'XiUsed'); fprintf('XiUsed = [%s]\n', num2str(C.XiUsed(:).', 8)); end
p = C.predBest;
fprintf('predBest.p1 = [%s]  pEnd = [%s]\n', num2str(p.p1, 6), num2str(p.pEnd, 6));
fprintf('predBest.rest = %g  tendon = %g  BPAcount = %g\n', p.rest, p.tendon, p.BPAcount);
fprintf('predBest torqueZ range [%.3f %.3f] over %d angles\n', ...
    min(p.TorqueZ), max(p.TorqueZ), numel(p.TorqueZ));
fprintf('predBest momentArm range [%.4f %.4f] m\n', min(p.momentArm), max(p.momentArm));
fprintf('ctxCamBest fields: %s\n', strjoin(fieldnames(C.ctxCamBest)', ', '));
if isfield(C.ctxCamBest, 'routeSeed')
    fprintf('ctxCamBest.routeSeed =\n%s\n', num2str(C.ctxCamBest.routeSeed, 6));
end
% human target actually used (OpenSim_Vasti_Results column 5 = vas_med)
h = C.ctxCamBest;
if isfield(h, 'humanTorque')
    fprintf('humanTorque (vas_med) range [%.3f %.3f] N*m over %d angles\n', ...
        min(h.humanTorque), max(h.humanTorque), numel(h.humanTorque));
end
fprintf('CAM PROBE DONE\n');
