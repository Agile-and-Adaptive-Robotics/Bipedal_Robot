%dbg_geo_20261026.m -- print the cam ctx geo constants + active-release
% schedule needed to design the schematic bone outlines and the 4-pose set.
root = 'D:/GitHub/Bipedal_Robot';
S = load(fullfile(root, 'Code/Matlab/Mesh_Optimization/Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
ctx = S.ctxCamBest; ri = S.predBest.routeInfo; phiD = ctx.phiD(:);
geo = ctx.geo; N = numel(phiD);
fprintf('phiD: %d poses [%.1f %.1f]\n', N, min(phiD), max(phiD));
fprintf('femurCylCenter: %s  R=%.4f\n', mat2str(geo.femurCylCenter(1:2)), geo.femurCylClearRadius);
fprintf('femurLineX=%.4f  femurLineY=%s\n', geo.femurLineX, mat2str(geo.femurLineY));
fprintf('femurOffsetBoundary:\n'); disp(geo.femurOffsetBoundary(:,1:2));
fprintf('condyleClipX=%.4f clipY=%s\n', geo.femurCondyleClipX, mat2str(geo.femurCondyleClipY));
fprintf('tibiaLowerCenter: %s  R=%.4f\n', mat2str(geo.tibiaLowerCenter(1:2)), geo.tibiaLowerClearRadius);
fprintf('tibiaUpperCenter: %s  R=%.4f\n', mat2str(geo.tibiaUpperCenter(1:2)), geo.tibiaUpperClearRadius);
transitionIdx = find(any(ri.active(:,1:end-1) & ~ri.active(:,2:end), 1));
for t = transitionIdx
    row = find(ri.active(:,t) & ~ri.active(:,t+1));
    fprintf('RELEASE: pose %d (%.2f deg) row p%d off\n', t, phiD(t), row);
end
rowsE = find(isfinite(ri.eliminatedAngleD));
for kk = 1:numel(rowsE)
    fprintf('eliminatedAngleD p%d @ %.2f deg\n', rowsE(kk), ri.eliminatedAngleD(rowsE(kk)));
end
% ICR + tibia axis sample at full flexion / extension for shaft drawing
for ii = [1 N]
    Tt1 = ctx.T_t1_ICR(:,:,ii); Tp = ctx.T_Pam(:,:,ii);
    fprintf('pose %d (%.1f deg): t1 origin in femur = ICR*pam...\n', ii, phiD(ii));
    oICR = RowVecTrans(Tt1, [0 0 0]); oF = RowVecTrans(Tp, oICR);
    xICR = RowVecTrans(Tt1, [0.1 0 0]); xF = RowVecTrans(Tp, xICR);
    fprintf('  t1 origin -> femur %s ; t1 +x(0.1m) -> %s\n', mat2str(oF(1:2)), mat2str(xF(1:2)));
end
disp('DONE_dbg_geo');
