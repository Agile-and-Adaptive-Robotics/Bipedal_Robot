% Read-only inspection of the flexor result of record for R1 prose facts.
% Nothing is written; records untouched.
R = 'D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Results/Bifemsh_20mm_Result_20260921_2037.mat';
S = load(R);
disp('=== top-level fields ==='); disp(fieldnames(S));
if isfield(S,'predBest')
    p = S.predBest;
    disp('=== predBest fields (route-relevant) ===');
    disp(intersect(fieldnames(p), {'p1';'p2';'p1B';'p2B';'rest';'tendon';'BPAcount';'bendMeasure';'bendMeasure2';'routeInfo';'routeInfo2';'Location';'Location2';'TorqueZ'}, 'stable'));
    fprintf('BPAcount = %d\n', p.BPAcount);
    fprintf('p1  = [%g %g %g]  p2 = [%g %g %g]\n', p.p1, p.p2);
    if isfield(p,'p1B')
        fprintf('p1B = [%g %g %g]  p2B = [%g %g %g]\n', p.p1B, p.p2B);
        fprintf('p2B(3) = %g vs -p2(3) = %g (mirrored distal z)\n', p.p2B(3), -p.p2(3));
    end
    if isfield(p,'routeInfo') && isstruct(p.routeInfo)
        ri = p.routeInfo;
        if isfield(ri,'active')
            fprintf('BPA1 wrap: nActive=%d of N=%d, activeAtExtension=%d\n', sum(ri.active), numel(ri.active), ri.activeAtExtension);
        end
        if isfield(ri,'releaseAngleD'); fprintf('BPA1 releaseAngleD = %g\n', ri.releaseAngleD); end
    end
    if isfield(p,'routeInfo2') && isstruct(p.routeInfo2)
        ri2 = p.routeInfo2;
        if isfield(ri2,'active')
            fprintf('BPA2 wrap: nActive=%d of N=%d, activeAtExtension=%d\n', sum(ri2.active), numel(ri2.active), ri2.activeAtExtension);
        end
        if isfield(ri2,'releaseAngleD'); fprintf('BPA2 releaseAngleD = %g\n', ri2.releaseAngleD); end
    end
    if isfield(p,'bendMeasure'); fprintf('bendMeasure1: min %g max %g (m)\n', min(p.bendMeasure), max(p.bendMeasure)); end
    if isfield(p,'bendMeasure2'); fprintf('bendMeasure2: min %g max %g (m)\n', min(p.bendMeasure2), max(p.bendMeasure2)); end
    if isfield(p,'TorqueZ'); fprintf('TorqueZ range %g .. %g N*m\n', min(p.TorqueZ), max(p.TorqueZ)); end
end
if isfield(S,'xBest'); fprintf('xBest = [%s]\n', num2str(S.xBest, '%g ')); end
if isfield(S,'ctx') && isstruct(S.ctx)
    c = S.ctx;
    if isfield(c,'BPAcount'); fprintf('ctx.BPAcount = %d\n', c.BPAcount); end
    if isfield(c,'Xi3'); fprintf('ctx.Xi3 = %g\n', c.Xi3); end
    if isfield(c,'Xi0'); fprintf('ctx.Xi0 = %g\n', c.Xi0); end
    if isfield(c,'requiredTorqueMargin'); fprintf('ctx.requiredTorqueMargin = %g\n', c.requiredTorqueMargin); end
end
disp('DONE');
