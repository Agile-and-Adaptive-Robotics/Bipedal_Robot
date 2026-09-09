%Probe: why is 40cm-tendon (kf(4)) returning NaN now?
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');
for tm = {'2trans','1trans'}
    [f, bpa] = minimizeFlxPin2brk(0, Inf, Inf, 4, true, tm{1});
    d = bpa(4);
    fprintf('== %s ==\n', tm{1});
    fprintf('  f = [%g %g %g]\n', f(1), f(2), f(3));
    fprintf('  NaN counts: L_p %d/%d | Lmt_p %d/%d | strain_p %d/%d | M_p %d/%d\n', ...
        sum(isnan(d.L_p(:))), numel(d.L_p), sum(isnan(d.Lmt_p)), numel(d.Lmt_p), ...
        sum(isnan(d.strain_p)), numel(d.strain_p), sum(isnan(d.M_p(:))), numel(d.M_p));
    fprintf('  rest=%.4f ten=%.4f Kmax=%.4f | Ak(1)=%.2f min|Ak|=%.2f | CP=%d size(Loc)=[%d %d %d]\n', ...
        d.rest, d.ten, d.Kmax, d.Ak(1), min(abs(d.Ak)), d.CP, size(d.Loc,1), size(d.Loc,2), size(d.Loc,3));
    fprintf('  strain (kf template) range: [%.3f %.3f]\n', min(d.strain), max(d.strain));
end
