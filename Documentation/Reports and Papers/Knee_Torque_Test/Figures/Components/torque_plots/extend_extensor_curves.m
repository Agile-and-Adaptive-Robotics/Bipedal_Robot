function extend_extensor_curves()
% Continue the original Figure 9 models through the measured zero in panel D.
here = fileparts(mfilename('fullpath'));
papers = fileparts(fileparts(fileparts(fileparts(here))));
repo = fileparts(fileparts(papers));
data = fullfile(repo,'Testing_Data','2022_02_Festo');
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(data,fullfile(repo,'Code','Matlab','Robot_Data'), ...
    fullfile(repo,'Code','Matlab','Functions'), ...
    fullfile(repo,'Code','Matlab','Functions','ModernRobotics'));
source = load(fullfile(data,'ExtPinBPASet.mat'),'ke');
published = load(fullfile(papers,'Dissertation','Figures','30-torque','Components', ...
    'extensor_figures','extensor_data.mat'),'tests','terms');
order = [published.tests.id];
trials = source.ke(order);
angles = trials(4).Aexp(trials(4).Mexp == 0 & trials(4).Aexp > 0);
assert(isscalar(angles) && angles > 35);
% Replay the entire saved grid and continue to the last measured angle.
extra = (35.1:0.1:angles)';
for i = 1:numel(trials)
    trial = trials(i);
    oldAngles = trial.Ak(:);
    assert(oldAngles(end) == 35);
    % Above the final contact release (27 degrees), the route is straight.
    % Replay the saved unwrapped geometry before extending the knee rotation.
    unwrapped = trial.Loc(:,:,end);
    assert(all(trial.Loc(:,:,oldAngles>27) == unwrapped,'all'));
    allAngles = [oldAngles;extra];
    T = repmat(trial.Tk(:,:,end),1,1,numel(allAngles));
    T(1,1,:) = cosd(allAngles); T(1,2,:) = -sind(allAngles);
    T(2,1,:) = sind(allAngles); T(2,2,:) = cosd(allAngles);
    assert(max(abs(T(:,:,1:numel(oldAngles))-trial.Tk),[],'all') < 1e-12);
    trial.Tk = cat(3,trial.Tk,T(:,:,numel(oldAngles)+1:end));
    trial.Loc = cat(3,trial.Loc,repmat(unwrapped,1,1,numel(extra)));
    trial.Ak = allAngles;
    rigid = MonoPamDataExplicit('Pinned extensor',trial.Loc,trial.CP, ...
        trial.dBPA,trial.Tk,trial.rest,trial.Kmax,trial.ten,trial.fitn,trial.P);
    trial.Lmt = rigid.MuscleLength;
    torque = rigid.Torque;
    trial.M = torque(:,3);
    trials(i) = trial;
end
p = published.terms;
[~,fitted] = minimizeExtX3(p(1),p(2),p(3),p(4),1:numel(trials),'2trans',trials);
for i = 1:numel(trials)
    saved = published.tests(i);
    n = numel(saved.angle);
    baseline = trials(i).M;
    corrected = fitted(i).M_p(:,3);
    assert(isequal(isnan(baseline(1:n)),isnan(saved.unoptimized(:))));
    assert(isequal(isnan(corrected(1:n)),isnan(saved.optimized(:))));
    assert(max(abs(baseline(1:n) - saved.unoptimized(:)),[],'omitnan') < 1e-9);
    assert(max(abs(corrected(1:n) - saved.optimized(:)),[],'omitnan') < 1e-9);
    extension(i) = struct('id',saved.id,'angle',extra,'unoptimized', ...
        baseline(n+1:end),'optimized',corrected(n+1:end)); %#ok<AGROW>
end
assert(extension(4).angle(end) == angles && extension(4).unoptimized(end) == 0 ...
    && extension(4).optimized(end) == 0);
save(fullfile(here,'extensor_extension.mat'),'extension','-v7');
end
