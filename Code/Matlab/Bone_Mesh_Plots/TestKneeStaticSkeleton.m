function tests = TestKneeStaticSkeleton
% Run with: runtests('TestKneeStaticSkeleton.m')
% Uses temporary synthetic meshes; no project data or optimization required.
tests = functiontests(localfunctions);
end

function setupOnce(tc)
plotterFile = which('AnimateKneeBoneMuscle');
tc.assertNotEmpty(plotterFile,'AnimateKneeBoneMuscle.m must be on the MATLAB path.');
plotterSource = fileread(plotterFile);
tc.assertTrue(contains(plotterSource,"addParameter(p,'SpineFile'") && ...
    contains(plotterSource,"addParameter(p,'FullSkeleton'"), ...
    sprintf(['MATLAB is using an older plotter: %s\n' ...
    'Replace THAT AnimateKneeBoneMuscle.m with the file in this package, ' ...
    'then run clear AnimateKneeBoneMuscle and rehash.'],plotterFile));
fprintf('Testing plotter: %s\n',plotterFile);
folder = tempname;
mkdir(folder);
tc.addTeardown(@() rmdir(folder,'s'));
tc.TestData.folder = folder;
oldVisible = get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
tc.addTeardown(@() set(groot,'defaultFigureVisible',oldVisible));
end

function testStaticBonesPreserveNativeXYZ(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,true);
fig = AnimateKneeBoneMuscle(args{:},'Static',true);
tc.addTeardown(@() close(fig));
verifyEqual(tc,points(fig,'PelvisMesh'),F.cloud*F.M,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'SacrumMesh'),F.cloud*F.M,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'SpineMesh'), ...
    (F.cloud+[-0.1007,0.0815,0])*F.M,'AbsTol',1e-12);
ankle = [0,-0.43,0];
subtalar = [-0.04877,-0.04195,0.00792];
mtp = [0.1788,-0.002,0.00108];
foot = [F.cloud+ankle;F.cloud+ankle+subtalar; ...
    F.cloud+ankle+subtalar+mtp];
verifyEqual(tc,points(fig,'FootMesh'), ...
    (transform(F.T(:,:,F.pos),foot)+F.hip)*F.M,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'TibiaMesh'), ...
    (transform(F.T(:,:,F.pos),F.cloud)+F.hip)*F.M,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,F.pos),'AbsTol',1e-12);
verifyEqual(tc,points(fig,'HumanRoute1'), ...
    human(F,F.pos),'AbsTol',1e-12);
h = findall(fig,'Tag','RobotRoute');
verifyTrue(tc,contains(h.Parent.Title.String,'Knee angle = 0.0 deg'));
verifyEqual(tc,getappdata(fig,'KneeBoneDisplayMap'),F.M,'AbsTol',1e-12);
verifyEqual(tc,h.Parent.XLabel.String,'X, m');
verifyEqual(tc,h.Parent.YLabel.String,'Y, m');
verifyEqual(tc,h.Parent.ZLabel.String,'Z, m');
end

function testExplicitMappingIsAppliedOnlyOnce(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,false);
M = [0,0,-1;0,1,0;-1,0,0];
fig = AnimateKneeBoneMuscle(args{:},'Static',true, ...
    'StaticFullSkeleton',false,'DisplayAxisMap',M);
tc.addTeardown(@() close(fig));
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,F.pos)*M,'AbsTol',1e-12);
end

function testAlreadyMappedRotationIsNotMappedAgain(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,false);
M = [0,0,-1;0,1,0;-1,0,0];
fig = AnimateKneeBoneMuscle(args{:},'Static',true, ...
    'StaticFullSkeleton',false,'DisplayRotation',M);
tc.addTeardown(@() close(fig));
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,F.pos)*M,'AbsTol',1e-12);
end

function testCameraOrbitDoesNotSwapPointCoordinates(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,false);
fig = AnimateKneeBoneMuscle(args{:},'Static',true, ...
    'FullSkeleton',false,'DisplayAxisMap',eye(3),'CameraOrbitDeg',90);
tc.addTeardown(@() close(fig));
% Literal native XYZ, independent of any proposed swap matrix.
verifyEqual(tc,points(fig,'FemurMesh'),F.cloud+F.hip,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,F.pos),'AbsTol',1e-12);
verifyEqual(tc,getappdata(fig,'KneeBoneDisplayMap'),eye(3),'AbsTol',1e-12);
end

function testFullSkeletonAnimationMovesFootAndFitsWholeSweep(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,true);
fig = AnimateKneeBoneMuscle(args{:},'FullSkeleton',true,'FrameIndices',1:100);
tc.addTeardown(@() close(fig));
verifyTrue(tc,getappdata(fig,'KneeBoneFullSkeleton'));
verifyEqual(tc,points(fig,'PelvisMesh'),F.cloud,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'SacrumMesh'),F.cloud,'AbsTol',1e-12);
verifyEqual(tc,points(fig,'SpineMesh'),F.cloud+[-0.1007,0.0815,0],'AbsTol',1e-12);
verifyEqual(tc,points(fig,'FootMesh'),foot(F,100),'AbsTol',1e-12);
verifyGreaterThan(tc,norm(foot(F,100)-foot(F,F.pos),'fro'),0.01);
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,100),'AbsTol',1e-12);
h = findall(fig,'Tag','RobotRoute');
ax = h.Parent;
lo = [ax.XLim(1),ax.YLim(1),ax.ZLim(1)];
hi = [ax.XLim(2),ax.YLim(2),ax.ZLim(2)];
for i = 1:100
    P = [foot(F,i);robot(F,i);human(F,i); ...
        transform(F.T(:,:,i),F.cloud)+F.hip];
    verifyTrue(tc,all(all(P > lo & P < hi)));
end
end

function testStaticAndFullAnimationMatchAtZero(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,true);
figStatic = AnimateKneeBoneMuscle(args{:},'Static',true);
tc.addTeardown(@() close(figStatic));
figMoving = AnimateKneeBoneMuscle(args{:},'FullSkeleton',true, ...
    'FrameIndices',[1,100,F.pos]);
tc.addTeardown(@() close(figMoving));
for tag = {'FemurMesh','TibiaMesh','PelvisMesh','SacrumMesh','SpineMesh', ...
        'FootMesh','RobotRoute','HumanRoute1'}
    verifyEqual(tc,points(figStatic,tag{1}),points(figMoving,tag{1}),'AbsTol',1e-12);
end
end

function testStaticLimitsContainFullSkeletonWithPadding(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,true);
% Deliberately too-small explicit limits must not clip the full skeleton.
fig = AnimateKneeBoneMuscle(args{:},'Static',true, ...
    'XLim',[-0.01,0.01],'YLim',[-0.01,0.01],'ZLim',[-0.01,0.01], ...
    'StaticPadding',[0.15,0.20,0.25]);
tc.addTeardown(@() close(fig));
h = findall(fig,'Tag','RobotRoute');
ax = h.Parent;
lo = [ax.XLim(1),ax.YLim(1),ax.ZLim(1)];
hi = [ax.XLim(2),ax.YLim(2),ax.ZLim(2)];
for tag = {'FemurMesh','TibiaMesh','PelvisMesh','SacrumMesh', ...
        'SpineMesh','FootMesh','RobotRoute','HumanRoute1'}
    P = points(fig,tag{1});
    verifyTrue(tc,all(all(P-lo >= [0.15,0.20,0.25]-1e-12)));
    verifyTrue(tc,all(all(hi-P >= [0.15,0.20,0.25]-1e-12)));
end
end

function testAnimationKeepsPerPoseRoutesAndDoesNotLoadExtraBones(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,false);
% Missing static files are intentional: animation must never open them.
fig = AnimateKneeBoneMuscle(args{:},'FullSkeleton',false,'FrameIndices',1:100);
tc.addTeardown(@() close(fig));
verifyEmpty(tc,findall(fig,'Tag','FootMesh'));
verifyEmpty(tc,findall(fig,'Tag','SpineMesh'));
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,100),'AbsTol',1e-12);
expected = robot(F,100);
verifyEqual(tc,points(fig,'RobotEnd'),expected(end,:),'AbsTol',1e-12);
verifyGreaterThan(tc,norm(robot(F,1)-robot(F,100),'fro'),0.1);
h = findall(fig,'Tag','RobotRoute');
verifyEqual(tc,h.Color,[205,52,181]/255,'AbsTol',1e-12);
hEnd = findall(fig,'Tag','RobotEnd');
verifyEqual(tc,hEnd.Color,[157,2,215]/255,'AbsTol',1e-12);
hHuman = findall(fig,'Tag','HumanRoute1');
verifyEqual(tc,hHuman.Color,[250,135,117]/255,'AbsTol',1e-12);
ax = h.Parent;
allP = (F.cloud+F.hip)*F.M;
for i = 1:100
    allP = [allP;robot(F,i);human(F,i); ... %#ok<AGROW>
        (transform(F.T(:,:,i),F.cloud)+F.hip)*F.M];
end
lo = min(allP,[],1);
hi = max(allP,[],1);
pad = max(0.005,0.05*(hi-lo));
verifyEqual(tc,ax.XLim,[lo(1)-pad(1),hi(1)+pad(1)+0.10],'AbsTol',1e-12);
verifyEqual(tc,ax.YLim,[lo(2)-pad(2),hi(2)+pad(2)],'AbsTol',1e-12);
verifyEqual(tc,ax.ZLim,[lo(3)-pad(3)-0.05,hi(3)+pad(3)+0.05],'AbsTol',1e-12);
end

function testKneeOnlyStaticDoesNotRequireExtraFiles(tc)
F = fixture(tc.TestData.folder);
args = plotArgs(F,false);
fig = AnimateKneeBoneMuscle(args{:},'Static',true,'FullSkeleton',false);
tc.addTeardown(@() close(fig));
verifyEmpty(tc,findall(fig,'Tag','PelvisMesh'));
verifyEqual(tc,points(fig,'RobotRoute'),robot(F,F.pos),'AbsTol',1e-12);
end

function F = fixture(folder)
F.M = eye(3);
F.hip = [-0.0707,-0.0661,0.0835];
F.cloud = [0,0,0;0.03,-0.10,0.02;-0.01,0.12,-0.04];
F.file = fullfile(folder,'mesh.csv');
writematrix(F.cloud,F.file);
F.missingFile = fullfile(folder,'intentionally_missing.csv');
F.phi = linspace(-2*pi/3,pi/18,100).';
F.pos = 93;
F.phi(F.pos) = 0;
F.T = repmat(eye(4),1,1,100);
F.TPam = F.T;
F.L = zeros(3,3,100);
F.H = struct('Location',[0.01,0.05,0.02;0.04,-0.08,0.01],'Cross',2);
for i = 1:100
    a = F.phi(i);
    R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
    F.T(:,:,i) = [R,[-0.0045;-0.3958;0];0,0,0,1];
    F.TPam(:,:,i) = [R,[0.02;-0.42;0.003];0,0,0,1];
    F.L(:,:,i) = [0.04,0.14,0.001;0.01+i/10000,0.02,0.003; ...
        0.034,-0.21,0.001];
end
end

function args = plotArgs(F,fullFiles)
extraFile = F.missingFile;
if fullFiles
    extraFile = F.file;
end
args = {F.T,F.T,F.phi,F.pos,F.L(1,:,F.pos),F.L(end,:,F.pos),{F.H}, ...
    'FemurFile',F.file,'TibiaFile',F.file, ...
    'SpineFile',extraFile,'SacrumFile',extraFile,'PelvisFile',extraFile, ...
    'TalusFile',extraFile,'CalcaneusFile',extraFile,'ToesFile',extraFile, ...
    'Location',F.L,'CrossPoint',2,'T_Pam',F.TPam, ...
    'Hip',F.hip,'DisplayRotation',eye(3), ...
    'XLim',[],'YLim',[],'ZLim',[],'Loop',false,'PauseTime',0};
end

function P = robot(F,i)
P = F.L(:,:,i);
P(2:end,:) = transform(F.TPam(:,:,i),P(2:end,:));
P = (P+F.hip)*F.M;
end

function P = human(F,i)
P = F.H.Location;
P(2:end,:) = transform(F.T(:,:,i),P(2:end,:));
P = (P+F.hip)*F.M;
end

function P = foot(F,i)
ankle = [0,-0.43,0];
subtalar = [-0.04877,-0.04195,0.00792];
mtp = [0.1788,-0.002,0.00108];
P = [F.cloud+ankle;F.cloud+ankle+subtalar;F.cloud+ankle+subtalar+mtp];
P = (transform(F.T(:,:,i),P)+F.hip)*F.M;
end

function P = transform(T,P)
P = (T*[P,ones(size(P,1),1)].').';
P = P(:,1:3);
end

function P = points(fig,tag)
h = findall(fig,'Tag',tag);
P = [h.XData(:),h.YData(:),h.ZData(:)];
end
