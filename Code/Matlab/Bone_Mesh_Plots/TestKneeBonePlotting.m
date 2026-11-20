function tests = TestKneeBonePlotting
% Synthetic graphics checks: no real meshes or BPA classes are needed.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
folder = tempname;
mkdir(folder);
testCase.addTeardown(@() rmdir(folder,'s'));
testCase.TestData.folder = folder;
testCase.TestData.package = fileparts(mfilename('fullpath'));
oldVisibility = get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
testCase.addTeardown(@() set(groot,'defaultFigureVisible',oldVisibility));
end

function setup(testCase)
F = fixture(testCase.TestData.folder);
testCase.TestData.F = F;
end

function testStaticSelectsZeroAndPreservesNativeFrames(testCase)
F = testCase.TestData.F;
args = plotArgs(F);
fig = AnimateKneeBoneMuscle(args{:},'Static',true);
testCase.addTeardown(@() closeFigure(fig));
expected = robotExpected(F,F.pos);
verifyEqual(testCase,pathData(fig,'RobotRoute'),expected,'AbsTol',1e-12);
verifyEqual(testCase,pathData(fig,'RobotOrigin'),expected(1,:),'AbsTol',1e-12);
verifyEqual(testCase,pathData(fig,'RobotEnd'),expected(end,:),'AbsTol',1e-12);
h = findobj(fig,'Tag','RobotRoute');
verifyTrue(testCase,contains(h.Parent.Title.String,'Knee angle = 0.0 deg'));
% The actual t1 endpoint is converted once into ICR, then once into femur.
endpoint = homogeneous(F.TPam(:,:,F.pos)*F.TICR(:,:,F.pos),F.pEnd);
verifyEqual(testCase,expected(end,:),displayPoints(F,endpoint),'AbsTol',1e-12);
end

function testAnimationUsesChangingWrapAndLastRow(testCase)
F = testCase.TestData.F;
args = plotArgs(F);
fig = AnimateKneeBoneMuscle(args{:},'FrameIndices',[1,2,3]);
testCase.addTeardown(@() closeFigure(fig));
expected = robotExpected(F,3);
verifyEqual(testCase,pathData(fig,'RobotRoute'),expected,'AbsTol',1e-12);
verifyEqual(testCase,pathData(fig,'RobotEnd'),expected(9,:),'AbsTol',1e-12);
atZero = robotExpected(F,2);
verifyGreaterThan(testCase,norm(expected(7,:)-atZero(7,:)),1e-3);
% Eliminated row 8 repeats row 9 at this frame; it is not a stale contact.
verifyEqual(testCase,expected(8,:),expected(9,:),'AbsTol',1e-12);
end

function testStaticAndAnimationMatchAtZero(testCase)
F = testCase.TestData.F;
args = plotArgs(F);
figStatic = AnimateKneeBoneMuscle(args{:},'Static',true);
testCase.addTeardown(@() closeFigure(figStatic));
figMoving = AnimateKneeBoneMuscle(args{:},'FrameIndices',[1,3,F.pos]);
testCase.addTeardown(@() closeFigure(figMoving));
for tag = {'RobotRoute','RobotOrigin','RobotEnd','TibiaMesh','HumanRoute1','HumanRoute2'}
    verifyEqual(testCase,pathData(figMoving,tag{1}),pathData(figStatic,tag{1}), ...
        'AbsTol',1e-12);
end
end

function testHumanPathsUseHumanNotRobotTransforms(testCase)
F = testCase.TestData.F;
args = plotArgs(F);
fig = AnimateKneeBoneMuscle(args{:},'FrameIndices',[1,3]);
testCase.addTeardown(@() closeFigure(fig));
for h = 1:numel(F.humans)
    H = F.humans{h};
    P = H.Location(:,:,3);
    P(H.Cross:end,:) = homogeneous(F.T(:,:,3),P(H.Cross:end,:));
    verifyEqual(testCase,pathData(fig,sprintf('HumanRoute%d',h)), ...
        displayPoints(F,P),'AbsTol',1e-12);
end
verifyEqual(testCase,pathData(fig,'TibiaMesh'), ...
    displayPoints(F,homogeneous(F.T(:,:,3),F.cloud)),'AbsTol',1e-12);
end

function testDisplayedLengthIsRouteNotChord(testCase)
F = testCase.TestData.F;
args = plotArgs(F);
fig = AnimateKneeBoneMuscle(args{:},'Static',true);
testCase.addTeardown(@() closeFigure(fig));
P = robotExpected(F,F.pos);
lengthRoute = sum(sqrt(sum(diff(P).^2,2)));
verifyGreaterThan(testCase,lengthRoute,norm(P(end,:)-P(1,:))+1e-3);
h = findobj(fig,'Tag','RobotRoute');
verifyTrue(testCase,contains(h.Parent.Title.String,sprintf('%.3f m',lengthRoute)));
end

function testBoneMarkersAndGlobalDefaults(testCase)
F = testCase.TestData.F;
oldLine = get(groot,'defaultLineMarkerEdgeColor');
oldScatter = get(groot,'defaultScatterMarkerEdgeColor');
testCase.addTeardown(@() set(groot,'defaultLineMarkerEdgeColor',oldLine, ...
    'defaultScatterMarkerEdgeColor',oldScatter));
set(groot,'defaultLineMarkerEdgeColor','none','defaultScatterMarkerEdgeColor','none');
args = plotArgs(F);
fig = AnimateKneeBoneMuscle(args{:},'Static',true);
testCase.addTeardown(@() closeFigure(fig));
verifyEqual(testCase,get(groot,'defaultLineMarkerEdgeColor'),'none');
verifyEqual(testCase,get(groot,'defaultScatterMarkerEdgeColor'),'none');
for tag = {'FemurMesh','TibiaMesh','RobotRoute','RobotOrigin','RobotEnd'}
    h = findobj(fig,'Tag',tag{1});
    verifyEqual(testCase,h.MarkerEdgeColor,'none');
    verifyTrue(testCase,isnumeric(h.MarkerFaceColor));
    verifyEqual(testCase,h.Marker,'o');
end
end

function testAutoLimitsIncludeAllRouteFrames(testCase)
F = testCase.TestData.F;
args = plotArgs(F);
fig = AnimateKneeBoneMuscle(args{:},'FrameIndices',[1,2,3]);
testCase.addTeardown(@() closeFigure(fig));
h = findobj(fig,'Tag','RobotRoute');
lo = [h.Parent.XLim(1),h.Parent.YLim(1),h.Parent.ZLim(1)];
hi = [h.Parent.XLim(2),h.Parent.YLim(2),h.Parent.ZLim(2)];
for i = 1:3
    P = robotExpected(F,i);
    verifyTrue(testCase,all(all(P > lo & P < hi)));
end
end

function testLegacyTwoPointCallsStillWork(testCase)
F = testCase.TestData.F;
fig = AnimateKneeBoneMuscle(F.T,F.TICR,F.phi,F.pos,F.p1,F.pEnd,F.humans{1}, ...
    'FemurFile',F.femurFile,'TibiaFile',F.tibiaFile, ...
    'Hip',F.hip,'DisplayRotation',F.R,'Loop',false,'PauseTime',0, ...
    'FrameIndices',[1,3]);
testCase.addTeardown(@() closeFigure(fig));
v2 = homogeneous(F.TICR(:,:,F.pos),F.pEnd);
P = [F.p1;homogeneous(F.T(:,:,3),v2)];
verifyEqual(testCase,pathData(fig,'RobotRoute'),displayPoints(F,P),'AbsTol',1e-12);
end

function testMuscleBonePlottingScriptSelectsStaticMode(testCase)
F = testCase.TestData.F;
bonePlotArgs = plotArgs(F); %#ok<NASGU>
run(fullfile(testCase.TestData.package,'MuscleBonePlotting.m'));
fig = gcf;
testCase.addTeardown(@() closeFigure(fig));
verifyEqual(testCase,pathData(fig,'RobotRoute'),robotExpected(F,F.pos),'AbsTol',1e-12);
end

function testRejectsMissingRobotTransforms(testCase)
F = testCase.TestData.F;
verifyError(testCase,@() AnimateKneeBoneMuscle( ...
    F.T,F.TICR,F.phi,F.pos,F.p1,F.pEnd,F.humans, ...
    'Location',F.Location,'CrossPoint',6,'Static',true), ...
    'AnimateKneeBoneMuscle:Transforms');
end

function testResultScriptRebuildsAndUsesSignedTarget(testCase)
source = fileread('Knee_Extensor_20mm.m');
verifyFalse(testCase,contains(source,'S.Location'));
verifyFalse(testCase,contains(source,'S.p1'));
verifyFalse(testCase,contains(source,'humanTorqueAbs'));
verifyTrue(testCase,contains(source,'buildDistalRingLocation20mm(p1,pEnd,tendon,ctx)'));
verifyTrue(testCase,contains(source,'ctx.humanAngleD,ctx.humanTorque,phiD'));
verifyTrue(testCase,contains(source,'p1 = xBest(1:3)'));
verifyTrue(testCase,contains(source,'pEnd = xBest(4:6)'));
verifyTrue(testCase,contains(source,'kmax = rest*(1-KMAX)'));
verifyTrue(testCase,contains(source,"'FrameIndices',frames"));
end

function F = fixture(folder)
F.phi = [-2*pi/3;0;pi/18];
F.pos = 2;
F.p1 = [0.03,0.148,0.001];
F.pEnd = [0.033,-0.217,0.001];
F.hip = [-0.0707,-0.0661,0.0835];
F.R = [1,0,0;0,0,1;0,-1,0];
F.T = repmat(eye(4),1,1,3);
F.TPam = F.T;
F.TICR = F.T;
F.Location = zeros(9,3,3);
human1 = struct('Location',zeros(4,3,3),'Cross',3);
human2 = struct('Location',zeros(3,3,3),'Cross',2);
for i = 1:3
    a = F.phi(i);
    R = [cos(a),-sin(a),0;sin(a),cos(a),0;0,0,1];
    F.T(:,:,i) = [R,[-0.005;-0.4;0];0,0,0,1];
    F.TPam(:,:,i) = [R,[0.021;-0.417;0.006*i];0,0,0,1];
    F.TICR(1:3,4,i) = [-0.027-0.001*i;-0.025;0.002*i];
    P = [F.p1;0.084,-0.275,0.001;0.066,-0.427,0.001; ...
        0.052,-0.447,0.001;0.038,-0.453,0.001; ...
        0.070,0.034,0.001;0.074+0.002*i,0.019,0.001+0.003*i; ...
        0.072,-0.012,0.001;F.pEnd];
    if i == 3
        P(8,:) = P(9,:);
    end
    P(6:end,:) = homogeneous(F.TICR(:,:,i),P(6:end,:));
    F.Location(:,:,i) = P;
    human1.Location(:,:,i) = [0.01,-0.2,0.02;0.03,-0.3,0.01; ...
        0.04,-0.01,0.005;0.03+0.001*i,-0.06,0.002];
    human2.Location(:,:,i) = [0.02,-0.22,0.01;0.05,-0.02,0.01;0.04,-0.07,0.01];
end
F.humans = {human1,human2};
F.cloud = [0,0,0;0.01,-0.15,0.01;-0.01,-0.4,-0.01];
F.femurFile = fullfile(folder,'femur.csv');
F.tibiaFile = fullfile(folder,'tibia.csv');
writematrix(F.cloud,F.femurFile);
writematrix(F.cloud,F.tibiaFile);
end

function args = plotArgs(F)
% These older checks exercise native-frame routing, not the new display map
% or extra static bones. Keep them independent of real skeleton mesh files.
args = {F.T,F.TICR,F.phi,F.pos,F.p1,F.pEnd,F.humans, ...
    'FemurFile',F.femurFile,'TibiaFile',F.tibiaFile, ...
    'Location',F.Location,'CrossPoint',6,'T_Pam',F.TPam, ...
    'HumanLabels',{'Human 1','Human 2'},'Hip',F.hip,'DisplayRotation',F.R, ...
    'DisplayAxisMap',eye(3),'StaticFullSkeleton',false, ...
    'XLim',[],'YLim',[],'ZLim',[],'Loop',false,'PauseTime',0};
end

function P = robotExpected(F,i)
P = F.Location(:,:,i);
P(6:end,:) = homogeneous(F.TPam(:,:,i),P(6:end,:));
P = displayPoints(F,P);
end

function points = homogeneous(T,points)
P = (T*[points,ones(size(points,1),1)].').';
points = P(:,1:3);
end

function points = displayPoints(F,points)
points = (points+F.hip)*F.R;
end

function P = pathData(fig,tag)
h = findobj(fig,'Tag',tag);
P = [h.XData(:),h.YData(:),h.ZData(:)];
end

function closeFigure(fig)
if isgraphics(fig)
    close(fig);
end
end
