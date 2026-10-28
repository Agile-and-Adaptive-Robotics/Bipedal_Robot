root='D:/Github/Bipedal_Robot';
addpath(fullfile(root,'Code/Matlab/Functions'));addpath(fullfile(root,'Code/Matlab/Functions/ModernRobotics'));addpath(fullfile(root,'Code/Matlab/Robot_Data'));addpath(fullfile(root,'Code/Matlab/Mesh_Optimization'));addpath(fullfile(root,'Code/Matlab/Mesh_Optimization/Results'));
cd(fullfile(root,'Code/Matlab/Mesh_Optimization'));warning('off','all');
%% Joint kinematics (copied from Knee_Flexor_Data_20mm.m)
positions = 100;

knee_angle_x = [-2.0944; -1.74533; -1.39626; -1.0472; -0.698132; -0.349066; -0.174533;  0.197344;  0.337395;  0.490178;   1.52146;   2.0944];
knee_x =       [-0.0032;  0.00179;  0.00411;  0.0041;   0.00212;    -0.001;   -0.0031; -0.005227; -0.005435; -0.005574; -0.005435; -0.00525];
fcn1 = fit(knee_angle_x,knee_x,'cubicspline');
knee_angle_y = [-2.0944; -1.22173; -0.523599; -0.349066; -0.174533;  0.159149; 2.0944];
knee_y =       [-0.4226;  -0.4082;    -0.399;   -0.3976;   -0.3966; -0.395264; -0.396];
fcn2 = fit(knee_angle_y,knee_y,'cubicspline');
knee_angle = [0.17; 0.09; 0.03; 0.00; -0.09; -0.17; -0.26; -0.52; -0.79; -1.05; -1.31; -1.57; -1.83; -2.09; -2.36; -2.62];
knee_x_Pam =     ([23.30	22.22	21.55	21.09	19.91	18.70	17.48	13.82	10.44	7.60	5.52	4.35	4.16	5.01	7.04	10.47]')/1000;
fcn3 = fit(knee_angle,knee_x_Pam,'cubicspline');
knee_y_Pam =     ([-416.65	-417.03	-417.19	-417.28	-417.41	-417.41	-417.30	-416.28	-414.36	-411.72	-408.62	-405.32	-402.08	-399.16	-396.85	-395.66]')/1000;
fcn4 = fit(knee_angle,knee_y_Pam,'cubicspline');

t1_ICR_x = ([29.66	28.54	27.86	27.40	26.23	25.03	23.81	20.03	16.17	12.34	8.67	5.24	2.04	-1.01	-4.1	-7.58]')/1000;
fcn13 = fit(knee_angle,t1_ICR_x,'cubicspline');
t1_ICR_y = ([25.97	25.74	25.61	25.53	25.35	25.19	25.03	24.57	24.04	23.39	22.66	21.93	22.32	21.99	21.2	22.33]')/1000;
fcn14 = fit(knee_angle,t1_ICR_y,'cubicspline');

kneeMin = -2.0943951;
kneeMax = 0.17453293;
phi = linspace(kneeMin, kneeMax, positions);
[val, pos] = min(abs(phi));
phi(pos) = 0;

T_Pam = zeros(4, 4, positions);
T_t1_ICR = zeros(4, 4, positions);
for i = 1:positions
    hipToKnee_Pam = [fcn3(phi(i)), fcn4(phi(i)), 0];
    R_Pam = [cos(phi(i)), -sin(phi(i)), 0;
             sin(phi(i)),  cos(phi(i)), 0;
             0, 0, 1];
    T_Pam(:, :, i) = RpToTrans(R_Pam, hipToKnee_Pam');
    t1toICR = [fcn13(phi(i)), fcn14(phi(i)), 0];
    T_t1_ICR(:, :, i) = RpToTrans(eye(3), t1toICR');
end
phiD = phi*180/pi;

%% Design of record: flexor pair route (copied from Knee_Flexor_Data_20mm.m)
Name = 'Bicep Femoris (Short Head)';
CrossPoint = 2;

resultData = load('Bifemsh_20mm_Result.mat', ...
    'xBest', 'routeCtx', 'Xi3', 'XiUsed');
xBest = resultData.xBest;
routeCtx = resultData.routeCtx;
routeCtx.geo = buildGeoExclusion();

p1 = xBest(1:3);
p2 = xBest(4:6);

[Location, bendMeasure] = buildKneeFlexorRoute20mm(p1, p2, xBest(8), routeCtx);
[p1B, p2B] = flexorBpa2Endpoints20mm(p1, p2);
[LocationB, bendMeasureB] = ...
    buildKneeFlexorRoute20mm(p1B, p2B, xBest(8), routeCtx);
LocationPair = {Location; LocationB};
bendMeasurePair = {bendMeasure; bendMeasureB};

Dia = 20;
rest = xBest(7);
tendon = xBest(8) - 0.02;
KMAX = 0.235;                      % maximum contraction fraction of record
kmax = (1 - KMAX)*rest;
fitting = 0.021;
pres3 = 620;

XiUsed = resultData.XiUsed;
Xi0 = XiUsed(1); Xi1 = XiUsed(2); Xi2 = XiUsed(3); Xi3 = XiUsed(4);
fprintf('Design of record: rest = %.4f m, tendon = %.4f m, Xi = [%.6g %.6g %.6g %.6g]\n', ...
    rest, tendon, Xi0, Xi1, Xi2, Xi3)

wraps = 6;
BPAcount = 2;

Bifemsh_Pam3 = MonoPam_mult(Name, LocationPair, CrossPoint, Dia, T_Pam, ...
    rest, kmax, tendon, fitting, pres3, Xi0, Xi1, Xi2, Xi3, wraps, phiD, ...
    BPAcount, bendMeasurePair);
TqCorr  = Bifemsh_Pam3.Torque_p(:,3);   % corrected (chi0-chi3) at 620 kPa
TqRigid = Bifemsh_Pam3.Torque(:,3);     % original rigid-body at 620 kPa
fprintf('Pair torque at 620 kPa (N*m): rigid [min max] = [%.2f %.2f], corrected [min max] = [%.2f %.2f]\n', ...
    min(TqRigid), max(TqRigid), min(TqCorr), max(TqCorr))
fprintf('NaN count: rigid %d, corrected %d\n', sum(isnan(TqRigid)), sum(isnan(TqCorr)))

K_deg=[-61.5 -83 5 -5 -52 -45 -42 -114 -119 -62 -78 -87 -105 -22.5 -14];
Presh=[517 628 226.5 242.5 406.7 410.3 412 634 633 400 489 635 638 275.3 273.7];
predAP=zeros(size(Presh));
for u=1:numel(Presh)
 v=MonoPam_mult(Name,LocationPair,CrossPoint,Dia,T_Pam,rest,kmax,tendon,fitting,Presh(u),Xi0,Xi1,Xi2,Xi3,wraps,phiD,BPAcount,bendMeasurePair);
 predAP(u)=interp1(phiD,v.Torque_p(:,3),K_deg(u),'pchip');
end
figfile=fullfile(root,'Documentation/Reports and Papers/Knee_Torque_Test/Figures/Components/narrative_revision_20261028/Redesigned_flexor_validation.fig');
f=openfig(figfile,'invisible');ax=findall(f,'Type','axes');ll=findall(ax,'Type','line');
err=Inf;
for j=1:numel(ll)
 if numel(ll(j).YData)==numel(TqCorr)
  err=min(err,max(abs(ll(j).YData(:)-TqCorr(:))));
 end
end
fprintf('FLEXOR current 620 curve verification max error %.12g N m\n',err);assert(err<1e-8,'Saved plot/model mismatch');
hold(ax,'on');h=scatter(ax,K_deg,predAP,55,'o','filled','MarkerFaceColor',[1 .843 0],'MarkerEdgeColor',[.45 .35 0],'DisplayName','Prediction at measured pressure');
lg=findall(f,'Type','legend');lg.AutoUpdate='on';legend(ax,'show');lg.AutoUpdate='off';
lg.Units='normalized';pos=lg.Position;pos(1)=pos(1)+.035;pos(2)=pos(2)+.065;lg.Position=pos;lg.FontSize=10;
ylabel(ax,'Torque, N m','Interpreter','none');
outdir=fullfile(root,'Documentation/Reports and Papers/Knee_Torque_Test/Figures');
exportgraphics(f,fullfile(outdir,'Redesigned_flexor_pressure_publication.pdf'),'ContentType','vector');
savefig(f,fullfile(root,'tmp/paper_bridge/figscripts/Redesigned_flexor_pressure_publication.fig'));
out=struct('angle',K_deg,'pressure_kPa',Presh,'prediction',predAP,'parameters',XiUsed,'curve620_max_error',err,'classification','model evaluated at measured pressure; not measured-length hybrid');
fid=fopen(fullfile(root,'tmp/paper_bridge/figscripts/flexor_pressure_publication.json'),'w');fprintf(fid,'%s',jsonencode(out));fclose(fid);close(f);


