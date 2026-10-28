function checkpoint()
here=fileparts(mfilename('fullpath')); R=here;
while ~strcmpi(string(lastpart(R)),"Bipedal_Robot"), R=fileparts(R); end
addpath(genpath(fullfile(R,'Code','Matlab')));
fid=fopen(fullfile(here,'checkpoint.txt'),'w'); cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'MATLAB %s\n',version); fprintf('MATLAB %s\n',version);
% Prevent parfor from spawning a large pool during this bounded diagnostic.
ps=parallel.Settings; ps.Pool.AutoCreate=false;
C=load(fullfile(R,'Code','Matlab','Mesh_Optimization','Results','Vas_Pam_20mm_Result_cam_20261014_0416.mat'));
ctx=C.ctxCamBest; x=C.xFullBest(:).'; Xi=C.XiUsed(:).';
fprintf(fid,'Saved fields: %s\n',strjoin(fieldnames(C),', '));
fprintf(fid,'ctx fields: %s\n',strjoin(fieldnames(ctx),', '));
fprintf(fid,'nBPA stored pred=%d ctx=%d rest=%.16g tendon=%.16g\n',C.predBest.BPAcount,ctx.BPAcount,x(7),x(8));
p=predictKneeExt20mm(x,ctx); assert(p.ok,p.failReason);
d=abs(p.TorqueZ(:)-C.predBest.TorqueZ(:));
fprintf(fid,'BASELINE max torque replay difference %.16g Nm; torque %.12g..%.12g\n',max(d),min(p.TorqueZ),max(p.TorqueZ));
a=ctx.phiD(:); hs=interp1(ctx.humanAngleD(:),ctx.humanTorque(:),a,'linear','extrap');
cap=fullfile(R,'tmp','paper_bridge','capability_table_20261022'); T=readtable(fullfile(cap,'vas_med_surface_dense_20261027.csv'));
hnew=interp1(T.knee_angle_deg,T.tau_h_Nm,a,'linear','extrap');
H=readmatrix(fullfile(R,'Testing_Data','2022_02_Festo','OpenSim_Vasti_Results.txt'),'FileType','text','NumHeaderLines',7);
hdisk=interp1(H(:,2),H(:,5),a,'linear','extrap');
writetable(table(a,hs,1.05*hs,hnew,1.05*hnew,hdisk,hnew-hs,'VariableNames',{'angle_deg','saved_plain_Nm','saved_required_Nm','surface_plain_Nm','surface_required_Nm','disk_legacy_plain_Nm','surface_minus_saved_Nm'}),fullfile(here,'target_comparison.csv'));
fprintf(fid,'TARGET saved %.12g..%.12g surface %.12g..%.12g max diff %.12g; disk-saved %.12g\n',min(hs),max(hs),min(hnew),max(hnew),max(abs(hnew-hs)),max(abs(hdisk-hs)));
writetable(table(a,C.predBest.TorqueZ(:),p.TorqueZ(:),vecnorm(p.bpa.F_p,2,2),p.activeLength_effective,p.strain_f,'VariableNames',{'angle_deg','stored_torque_Nm','replayed_torque_Nm','total_force_N','effective_length_m','strain_f'}),fullfile(here,'saved_baseline.csv'));
% G1: previous alternative closure vs actual X3 class, identical n=1.
ctx1=ctx;ctx1.BPAcount=1; q=predictKneeExt20mm(x,ctx1); assert(q.ok,q.failReason);
L=q.Location; N=numel(a); mL=q.Lmt0;
comp=max(0,1-((x(7)-(mL-x(8)-2*ctx.fitting-Xi(1)))/x(7))/ctx.KMAX);
delta0=Xi(4).*q.bendMeasure(:).*comp.^2;
s0=x(7)-mL+Xi(1)+x(8)+2*ctx.fitting+delta0;
cb=bracketCb(L,ctx,Xi,N);
P=psiSolve(x(7),1,1,s0,zeros(N,1),cb,ctx.KMAX,ctx.Dia,ctx.targetPressure,x(8),ctx.wraps,@(r)maxBPAforce(r,'20'));
Fc=vecnorm(q.bpa.F_p,2,2); Falt=P.Fsingle;
err=abs(Falt-Fc)./max(1,abs(Fc)); [em,iw]=max(err);
fprintf(fid,'G1 FAIL max force rel=%.16g abs=%.16g N at %.12g deg; alternative %.16g canonical %.16g\n',em,max(abs(Falt-Fc)),a(iw),Falt(iw),Fc(iw));
fprintf(fid,'G1 length max=%.16g m; delta_L initial-vs-force-effective max=%.16g m\n',max(abs((x(7)-P.sCon)-q.activeLength_effective)),max(abs(delta0-x(7)*(q.strain_f-q.strain_p))));
% Use canonical moment arm only to isolate force contribution to torque.
talt=Falt.*cross(q.bpa.mA_p,q.bpa.F_p./max(Fc,eps),2); talt=talt(:,3);
fprintf(fid,'G1 force-only torque discrepancy using SAME canonical arm %.16g Nm\n',max(abs(talt-q.TorqueZ)));
writetable(table(a,Fc,Falt,err,q.activeLength_effective,x(7)-P.sCon,q.delta_L,delta0,q.TorqueZ,talt),fullfile(here,'gate1.csv'));
% Align bracket centroid, pre-force transform, spring, and total force exactly.
pa=L(1,:,92); pb=[9.48 -33.38 30.86]/1000; v=pa-pb;
rz=rotzLocal(atan2(v(2),v(1))); vv=rz'*v'; yy=atan2(vv(3),vv(1)); ry=[cos(yy) 0 sin(yy);0 1 0;-sin(yy) 0 cos(yy)]; tr=RpToTrans(rz*ry',pb');
Fu=zeros(N,3);
for k=1:N
 for j=2:size(L,1)
  if norm(L(j,:,k)-L(1,:,k))>1e-6
   pt=L(j,:,k);if j==ctx.CrossPoint,pt=RowVecTrans(ctx.T_Pam(:,:,k),pt);end
   Fu(k,:)=(pt-L(1,:,k))/norm(pt-L(1,:,k));break
  end
 end
end
pre=festo4(20,(s0/x(7))/ctx.KMAX,620)*q.bpa.Fmax; cbExact=zeros(N,1);
for k=1:N
 uv=RowVecTrans(tr\eye(4),Fu(k,:)*pre(k));uv=uv/norm(uv);cbExact(k)=uv*diag(1./[Xi(2),Xi(3),Xi(3)])*uv';
end
A=psiSolve(x(7),1,1,s0,zeros(N,1),cbExact,ctx.KMAX,20,620,x(8),ctx.wraps*6,@(r)maxBPAforce(r,'20'));
fprintf(fid,'G1 fully aligned scalar closure max force rel %.16g abs %.16g N effective length %.16g m\n',max(abs(A.Fsingle-Fc)./max(1,Fc)),max(abs(A.Fsingle-Fc)),max(abs((x(7)-A.sCon)-q.activeLength_effective)));
fprintf(fid,'G1 force-effective delta_L differs from initial by %.16g m; local normalize is row-unit normalization, not MATLAB global normalize.\n',max(abs(delta0-x(7)*(q.strain_f-q.strain_p))));
writetable(table(a,Fc,A.Fsingle,cb,cbExact,q.Lmt_p,q.bpa.gama,A.sCon),fullfile(here,'gate1_aligned.csv'));
MA=readmatrix(fullfile(R,'Testing_Data','2022_02_Festo','OpenSim_Vasti_MomentArm.txt'),'FileType','text','NumHeaderLines',7);
oldarm=interp1(MA(:,2),MA(:,5),a,'linear','extrap');newarm=interp1(T.knee_angle_deg,T.arm_m,a,'linear','extrap');
oldForce=hs./oldarm;
writetable(table(a,oldarm,newarm,oldForce,repmat(1294,N,1),hs,oldarm*1294,hnew),fullfile(here,'target_force_arm_audit.csv'));
fprintf(fid,'TARGET legacy inferred force torque/arm %.12g..%.12g N; legacy arm vs new arm max %.12g m; fixed-Fmax using old arm vs new torque max %.12g Nm\n',min(oldForce),max(oldForce),max(abs(oldarm-newarm)),max(abs(oldarm*1294-hnew)));
% G2 fixed-row two-crossing route. Fourth dimension indexes joints AND
% determines N2; use 4x4x2x2 with two poses per joint (no ambiguity).
TT=zeros(4,4,2,2); TT(:,:,1,1)=eye(4);TT(:,:,2,1)=RpToTrans(rotzLocal(.15),[0;0;0]);
TT(:,:,1,2)=RpToTrans(rotzLocal(-.2),[0;-.3;0]); TT(:,:,2,2)=RpToTrans(rotzLocal(.1),[0;-.3;0]);
loc=[.01 .2 .03;.01 .1 .03;.04 -.05 .03;.05 -.2 .03;.06 -.3 .03];
cfg=struct('nPulleyBPA',[2 2],'gain',[1 2],'tackleLineParts',[1 2],'pulleyExitIndex',[2 4]);
try
 b=BiPam_pulley2_20261022('fixed G2',loc,[3 5],20,TT,.7,.7*(1-ctx.KMAX),.1,ctx.fitting,620,Xi(1),Xi(2),Xi(3),ctx.wraps,cfg);
 fprintf(fid,'G2 actual T size=%s Fsingle=%s Torque=%s\n',mat2str(size(b.TransformationMat)),mat2str(size(b.FsingleBPA)),mat2str(size(b.Torque_p)));
 % Independent scalar root using canonical geometric span and cb outputs.
 % This isolates closure algebra, not the independent geometry adapter.
 sb=.7-b.MuscleLength(:)+Xi(1)+.1+2*ctx.fitting; dl=b.deltaL(:,:,b.ActiveCrossing);dl=dl(:);
 cb2=b.DeflClosure_cmB; %#ok<NASGU>
 % cb reconstructed independently with the previous projection on flattened
 % fixed-row geometry, matching its coordinate mapping.
 cbout=zeros(4,1); Fpre=b.Funit; pa=loc(1,:); pb=[9.48 -36.21 30.86]/1000; v=pa-pb;
 rz=rotzLocal(atan2(v(2),v(1))); vv=rz'*v'; yy=atan2(vv(3),vv(1)); ry=[cos(yy) 0 sin(yy);0 1 0;-sin(yy) 0 cos(yy)]; tr=RpToTrans(rz*ry',pb');
 pre=festo4(20,(.7-(b.MuscleLength(:)-.1-2*ctx.fitting-Xi(1)))/.7/ctx.KMAX,620)*b.Fmax;
 for k=1:4
  uv=RowVecTrans(tr\eye(4),Fpre(k,:)*pre(k));uv=uv/norm(uv);cbout(k)=uv*diag(1./[Xi(2),Xi(3),Xi(3)])*uv';
 end
 B=psiSolve(.7,2,2,sb, -dl,cbout,ctx.KMAX,20,620,.1,ctx.wraps,@(r)maxBPAforce(r,'20'));
 fprintf(fid,'G2 max force residual %.16g N; contraction %.16g m; slack matches %d infeasible matches %d\n',max(abs(B.Fsingle-b.FsingleBPA(:))),max(abs(B.sCon-b.sContraction(:))),isequal(B.slack,b.PulleySlack(:)),isequal(B.infeas,b.PulleyInfeasible(:)));
 writetable(table(b.FsingleBPA(:),B.Fsingle,b.sContraction(:),B.sCon,b.Ftendon(:),B.Ftendon,b.PulleySlack(:),B.slack),fullfile(here,'gate2.csv'));
 save(fullfile(here,'gate2_actual_class.mat'),'b','B','TT','loc','cfg');
catch ME
 fprintf(fid,'G2 ERROR %s\n%s\n',ME.identifier,getReport(ME,'extended','hyperlinks','off'));
end
save(fullfile(here,'checkpoint.mat'),'p','q','P','x','Xi','hs','hnew','delta0','ctx');
fprintf(fid,'No optimization: G1 equivalence unestablished. Tolerances: force relative 1e-3, length 1e-6 m, torque 1e-3 Nm.\n');
checkcode(mfilename('fullpath')); disp(fileread(fullfile(here,'checkpoint.txt')));
end
function s=lastpart(p)
[~,s]=fileparts(p);
end
function z=rotzLocal(t)
z=[cos(t) -sin(t) 0;sin(t) cos(t) 0;0 0 1];
end
function P = psiSolve(rest, nB, G, s0x3, span, cb, KMAX, Dia, Pres, ...
    tendon, wraps, FmaxOf)
N = numel(s0x3);
[~, iRef] = max(span);
DL = span(iRef) - span;
Fmax = FmaxOf(rest);
cSpr = 1 / (wraps * 1.51e-6 * 193e9 / tendon);
Ffun = @(s) festo4(Dia, (s / rest) / KMAX, Pres) * Fmax;
sHi = KMAX * rest;
sCon = zeros(N, 1); Fsingle = zeros(N, 1);
slack = false(N, 1); infeas = false(N, 1);
for i = 1:N
    psi = @(s) G * s - G * cb(i) * Ffun(s) - nB * Ffun(s) * cSpr / G ...
        - DL(i) - G * (s0x3(i) - DL(i));
    pl = psi(0); ph = psi(sHi);
    if pl > 0
        sCon(i) = 0; Fsingle(i) = Ffun(0); slack(i) = true;
    elseif ph < 0
        sCon(i) = sHi; Fsingle(i) = Ffun(sHi); infeas(i) = true;
    else
        try
            sCon(i) = max(fzero(psi, [0, sHi]), 0);
        catch
            sCon(i) = 0; slack(i) = true;
        end
        Fsingle(i) = Ffun(sCon(i));
    end
end
Ftendon = nB * Fsingle / G;
Ftendon(slack) = 0;
Ftendon(infeas) = NaN;
P = struct('Ftendon', Ftendon, 'Fsingle', Fsingle, 'sCon', sCon, ...
    'slack', slack, 'infeas', infeas, 'DL', DL);
end

function cb = bracketCb(Location, ctx, Xi, N)
cb = zeros(N, 1);
C = diag(1 ./ [Xi(2), Xi(3), Xi(3)]);
pA = Location(1, :, 1);
Pbr = [9.48, -36.21, 30.86] / 1000;
phbrA = pA - Pbr;
thbrA = atan2(phbrA(2), phbrA(1));
RhbrZ = [cos(thbrA), -sin(thbrA), 0; sin(thbrA), cos(thbrA), 0; 0, 0, 1];
pbrhA = RhbrZ' * phbrA';
thetaY = atan2(pbrhA(3), pbrhA(1));
Ry = [cos(thetaY), 0, sin(thetaY); 0, 1, 0; -sin(thetaY), 0, cos(thetaY)];
Thbr = RpToTrans(RhbrZ * Ry', Pbr');
Tinv = Thbr \ eye(4);
for ii = 1:N
    pt1 = Location(1, :, ii);
    pt2 = pt1;
    for r = 2:size(Location, 1)
        if norm(Location(r, :, ii) - pt1) > 1e-6
            pt2 = Location(r, :, ii);
            if r == ctx.CrossPoint
                pt2 = RowVecTrans(ctx.T_Pam(:, :, ii), pt2);
            end
            break
        end
    end
    u = pt2 - pt1;
    if norm(u) < 1e-9
        cb(ii) = 1 / Xi(2);
        continue
    end
    ub = RowVecTrans(Tinv, u / norm(u));
    cb(ii) = (ub * C) * ub.';
end
end





