function fr_solveRectusCamPilot(outputDir)
% Reduced-order common-cable closure using the saved robot knee cam contacts.
C = fr_setup();
kin = fr_kinematics(C.osimStock);
file = fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S = load(file,'predBest','ctxCamBest');
source = fullfile(C.camp,'run_20261105_rectus_sizing');
A = jsondecode(fileread(fullfile(source,'rectus_route.json')));
targets = readtable(fullfile(source,'rectus_pose_grid.csv'));
if ~isfolder(outputDir); mkdir(outputDir); end
ctx = S.ctxCamBest;
pred = S.predBest;
hip = unique(targets.hip_deg); hip = hip(1:2:end);
knee = unique(targets.knee_deg); knee = knee(1:2:end);
[H,K] = ndgrid(hip,knee);
n = numel(H);
[~,home] = min(abs(ctx.phiD));
% The exported pose array can collapse unused routing rows. Use the
% physical femoral mounting-base guide in the cam's native seed instead.
entry = pred.routeInfo.seed20(2,:);
fprintf('Femoral entry guide [mm]: %.1f %.1f %.1f\n',1000*entry);
anchor = A.points_m(1,:);          % saved pelvis attachment
q = kin.qDefault;
q.knee_angle_r = 0;
% Provisional home transform for checking the displayed OpenSim tibia.
% The cam solve itself uses the saved robot ICR and native t1 geometry.
homeTibia = fr_chainT(kin,q,'tibia_r','femur_r');
adapter = mldivide(homeTibia,ctx.T_f_t1(:,:,home));
path = cell(n,1);
lengths = zeros(n,1); bay = lengths; hipArm = lengths; kneeArm = lengths;
clipMargin = inf(n,1);
hipCompliance = lengths; kneeCompliance = lengths;
fit = ctx.fitting;
if ~isscalar(fit); fit = fit(1); end
for j = 1:numel(knee)
    [~,ix] = min(abs(ctx.phiD-knee(j)));
    % Follow the saved cam's active contacts and its robot knee transform.
    active = logical(pred.routeInfo.active(:,ix));
    active(1:2) = false; active(end) = true;
    rows = find(active);
    downstream = pred.Location(rows,:,ix);
    t = ctx.T_Pam(:,:,ix);
    distal = rows >= ctx.CrossPoint;
    downstream(distal,:) = apply(t,downstream(distal,:));
    first = find(distal,1);
    cam = [entry;downstream];
    assert(~isempty(first),'Rectus cam has no distal point');
    f = cam(first,:)-cam(first+1,:);
    f = f/norm(f);
    moment = cross(cam(first+1,:)-t(1:3,4)',f);
    % Distal bracket assignment in the robot tibial/ICR axes.
    fd = (t(1:3,1:3)'*f')';
    ck = sum(fd.^2 ./ [C.Xi.Xi1,C.Xi.Xi2,C.Xi.Xi1]);
    ca = clipSamples(cam,.004);
    t1 = ctx.T_f_t1(:,:,ix);
    tib = apply(adapter/t1,ca);
    % Femur contacts upstream and terminal attachment can touch by design.
    % This deliberately checks the free spans, excluding the end 5 mm.
    cm = min([fr_surfaceDistance(C,'femur_r',ca); ...
              fr_surfaceDistance(C,'tibia_r',tib)])-.002;
    for i = 1:numel(hip)
        id = i+(j-1)*numel(hip);
        q = kin.qDefault; q.hip_flexion_r = deg2rad(hip(i));
        pelvisToFemur = fr_chainT(kin,q,'pelvis','femur_r');
        start = apply(pelvisToFemur,anchor);
        points = [start;cam];
        lead = norm(entry-start);
        path{id} = points;
        lengths(id) = sum(vecnorm(diff(points),2,2));
        bay(id) = lead-2*fit;
        dq = 1e-5; qp=q; qm=q;
        qp.hip_flexion_r = q.hip_flexion_r+dq; qm.hip_flexion_r = q.hip_flexion_r-dq;
        plus = apply(fr_chainT(kin,qp,'pelvis','femur_r'),anchor);
        minus = apply(fr_chainT(kin,qm,'pelvis','femur_r'),anchor);
        hipArm(id) = -(norm(entry-plus)-norm(entry-minus))/(2*dq);
        kneeArm(id) = moment(3);
        dir = (entry-start)/lead;
        hipCompliance(id) = sum(dir.^2 ./ [C.Xi.Xi1,C.Xi.Xi2,C.Xi.Xi2]);
        kneeCompliance(id) = ck;
        clipMargin(id) = cm;
    end
end
[ok,idx] = ismember([H(:),K(:)],[targets.hip_deg,targets.knee_deg],'rows');
assert(all(ok),'Missing signed human target poses');
humanH = targets.human_hip_Nm(idx); humanK = targets.human_knee_Nm(idx);
trial = {}; bestByBranch = cell(1,2);
for b = 1:2
    branch = {'rigid','compliant'}; branch = branch{b};
    xi0 = 0; compliance = zeros(n,1);
    if b==2
        xi0 = C.Xi.Xi0_ext;
        compliance = hipCompliance+kneeCompliance;
    end
    best = [];
    for gain = 1:3
        r = fr_travelReeving(gain);
        cap = 1-C.LminFraction.d20;
        minL = range(lengths)/(gain*cap);
        maxL = min(bay+(max(lengths)-lengths)/gain);
        trial{end+1} = struct('branch',branch,'gain',gain,'minRest_mm',minL*1000, ...
            'maxRestProxy_mm',maxL*1000,'axialPossible',maxL>=minL); %#ok<AGROW>
        % Preserve the smallest-stroke result even if no packaged length fits.
        for L0 = unique([minL*1.02, linspace(max(minL*1.02,maxL*.8), ...
                max(minL*1.02,maxL),5)])
            for count = 1:3
                F0 = maxBPAforce(L0,'20',C.pressure);
                for offset = [0 .002 .005]
                    tendon = max(.001,max(lengths)-L0-2*fit-xi0+offset);
                    demand = L0+tendon+2*fit+xi0-lengths;
                    lo = zeros(n,1); hi = ones(n,1)*cap*L0;
                    % Separate tendons and bracket sets on the representative
                    % cam path. Solve each line's compliance before summing.
                    cable = @(s) F0*festo4(20,s/(cap*L0),C.pressure)/gain;
                    slack = demand+compliance.*cable(lo)<0;
                    stroke = gain*hi-demand-compliance.*cable(hi)<0;
                    for it = 1:28
                        mid = (lo+hi)/2;
                        lower = gain*mid-demand-compliance.*cable(mid)<0;
                        lo(lower)=mid(lower); hi(~lower)=mid(~lower);
                    end
                    s = (lo+hi)/2; tension = count*cable(s);
                    tension(slack|stroke) = 0;
                    tauH = tension.*hipArm; tauK = tension.*kneeArm;
                    misses = @(x,y) max(0,abs(y)-x.*sign(y));
                    dH = misses(tauH,humanH); dK = misses(tauK,humanK);
                    axial = bay-(L0-s);
                    % Package gate for the sleeve body and a provisional
                    % open-mesh check for the free cable; tackle geometry unresolved.
                    clearance = min(axial);
                    score = 1e4*(mean(slack|stroke)+max(0,-clearance)) + ...
                        sum(dH)/sum(abs(humanH))+sum(dK)/sum(abs(humanK));
                    if isempty(best) || score<best.score
                        best = struct('branch',branch,'gain',gain,'reeving',r, ...
                            'L0_m',L0,'tendon_m',tendon,'count',count,'score',score, ...
                            'slack',slack,'stroke',stroke,'s_m',s,'tension_N',tension, ...
                            'hip_Nm',tauH,'knee_Nm',tauK,'axialMargin_m',axial, ...
                            'hipMeets',dH<1e-6 & ~slack & ~stroke, ...
                            'kneeMeets',dK<1e-6 & ~slack & ~stroke);
                    end
                end
            end
        end
        if ~isempty(best) && ~any(best.slack|best.stroke) && min(best.axialMargin_m)>=0
            break
        end
    end
    bestByBranch{b}=best;
    sleeveMargin = inf(n,1);
    for id = 1:n
        points = path{id};
        u = points(2,:)-points(1,:); u = u/norm(u);
        first = points(1,:)+fit*u;
        last = first+(best.L0_m-best.s_m(id))*u;
        sample = clipSamples([first;last],.004);
        if ~isempty(sample)
            radius = fr_diameter(C,'20',best.s_m(id)/best.L0_m)/2000;
            sleeveMargin(id) = min(fr_surfaceDistance(C,'femur_r',sample))-radius;
        end
    end
    bestByBranch{b}.centerSleeveFemurMargin_m = sleeveMargin;
    save(fullfile(outputDir,'rectus_cam_pilot.mat'),'bestByBranch','trial','H','K', ...
        'humanH','humanK','path','lengths','bay','clipMargin','-v7');
    try
        emit(outputDir,hip,knee,H,K,humanH,humanK,best,path,clipMargin,sleeveMargin,C,fit,homeTibia);
    catch ME
        warning('fr_solveRectusCamPilot:Export','%s',getReport(ME,'basic'));
    end
end
save(fullfile(outputDir,'rectus_cam_pilot.mat'),'bestByBranch','trial','H','K', ...
    'humanH','humanK','path','lengths','bay','clipMargin','-v7');
end

function p = apply(T,x)
p = (T(1:3,1:3)*x'+T(1:3,4))';
end

function points = clipSamples(path,step)
points = zeros(0,3);
for i = 1:size(path,1)-1
    a=path(i,:); b=path(i+1,:); L=norm(b-a);
    if L<=.010; continue; end
    t=linspace(.005/L,1-.005/L,max(2,ceil(L/step)))';
    points=[points;a+t.*(b-a)]; %#ok<AGROW>
end
end

function emit(dir,hip,knee,H,K,hH,hK,b,path,clip,sleeveMargin,C,fit,homeTibia)
fprintf(['%s: G%d, n%d, L0 %.1f mm, slack %d, stroke %d, minimum axial %.1f mm, ' ...
    'cam cable SDF margin %.1f mm, knee %d/%d, hip %d/%d, both %d/%d\n'], ...
    b.branch,b.gain,b.count,b.L0_m*1000,nnz(b.slack),nnz(b.stroke), ...
    min(b.axialMargin_m)*1000,min(clip)*1000,nnz(b.kneeMeets),numel(H), ...
    nnz(b.hipMeets),numel(H),nnz(b.kneeMeets&b.hipMeets),numel(H));
fprintf('Central sleeve to femur SDF margin %.1f mm. Full hardware clearance unresolved.\n', ...
    min(sleeveMargin)*1000);
T = table(H(:),K(:),hH,b.hip_Nm,hK,b.knee_Nm,b.tension_N,b.slack,b.stroke, ...
    b.axialMargin_m,clip,sleeveMargin,'VariableNames',{'hip_deg','knee_deg','humanHip_Nm', ...
    'robotHip_Nm','humanKnee_Nm','robotKnee_Nm','tension_N','slack','stroke', ...
    'axialMargin_m','camCableProxyMargin_m','centerSleeveFemurMargin_m'});
writetable(T,fullfile(dir,[b.branch '_pose_grid.csv']));
f=figure('Visible','off','Color','w');
tiledlayout(f,1,2);
[~,i]=min(abs(hip)); [~,j]=min(abs(knee));
humanK=reshape(hK,size(H)); robotK=reshape(b.knee_Nm,size(H));
humanH=reshape(hH,size(H)); robotH=reshape(b.hip_Nm,size(H));
nexttile; plot(knee,humanK(i,:),'--',knee,robotK(i,:),'-','LineWidth',1.6);
xlabel('Knee angle (deg)'); ylabel('Signed torque (N m)'); title('Knee at neutral hip'); grid on
nexttile; plot(hip,humanH(:,j),'--',hip,robotH(:,j),'-','LineWidth',1.6);
xlabel('Hip flexion angle (deg)'); title('Hip at extended knee'); grid on
legend({'Stock human','Rectus cam seed'},'Location','best');
sgtitle(sprintf('Rectus femoris %s | common-cable equilibrium',b.branch));
exportgraphics(f,fullfile(dir,[b.branch '_signed_torque.png']),'Resolution',170);
close(f);

f=figure('Visible','off','Color','w');
ax=axes(f); hold(ax,'on');
field=load(fullfile(C.camp,'audit_20261105','surface_fields','femur_r.mat'));
[X,Y,Z]=meshgrid(field.x,field.y,field.z);
bone=isosurface(X,Y,Z,permute(double(field.distance),[2 1 3]),0);
patch(ax,'Faces',bone.faces,'Vertices',1000*bone.vertices, ...
    'FaceColor',[.68 .7 .75],'FaceAlpha',.3,'EdgeColor','none');
field=load(fullfile(C.camp,'audit_20261105','surface_fields','tibia_r.mat'));
[X,Y,Z]=meshgrid(field.x,field.y,field.z);
bone=isosurface(X,Y,Z,permute(double(field.distance),[2 1 3]),0);
patch(ax,'Faces',bone.faces,'Vertices',1000*apply(homeTibia,bone.vertices), ...
    'FaceColor',[.52 .59 .66],'FaceAlpha',.3,'EdgeColor','none');
id=i+(j-1)*numel(hip);
P=path{id}; p=1000*P;
route=plot3(ax,p(:,1),p(:,2),p(:,3),'-o','Color',[0 .35 .8],'LineWidth',1.8);
u=(P(2,:)-P(1,:))/norm(P(2,:)-P(1,:));
sleeve=1000*[P(1,:)+fit*u;P(1,:)+(fit+b.L0_m-b.s_m(id))*u];
axisLine=plot3(ax,sleeve(:,1),sleeve(:,2),sleeve(:,3),'-','Color',[.88 .37 .1],'LineWidth',5);
legend(ax,[route axisLine],{'Cable and cam contacts','BPA sleeve axis'},'Location','best');
xlabel(ax,'Anterior x (mm)'); ylabel(ax,'Superior y (mm)'); zlabel(ax,'Right z (mm)');
axis(ax,'equal'); grid(ax,'on'); view(ax,2);
xlim(ax,[min([p(:,1);sleeve(:,1)])-75,max([p(:,1);sleeve(:,1)])+75]);
title(ax,sprintf('Rectus femoris %s | neutral hip and extended knee',b.branch));
exportgraphics(ax,fullfile(dir,[b.branch '_route.png']),'Resolution',180);
close(f);
end
