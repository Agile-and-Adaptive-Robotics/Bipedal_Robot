function fr_screenVastiCamBranches(out)
% Size three separate 20 mm load paths on the saved robot cam contour.
C=fr_setup();
file=fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S=load(file,'ctxCamBest','predBest');
ctx=S.ctxCamBest;
pred=S.predBest;
H=jsondecode(fileread(fullfile(out,'targets','vasti_individual.json')));
mask=ctx.phiD>=-120-1e-2 & ctx.phiD<=0;
idx=find(mask);
q=ctx.phiD(idx(:)); q=q(:);
parts=cell(size(q));
entry=pred.routeInfo.seed20(2,:);
for k=1:numel(q)
    i=idx(k);
    P=pred.Location(:,:,i);
    distal=ctx.CrossPoint:size(P,1);
    T=ctx.T_Pam(:,:,i);
    P(distal,:)=(T(1:3,1:3)*P(distal,:)'+T(1:3,4))';
    active=logical(pred.routeInfo.active(:,i));
    active(1:2)=false;
    active(end)=true;
    parts{k}=[P(1,:);entry;P(active,:)];
end
fit=ctx.fitting;
if ~isscalar(fit); fit=fit(1); end
names={'vas_int_r','vas_lat_r','vas_med_r'};
rows={};
for j=1:numel(names)
    h=H.(names{j});
    target=interp1(h.knee_deg,h.target_Nm,q,'pchip','extrap');
    patch=fr_patch(C,'Femur','femur_r',h.origin_m,60);
    for br={'rigid','compliant'}
        name=br{1};
        path=fullfile(out,sprintf('%s_cam_femoral_bay_%s.mat',names{j},name));
        if isfile(path)
            R=load(path,'best'); rows{end+1}=R.best; continue %#ok<AGROW>
        end
        best=struct('score',inf);
        window=[inf,NaN,NaN,NaN];
        [U,V]=ndgrid(-.65:.65:.65);
        uv=[U(:),V(:)];
        for placement=[uv,repmat([.055,-.125],size(uv,1),1); ...
                uv,repmat([.075,-.145],size(uv,1),1)]'
            p=patch.project(placement(1),placement(2),C.standoff_mm(2));
            length=zeros(size(q)); bay=length;
            T=zeros(4,4,numel(q),2);
            for k=1:numel(q)
                P=parts{k}; P(1,:)=p; P(2,:)=P(2,:)+[placement(3),placement(4),0];
                step=diff(P); length(k)=sum(vecnorm(step,2,2));
                bay(k)=norm(step(1,:))-2*fit;
                u=step(1,:)/max(norm(step(1,:)),eps);
                T(:,:,k,1)=eye(4);
                T(1:3,4,k,1)=p'+length(k)*u';
            end
            T(:,:,1,2)=eye(4); T(:,:,2,2)=eye(4);
            % Cable torque follows virtual work on the actual cam contact path.
            arm=-gradient(length)./gradient(deg2rad(q));
            minimum=range(length)/(3*.24);
            if minimum-min(bay)<window(1)
                window=[minimum-min(bay),min(bay),range(length),minimum];
            end
            for gain=1:3
                lo=range(length)/(gain*.24);
                hi=min(bay);
                if lo>=hi || hi<=.03; continue; end
                for L=linspace(max(.03,1.03*lo),hi,4)
                    tendon0=max(.025,max(length)-2*fit-L);
                    for extra=[-.003,0,.003,.009]
                        tendon=max(.025,tendon0+extra);
                        cfg=struct('nPulleyBPA',[1 1],'gain',[gain 1], ...
                            'tackleLineParts',[gain 1]);
                        if strcmp(name,'rigid')
                            xi=[0 Inf Inf]; wraps=1e4;
                        else
                            xi=[C.Xi.Xi0_ext,C.Xi.Xi1,C.Xi.Xi2]; wraps=6;
                            cfg.insertionBracketStiffness=[C.Xi.Xi1,C.Xi.Xi2,C.Xi.Xi1];
                        end
                        m=BiPam_pulley2_xi3pose(names{j},[p;0 0 0],[2 2],20, ...
                            T,L,.76*L,tendon,fit,C.pressure,xi(1),xi(2),xi(3),wraps, ...
                            cfg,0,zeros(numel(q),2));
                        force=m.Ftendon(:,1);
                        bad=m.PulleySlack(:,1)|m.PulleyInfeasible(:,1)| ...
                            ~isfinite(force)|force<=1e-6;
                        force(bad)=0;
                        radius=max(fr_diameter(C,'20',m.sContraction(:,1)/L))/2000;
                        guide=parts{1}(2,:)+[placement(3),placement(4),0];
                        A=p+fit*(guide-p)/norm(guide-p);
                        B=A+max(L-m.sContraction(:,1))*(A-p)/norm(A-p);
                        sample=A+linspace(0,1,37)'.*(B-A);
                        margin=min(fr_surfaceDistance(C,'femur_r',sample))-radius-.001;
                        for n=2:5
                            tau=n*force.*arm;
                            deficit=max(0,abs(target)-tau.*sign(target));
                            score=1e3*mean(bad)+1e4*max(0,-margin)+ ...
                                100*sum(deficit)/max(sum(abs(target)),eps)+ ...
                                40*mean(deficit>1e-6)+.05*n;
                            if score<best.score
                                reeving=fr_travelReeving(gain);
                                best=struct('name',names{j},'branch',name,'score',score, ...
                                    'count',n,'gain',gain,'idealSheaves',n*reeving.sheaveCount, ...
                                    'rest_mm',1000*L,'tendon_mm',1000*tendon,'origin_m',p, ...
                                    'femoralGuide_m',guide, ...
                                    'knee_deg',q,'human_Nm',target,'robot_Nm',tau, ...
                                    'fractionMet',mean(deficit<=1e-6 & ~bad), ...
                                    'noForceFraction',mean(bad),'femurSleeveMargin_mm',1000*margin);
                            end
                        end
                    end
                end
                if isfield(best,'noForceFraction') && best.noForceFraction==0 && best.femurSleeveMargin_mm>=0
                    break
                end
            end
        end
        if ~isfield(best,'name')
            best.name=names{j}; best.branch=name;
            best.failure='No nominal length window within the sampled femoral bay';
        end
        best.minimumBayGap_mm=1000*window(1);
        best.femoralBay_mm=1000*window(2);
        best.camExcursion_mm=1000*window(3);
        best.minimumRestAtGain3_mm=1000*window(4);
        save(path,'best','-v7');
        rows{end+1}=best; %#ok<AGROW>
        fprintf('%s %s fixed-cam sizing complete.\n',names{j},name);
    end
end
fid=fopen(fullfile(out,'vasti_cam_femoral_bay_screen.json'),'w');
fprintf(fid,'%s\n',jsonencode(rows)); fclose(fid);
end
