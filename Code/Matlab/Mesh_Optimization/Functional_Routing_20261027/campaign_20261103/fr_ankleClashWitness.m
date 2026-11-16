function fr_ankleClashWitness(old,out)
% Identify the governing tibialis anterior–EHL envelope pair and pose.
C=fr_setup(); kin=fr_kinematics(C.osimStock);
if ~isfolder(out); mkdir(out); end
for branch={'compliant','rigid'}
    S=load(fullfile(old,['assembled_' branch{1} '.mat']),'caseData','keys','T');
    a=S.caseData{3}; b=S.caseData{9};
    best=inf;
    for i=1:numel(a.poses)
        A=primitives(a,i); B=primitives(b,i);
        for j=1:numel(A)
            for k=1:numel(B)
                x=A(j); y=B(k);
                if ~x.rigid && ~y.rigid; continue; end
                g=fr_segmentDistance(x.a,x.b,y.a,y.b)-x.r-y.r-C.clearanceEnvelope_m;
                if g<best
                    best=g;
                    at=i;
                    first=x;
                    second=y;
                end
            end
        end
    end
    match=strcmp(S.T.group_a,'tib_ant') & strcmp(S.T.group_b,'ext_hal');
    assert(abs(1000*best-S.T.rigid_margin_mm(match))<1e-7);
    T=fr_bodyT(kin,a.poses{at},'tibia_r');
    w=struct('margin_mm',1000*best,'angles_deg',S.keys(at,:), ...
        'tib_ant',local(first,T),'ehl',local(second,T));
    w.tib_antParts=arrayfun(@(x)local(x,T),primitives(a,at));
    w.ehlParts=arrayfun(@(x)local(x,T),primitives(b,at));
    fid=fopen(fullfile(out,['tib_ant_ehl_' branch{1} '_clash.json']),'w');
    fprintf(fid,'%s\n',jsonencode(w)); fclose(fid);
    fprintf('%s clash %.2f mm at ankle/subtalar/MTP %.2f / %.2f / %.2f deg: %s against %s\n', ...
        branch{1},1000*best,w.angles_deg,first.name,second.name);
end
end

function Q=primitives(d,i)
Q=struct('a',{},'b',{},'r',{},'name',{},'rigid',{},'channel',{});
for v=d.sleeves{i}
    Q(end+1)=part(v.p1,v.p2,v.r,'BPA vessel',false,v.branch); %#ok<AGROW>
end
[p,r,channel]=fr_fittingSpheres(d.sleeves{i},d.fit);
for j=1:size(p,1)
    name='proximal fitting';
    if mod(j,2)==0; name='distal fitting'; end
    Q(end+1)=part(p(j,:),p(j,:),r(j),name,true,channel(j)); %#ok<AGROW>
end
for j=1:size(d.hardware{i},1)
    name='tackle envelope';
    if j==1; name='shank stop'; end
    h=d.hardware{i}(j,:);
    Q(end+1)=part(h(1:3),h(1:3),h(4),name,true,j); %#ok<AGROW>
end
end

function p=part(a,b,r,name,rigid,channel)
p=struct('a',a,'b',b,'r',r,'name',name,'rigid',rigid,'channel',channel);
end

function p=local(p,T)
p.a=1000*(T(1:3,1:3)'*(p.a'-T(1:3,4)))';
p.b=1000*(T(1:3,1:3)'*(p.b'-T(1:3,4)))';
p.r=1000*p.r;
end
