function fr_checkEDLMatchingPose(outputDir)
% Reproduce the historical clearance result and replay the same full pose.
C=fr_setup();
C.exactSegmentDistance=true;
C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);
stage=fullfile(C.camp,'run_20261105_distal_stages');
S=load(fullfile(stage,'ankle_guide_arrow','ext_dig_compliant.mat'),'result');
c=S.result.chosen;
branch='compliant';
cs=fr_densePilotCase(C,kin,c.caseDef);
legacy=cell(1,numel(cs.passes));
for p=1:numel(cs.passes)
    pass=cs.passes{p};
    poses=pass.qStructs;
    [m,d]=fr_mechBuild(C,cs,c.design,kin,branch,struct('passIdx',p,'qList',{poses}));
    sleeves=fr_sleevePrimitives(C,cs,c.design,kin,m,d);
    cc=cs; cc.cross=cs.crossFor{p};
    ob=fr_obstacles(C,cc,d.loc,d.rowBodies,kin,poses,{},sleeves);
    b=ob.boneWitness;
    i=b.poseIndex;
    q=poses{i};
    legacy{p}=struct('name',pass.name,'pose',q,'bone',witness(kin,b,q,sleeves{i}), ...
        'force_N',m.Ftendon(i),'strain',m.strain_p(i),'length_m',m.sleeveLength_m(i), ...
        'keptBodies',{cs.rowOwning(cs.rowUse{p})},'omittedBodies',{cs.rowOwning(~cs.rowUse{p})});
end
[~,p]=min(cellfun(@(r) r.bone.margin_m,legacy));
old=legacy{p}; q=old.pose;
ankle=rad2deg(q.ankle_angle_r); mtp=rad2deg(q.mtp_angle_r);
G=fr_ankleCoupledSetup(C,kin,c.caseDef,2,{[ankle,ankle-1e-3],[mtp,mtp-1e-3]});
[e,d]=fr_evalCoupledAnkle(C,kin,c.caseDef,c.design,branch,2,G);
assert(abs(e.ankle_deg(1)-ankle)<1e-9 && abs(e.second_deg(1)-mtp)<1e-9);
assert(isequaln(q,e.poses{1}),'Full poses differ');
ob=fr_obstacles(C,G.cs,d.loc,d.frames,kin,{q},{},{d.sleeves{1}});
b=witness(kin,ob.boneWitness,q,d.sleeves{1});
full=struct('pose',q,'bone',b,'force_N',e.tendonForce_N(1), ...
    'strain',d.strain(1),'length_m',d.length_m(1),'keptBodies',{d.owning});
% Restore the missing foot guides in the old single-crossing representation
% at this exact pose. Express each calcaneus guide in the toes frame.
restored=cs;
for row=find(~cs.rowUse{p})
    restored.rowUse{p}(row)=true;
    restored.rowFrame{p}{row}='toes_r';
    restored.rowExpr{p}{row}=fr_chainT(kin,q,cs.rowOwning{row},'toes_r');
end
[m,a]=fr_mechBuild(C,restored,c.design,kin,branch,struct('passIdx',p,'qList',{{q,q}}));
caps=fr_sleevePrimitives(C,restored,c.design,kin,m,a);
restored.cross=restored.crossFor{p};
ob=fr_obstacles(C,restored,a.loc,a.rowBodies,kin,{q},{},caps(1));
guideProbe=struct('boneMargin_m',ob.boneMargin,'strain',m.strain_p(1), ...
    'force_N',m.Ftendon(1),'length_m',m.sleeveLength_m(1));
rows={struct('name','Old one-coordinate route','second',e.second,'bone',old.bone), ...
    struct('name','Full route at the same pose','second',e.second,'bone',full.bone)};
if ~isfolder(outputDir); mkdir(outputDir); end
writeJson(fullfile(outputDir,'CLEARANCE_WITNESSES.json'),rows);
writeJson(fullfile(outputDir,'MATCHED_POSE.json'), ...
    struct('old',old,'full',full,'legacy',{legacy},'guideProbe',guideProbe));
save(fullfile(outputDir,'matched_pose.mat'),'legacy','full','old','e','guideProbe','-v7');
fprintf('Old %s: ankle %+.3f, MTP %+.3f, bone %+.3f mm, %.3f mm radius, strain %.4f, %.2f N\n', ...
    old.name,ankle,mtp,1000*old.bone.margin_m,1000*old.bone.radius_m,old.strain,old.force_N);
fprintf('Full same pose: bone %+.3f mm, %.3f mm radius, strain %.4f, %.2f N\n', ...
    1000*b.margin_m,1000*b.radius_m,full.strain,full.force_N);
fprintf('Old model with calcaneus guides restored: bone %+.3f mm, strain %.4f, %.2f N\n', ...
    1000*guideProbe.boneMargin_m,guideProbe.strain,guideProbe.force_N);
end

function b=witness(kin,b,q,capsules)
cap=capsules(b.capsuleIndex);
T=fr_bodyT(kin,q,b.body);
b.sleeve_m=(T(1:3,1:3)'*([cap.p1;cap.p2]'-T(1:3,4)))';
direction=diff(b.sleeve_m);
along=dot(b.closest_m-b.sleeve_m(1,:),direction)/norm(direction);
b.part='inflated vessel';
if along < -1e-8; b.part='proximal fitting envelope'; end
if along > norm(direction)+1e-8; b.part='distal fitting envelope'; end
b.ankle_deg=rad2deg(q.ankle_angle_r);
b.second_deg=rad2deg(q.mtp_angle_r);
end

function writeJson(file,record)
fid=fopen(file,'w'); assert(fid>=0);
close=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(record,'PrettyPrint',true));
end
