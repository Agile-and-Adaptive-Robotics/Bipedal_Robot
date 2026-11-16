function G=fr_proximalCoupled(C,kin,cs,targets)
% Same-pose two-coordinate equilibrium for thigh and gastroc routes.
G=[];
if isempty(targets); return; end
if startsWith(targets.coordinate1,'hip_')
    body={'pelvis','femur_r','tibia_r'};
else
    body={'femur_r','tibia_r','calcn_r'};
end
frames=cs.rowOwning;
if ~all(ismember(frames,body)); return; end
first=find(~strcmp(frames,body{1}),1);
second=find(strcmp(frames,body{3}),1);
if isempty(first) || isempty(second) || any(~strcmp(frames(first:second-1),body{2}))
    return
end
a=targets.axis1_deg(:)';
b=targets.axis2_deg(:)';
q1=cell(1,numel(a)); q2=cell(1,numel(b));
for i=1:numel(a); q1{i}=kin.qDefault; q1{i}.(targets.coordinate1)=deg2rad(a(i)); end
for i=1:numel(b); q2{i}=kin.qDefault; q2{i}.(targets.coordinate2)=deg2rad(b(i)); end
[I,J]=ndgrid(1:numel(a),1:numel(b));
poses=q1(I(:));
for i=1:numel(poses); poses{i}.(targets.coordinate2)=q2{J(i)}.(targets.coordinate2); end
pass=cs.passes{1};
pass.sweptCoord=targets.coordinate1;
pass.T1=@(k,q)fr_chainT(k,q,body{2},body{1});
pass.T2=@(k,q)fr_chainT(k,q,body{3},body{2});
pass.qStructs=poses;
cs.passes={pass};
cs.rowUse={true(size(frames))};
cs.rowFrame={frames};
cs.rowExpr={cell(size(frames))};
cs.cross=[first,second];
cs.crossFor={cs.cross};
G=struct('cs',cs,'target',targets,'indices',[I(:),J(:)], ...
    'qg',struct('passIdx',1,'qList',{q1},'qList2',{q2}));
end
