function G = fr_ankleCoupledSetup(C,kin,cs,n,angles)
% Ankle is crossing 1. Subtalar or MTP is crossing 2.
map=fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Routing27_v2','actuator_map_27.json');
rom=jsondecode(fileread(map)).rom_table_deg;
q0=kin.qDefault;
if any(ismember(cs.rowOwning,{'toes_r'}))
    distal='toes_r'; second='mtp_angle_r';
else
    distal='calcn_r'; second='subtalar_angle_r';
end
cc=cs;
frames=cs.rowOwning;
expr=cell(size(frames));
for r=1:numel(frames)
    if ~ismember(frames{r},{'tibia_r','talus_r',distal})
        expr{r}=fr_chainT(kin,q0,frames{r},'talus_r');
        frames{r}='talus_r';
    end
end
crossing=[find(~strcmp(frames,'tibia_r'),1),find(strcmp(frames,distal),1)];
assert(numel(crossing)==2 && all(diff(crossing)>=0),'Missing distal crossing');
assert(all(strcmp(frames(1:crossing(1)-1),'tibia_r')) && ...
    all(strcmp(frames(crossing(1):crossing(2)-1),'talus_r')) && ...
    all(strcmp(frames(crossing(2):end),distal)),'Noncontiguous ankle route frames');
ax=linspace(rom.ankle_angle_r(1),rom.ankle_angle_r(2),n);
bx=linspace(rom.(second)(1),rom.(second)(2),n);
if nargin>=5
    ax=angles{1}(:)'; bx=angles{2}(:)';
end
q1=cell(1,numel(ax)); q2=cell(1,numel(bx));
for i=1:numel(ax)
    q1{i}=q0; q1{i}.ankle_angle_r=deg2rad(ax(i));
end
for i=1:numel(bx)
    q2{i}=q0; q2{i}.(second)=deg2rad(bx(i));
end
[I,J]=ndgrid(1:numel(ax),1:numel(bx));
poses=q1(I(:));
for i=1:numel(poses); poses{i}.(second)=q2{J(i)}.(second); end
pass=cs.passes{1};
pass.sweptCoord='ankle_angle_r';
pass.T1=@(k,q) fr_chainT(k,q,'talus_r','tibia_r');
pass.T2=@(k,q) fr_chainT(k,q,distal,'talus_r');
pass.qStructs=poses;
cc.passes={pass};
cc.rowUse={true(size(frames))};
cc.rowFrame={frames}; cc.rowExpr={expr};
cc.cross=crossing; cc.crossFor={crossing}; cc.envelopePassIdx=1;
G=struct('cs',cc,'indices',[I(:),J(:)],'second',second, ...
    'ankle_deg',ax(I(:))','second_deg',bx(J(:))', ...
    'qg',struct('passIdx',1,'qList',{q1},'qList2',{q2}));
end
