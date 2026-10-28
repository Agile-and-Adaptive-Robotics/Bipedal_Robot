function audit_builder
here=fileparts(mfilename('fullpath')); root=fileparts(fileparts(here)); cd(here);
diary('builder_comparison.txt'); addpath(root);
load_system(fullfile(root,'demos','KneeReflexDemo.slx')); original=signature('KneeReflexDemo'); close_system('KneeReflexDemo',0);
mkdir('builder_snapshot'); copyfile(fullfile(root,'demos','sns_build_demo.m'),fullfile(here,'builder_snapshot','sns_build_demo.m'));
run(fullfile(here,'builder_snapshot','sns_build_demo.m'));
built=signature('KneeReflexDemo'); cd(here);
same=isequaln(original,built); fprintf('Exact functional block-parameter and input-edge comparison: %d\n',same);
save('builder_comparison.mat','original','built','same');
for k=1:numel(original)
 if k>numel(built)||~isequaln(original(k),built(k)), disp(original(k)); if k<=numel(built),disp(built(k));end,end
end
for nm={'sns_demo_results.mat','knee_twin_ref.mat'}
 fprintf('Existing result %s\n',nm{1}); disp(whos('-file',fullfile(root,'results',nm{1})));
end
close_system('KneeReflexDemo',0); diary off;
end
function s=signature(m)
b=sort(find_system(m,'LookUnderMasks','all','FollowLinks','on','Type','Block')); s=struct([]);
for k=1:numel(b)
 s(k).path=erase(b{k},[m '/']); p=get_param(b{k},'DialogParameters'); q=struct;
 if ~isempty(p), n=sort(fieldnames(p));for j=1:numel(n),try,q.(n{j})=get_param(b{k},n{j});catch,end,end,end
 s(k).parameters=q; ph=get_param(b{k},'PortHandles'); edges={};
 for j=1:numel(ph.Inport)
 l=get_param(ph.Inport(j),'Line'); if l<0,edges{j}='unconnected';else
 src=get_param(l,'SrcPortHandle'); edges{j}=sprintf('%s:%d',erase(getfullname(get_param(src,'Parent')),[m '/']),get_param(src,'PortNumber'));end
 end
 s(k).inputs=edges;
end
end
