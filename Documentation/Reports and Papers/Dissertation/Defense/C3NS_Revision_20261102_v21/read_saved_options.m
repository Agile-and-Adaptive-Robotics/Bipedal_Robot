rootDir = fileparts(mfilename('fullpath'));
repoRoot = rootDir;
while ~endsWith(repoRoot, 'Bipedal_Robot')
    repoRoot=fileparts(repoRoot);
end
names={'minimizeFlxPin10_results.mat','minimizeExt10mmX3_results.mat'};
out=struct;
for k=1:numel(names)
    source=fullfile(repoRoot,'Testing_Data','2022_02_Festo',names{k});
    S=load(source);
    info=struct;
    props={'PopulationSize','MaxGenerations','MutationFcn','CrossoverFraction','CrossoverFcn','FunctionTolerance','ParetoFraction','MaxStallGenerations'};
    if isfield(S,'opts')
        for j=1:numel(props)
            val=S.opts.(props{j});
            if isa(val,'function_handle'),val=func2str(val);end
            if iscell(val)
                for q=1:numel(val)
                    if isa(val{q},'function_handle'),val{q}=func2str(val{q});end
                end
            end
            info.(props{j})=val;
        end
    end
    if isfield(S,'output'),info.saved_last_output=S.output;end
    info.source=source;
    out.(sprintf('case%d',k))=info;
end
fid=fopen(fullfile(rootDir,'build','saved_solver_options.json'),'w');
fwrite(fid,jsonencode(out,PrettyPrint=true),'char');fclose(fid);
disp(jsonencode(out));
