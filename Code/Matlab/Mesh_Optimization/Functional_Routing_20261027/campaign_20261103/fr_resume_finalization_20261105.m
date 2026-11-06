function fr_resume_finalization_20261105()
% Recover the completed 27-group sweep without repeating solved work.
C=fr_setup();root=fullfile(C.camp,'audit_20261105','reoptimization');
if ~isfile(fullfile(root,'finalization_complete.json'))&&isfile(fullfile(root,'coordinated_progress.mat'))
 P=load(fullfile(root,'coordinated_progress.mat'),'branch','start');
 if strcmp(P.branch,'compliant')&&P.start==25
  fr_apply_vasti_rule_20261105();fr_export_repaired_20261105();return
 end
end
if ~isfile(fullfile(root,'finalization_complete.json'))
 S=load(fullfile(root,'progress.mat'),'results','wave');results=S.results;
 assert(numel(results)==27&&all(~cellfun(@isempty,results)),'Incomplete initial sweep');
 for i=1:27;assert(isfield(results{i},'rigid')&&isfield(results{i},'compliant'));end
 save(fullfile(root,'final_results.mat'),'results','-v7');
 fprintf('Recovered %d groups at wave %d; proceeding to final-model refinement.\n',numel(results),S.wave);
end
fr_export_repaired_20261105();
end
