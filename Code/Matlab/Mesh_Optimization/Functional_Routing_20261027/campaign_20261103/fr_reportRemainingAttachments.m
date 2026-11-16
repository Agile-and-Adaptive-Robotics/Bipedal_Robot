function fr_reportRemainingAttachments(out)
% Export winning nominal attachment displacements in their owning frames.
files=[dir(fullfile(out,'*_rigid.mat'));dir(fullfile(out,'*_compliant.mat'))];
rows=cell(0,11);
for k=1:numel(files)
    [~,stem]=fileparts(files(k).name);
    if ~isfile(fullfile(out,[stem '.json'])); continue; end
    file=fullfile(out,'paired_refit',files(k).name);
    if ~isfile(file); file=fullfile(out,files(k).name); end
    S=load(file,'chosen','record');
    cs=S.chosen.caseDef;
    r=S.record;
    ids=find(cs.rowUse{1});
    assert(isequal(r.bodies,cs.rowOwning(ids)));
    for j=1:numel(ids)
        if ~any([cs.movable.row]==ids(j)); continue; end
        start=1000*cs.seedLoc(ids(j),:);
        finish=1000*r.points_m(j,:);
        rows(end+1,:)={r.group,r.branch,ids(j),r.bodies{j}, ...
            start(1),start(2),start(3),finish(1),finish(2),finish(3),norm(finish-start)}; %#ok<AGROW>
    end
end
assert(size(unique(string(rows(:,1:2)),'rows'),1)==34);
T=cell2table(rows,'VariableNames',{'group','branch','route_row','body', ...
    'seed_x_mm','seed_y_mm','seed_z_mm','selected_x_mm','selected_y_mm','selected_z_mm','movement_mm'});
writetable(T,fullfile(out,'remaining20_attachment_movements.csv'));
end
