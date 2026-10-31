function regrade_capability_active_20261028
R=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(genpath(fullfile(R,'Code','Matlab')));
cap=fullfile(R,'tmp','paper_bridge','capability_table_20261022');
res=fullfile(R,'Code','Matlab','Mesh_Optimization','Results','BiPulley_surveys');
T=readtable(fullfile(cap,'BiPulley_survey2_perMuscle_20261023.csv'),'TextType','string');
forty=["ercspn_r","semimem_r","vas_int_r","soleus_r","med_gas_r","tib_post_r","glut_max_r","add_mag_r"];
rows={};
for i=1:height(T)
 site=char(T.site(i));
 [~,bn,bx]=fileparts(T.matFile(i)); up=fullfile(res,bn+bx);
 if any(T.site(i)==forty), pat=['BiPulley_survey3_orig_' site '_nB4_G1_*.mat']; else, pat=['BiPulley_survey2_' site '_nB1_G1_*.mat'];end
 ds=dir(fullfile(res,pat));[~,order]=sort({ds.name}); ds=ds(order); orig=fullfile(ds(end).folder,ds(end).name);
 for side=1:2
  if side==1,p=orig;tag='original';else,p=up;tag='updated';end
  S=load(p); o=S.g.obj; dom=logical(S.dom); U=load(up,'surveyMeta'); slot=U.surveyMeta.slot;
  axisU=S.spec.kin.axes(slot,:);axisU=axisU/norm(axisU);
  % Replay the frozen design through the class of record. No optimizer or replacement mechanics.
  c=BiPam_pulley2_20261022(o.Name,o.Location,o.Cross,o.Diameter,o.TransformationMat,o.RestingL,o.Kmax,o.TendonL,o.FittingLength,o.Pressure,o.Xi0,o.Xi1,o.Xi2,o.Wraps,S.spec.pulley);
  tau=reshape(sum(c.Torque_ins(:,:,:,slot).*reshape(axisU,1,3,1),2),size(dom));
  assert(isequal(isfinite(tau),isfinite(S.g.tauB)),'finite mask changed');
  delta=max(abs(tau(isfinite(tau))-S.g.tauB(isfinite(tau))));
  assert(delta<1e-9,'saved torque replay changed: %s %g',site,delta);
  assert(isequal(c.PulleySlack,S.g.slack)&isequal(c.PulleyInfeasible,S.g.infeas),'flags changed');
  [active,fs]=biPulleyActiveForceMask_20261028(c,slot);

  meet=logical(S.g.meet);
  assert(~any(meet(dom)&~active(dom)),'meet without active');
  row=struct('site',site,'side',tag,'mat',char(p),'old_active',100*sum(S.g.active(dom))/sum(dom(:)),'active',100*sum(active(dom))/sum(dom(:)),'meet',100*sum(meet(dom))/sum(dom(:)),'domain',sum(dom(:)),'active_cells',sum(active(dom)),'meet_cells',sum(meet(dom)),'max_torque_delta',delta,'min_force',min(fs(dom)),'max_force',max(fs(dom)),'min_legacy_strain',min(c.strain_p),'max_legacy_strain',max(c.strain_p),'min_solved_contraction',min(c.sContraction(:)),'max_solved_contraction',max(c.sContraction(:)),'rest',c.RestingL,'nBPA',S.surveyMeta.nBPA,'G',S.surveyMeta.G,'diameter',c.Diameter);
  rows{end+1}=row;
  fprintf('%s %s active %.4f -> %.4f meet %.4f delta %.3g\n',site,tag,row.old_active,row.active,row.meet,delta);
 end
end
out=fullfile(R,'tmp','torque_comparison_20261028');if ~isfolder(out),mkdir(out);end
rows=[rows{:}]; writetable(struct2table(rows),fullfile(out,'canonical_replay.csv'));
fid=fopen(fullfile(out,'canonical_replay.json'),'w');fprintf(fid,'%s',jsonencode(rows));fclose(fid);
fprintf('ALL 54 frozen designs replayed.\n');
end




