function fr_validate_repaired_layouts_20261105(results)
% Branch-consistent all-group validation; cells avoid schema-dependent appends.
C=fr_setup();C.exactSegmentDistance=true;C.fullFittingEnvelopes=true;
kin=fr_kinematics(C.osimStock);root=fullfile(C.camp,'audit_20261105','reoptimization');
summary={};dense={};
for br={'rigid','compliant'}
 branch=br{1};envs=cellfun(@(r)r.(branch).envelope,results,'UniformOutput',false);
 for i=1:numel(results)
  chosen=results{i}.(branch);cs=chosen.caseDef;other=envs;other(i)=[];
  opts=struct('quiet',true,'envelopes',{other});
  [~,rep]=fr_objective(C,cs,kin,fr_packX(cs,chosen.design),branch,opts);
  chosen.globalRep=rep;results{i}.(branch)=chosen;
  save(fullfile(root,cs.name,['coordinated_' branch '.mat']),'chosen','-v7');
  summary{end+1}=struct('group',cs.name,'branch',branch,'nBPA',cs.seed.nBPA,'diameter',cs.diaStr, ...
   'gain',max(cs.txConfig.gain),'meetFraction',rep.fractionPosesMet,'violation',rep.violation, ...
   'feasible',rep.feasible,'topViolator',rep.topViolator); %#ok<AGROW>
  csD=cs;
  for p=1:numel(csD.passes)
   pass=csD.passes{p};old=pass.th1(:);v=linspace(min(old),max(old),25)';pass.th1=v';
   if numel(pass.th2)==numel(old);pass.th2=interp1(old,pass.th2(:),v,'linear','extrap')';
   else;pass.th2=repmat(pass.th2(1),1,25);end
   pass.qStructs=cell(1,25);
   for k=1:25
    q=kin.qDefault;q.(pass.sweptCoord)=v(k);
    if ~isempty(pass.sweptCoord2);q.(pass.sweptCoord2)=pass.th2(k);end
    pass.qStructs{k}=q;
   end
   csD.passes{p}=pass;
  end
  csD.target=fr_targets(C,kin,csD);
  [~,rd]=fr_objective(C,csD,kin,fr_packX(csD,chosen.design),branch,opts);
  [m,~]=fr_mechBuild(C,csD,chosen.design,kin,branch,struct('passIdx',1,'qList',{csD.passes{1}.qStructs}));
  dense{end+1}=struct('group',cs.name,'branch',branch,'gridN',25,'targetIsResampled',true, ...
   'meetFraction',rd.fractionPosesMet,'feasible',rd.feasible,'violation',rd.violation,'topViolator',rd.topViolator, ...
   'tauJumpMax_Nm',max(abs(diff(rd.tauPrimary))),'sleeveLengthJumpMax_mm',1000*max(abs(diff(m.sleeveLength_m))), ...
   'contractionMin',min(m.strain_p),'contractionMax',max(m.strain_p), ...
   'invalidTorqueFraction',mean(rd.invalidTorqueCells)); %#ok<AGROW>
  fprintf('FINAL VALIDATE %s %s coarse %.3f dense %.3f violation %.4g\n',cs.name,branch,rep.fractionPosesMet,rd.fractionPosesMet,rd.violation);
 end
end
save(fullfile(root,'final_results.mat'),'results','summary','dense','-v7');
writejson(fullfile(root,'summary.json'),summary);writejson(fullfile(root,'dense_verification.json'),dense);
writejson(fullfile(root,'finalization_complete.json'),struct('complete',true,'coordinatedRounds',1,'denseGrid',25, ...
 'branchesUseMatchingOtherGroups',true,'vas_int_20mm_only',true));
end
function writejson(file,data)
fid=fopen(file,'w');assert(fid>0);fwrite(fid,jsonencode(data));fclose(fid);
end
