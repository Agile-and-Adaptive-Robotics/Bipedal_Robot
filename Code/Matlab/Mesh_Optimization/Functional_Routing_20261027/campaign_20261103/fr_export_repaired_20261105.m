function fr_export_repaired_20261105()
% Export the exact repaired checkpoints and solved vessels; no guessed CAD map.
C=fr_setup();C.exactSegmentDistance=true;C.fullFittingEnvelopes=true;kin=fr_kinematics(C.osimStock);
fr_verify_repairs_20261105();
fr_finalize_repaired_20261105();
fr_apply_vasti_rule_20261105();
outdir=fullfile(C.camp,'audit_20261105','reoptimization');
S=load(fullfile(outdir,'final_results.mat'),'results');rows={};
qList=cell(1,25);
baseQ=kin.qDefault;coordinateNames=keys(kin.coordDefault);
for j=1:numel(coordinateNames);name=coordinateNames{j};baseQ.(name)=kin.coordDefault(name);end
for k=1:25
 q=baseQ;q.knee_angle_r=deg2rad(-5*(k-1));q.knee_angle_l=q.knee_angle_r;qList{k}=q;
end
kinRows=cell(1,25);
for k=1:25
 W=struct();for b=1:numel(kin.bodies);body=kin.bodies{b};W.(body)=fr_bodyT(kin,qList{k},body);end
 kinRows{k}=struct('knee_deg',-5*(k-1),'W',W);
end
fid=fopen(fullfile(outdir,'matlab_body_frames.json'),'w');assert(fid>0);fwrite(fid,jsonencode(kinRows));fclose(fid);
for i=1:numel(S.results)
 for br={'compliant','rigid'}
  chosen=S.results{i}.(br{1});cs=chosen.caseDef;rep=chosen.globalRep;
  qg=struct('passIdx',1,'qList',{qList});
  [m,meta]=fr_mechBuild(C,cs,chosen.design,kin,br{1},qg);
  caps=fr_sleevePrimitives(C,cs,chosen.design,kin,m,meta);
  env=fr_envelope(C,cs,chosen.design,kin,qList,br{1});
  status='best-infeasible';if rep.feasible;status='constraint-feasible';end
  ids=find(cs.rowUse{1});
  record=struct('group',cs.name,'branch',br{1},'status',status, ...
   'checkpoint',fullfile(outdir,cs.name,['coordinated_' br{1} '.mat']), ...
   'nominal_diameter_mm',str2double(cs.diaStr),'nBPA',cs.seed.nBPA, ...
   'family',cs.txConfig.tackleFamily,'gain',max(cs.txConfig.gain), ...
   'design',chosen.design,'recordedLoc',meta.loc,'assignedBodies',{cs.rowFrame{1}(ids)}, ...
   'fractionMet',rep.fractionPosesMet,'recordedViolation',rep.violation, ...
   'seedSignatureMatches',max(abs(meta.loc-rep.loc),[],'all')<1e-9, ...
   'constraints',rep.cinfo,'terms',rep.terms,'tauRobot',rep.tauPrimary,'tauTarget',rep.tauTarget, ...
   'invalidTorqueCells',rep.invalidTorqueCells, ...
   'primaryDof',cs.passes{1}.sweptCoord,'primaryAngles_rad',cs.passes{1}.th1, ...
   'sleeves',{caps},'fittings',{env.spheres}, ...
   'kneeContraction',m.strain_p,'kneeLength_m',m.sleeveLength_m,'kneeForce_N',m.FsingleBPA, ...
   'maxContraction',1-m.kmax/m.rest,'scope','Repaired bounded candidate, not a qualified hardware design');
  % Independent primary-ROM state makes absence of motion in a knee-only
  % animation distinguishable from zero contraction in the force solver.
  [pm,~]=fr_mechBuild(C,cs,chosen.design,kin,br{1},struct('passIdx',1,'qList',{cs.passes{1}.qStructs}));
  record.primaryContraction=pm.strain_p;record.primaryLength_m=pm.sleeveLength_m;
  record.primaryForce_N=pm.FsingleBPA;
  loadT=max(pm.Ftendon);loadR=max(pm.ReactionFmag);
  for p=2:numel(cs.passes)
   [mp,~]=fr_mechBuild(C,cs,chosen.design,kin,br{1},struct('passIdx',p,'qList',{cs.passes{p}.qStructs}));
   loadT=max(loadT,max(mp.Ftendon));loadR=max(loadR,max(mp.ReactionFmag));
  end
  record.sampledTendonLoadMax_N=loadT;record.sampledBracketReactionMax_N=loadR;
  record.reversingSheaveWorstCase2T_N=2*loadT;
  rows{end+1}=record; %#ok<AGROW>
  fprintf('EXPORT REPAIRED %s %s n%d G%d meet %.3f\n',cs.name,br{1},cs.seed.nBPA,max(cs.txConfig.gain),rep.fractionPosesMet);
 end
end
tempfile=fullfile(outdir,'repaired_routes.json.tmp');
fid=fopen(tempfile,'w');assert(fid>0);fwrite(fid,jsonencode(rows));fclose(fid);
movefile(tempfile,fullfile(outdir,'repaired_routes.json'));
fprintf('Exported %d repaired branches.\n',numel(rows));
end
