function env=fr_retainedEnvelope(C,record,kin,qList)
% Evaluate retained hardware at the CURRENT candidate's actual poses.
% Fixed designs and fixed pose grids make caching exact across iterations.
persistent cache
if isempty(cache);cache=containers.Map('KeyType','char','ValueType','any');end
cs=record.caseDef;
vary='';names=fieldnames(qList{1});
for j=1:numel(names)
 if abs(qList{1}.(names{j})-qList{end}.(names{j}))>1e-10;vary=names{j};break;end
end
for p=1:numel(cs.passes)
 dof=cs.passes{p}.sweptCoord;
 if strcmp(dof,vary)||(startsWith(vary,'hip_')&&startsWith(dof,'hip_'))
  cs.envelopePassIdx=p;break
 end
end
key=jsonencode(struct('name',cs.name,'configuration',cs.txConfig,'diameter',cs.diaStr,'count',cs.seed.nBPA,'design',record.design,'branch',record.xiBranch,'pass',cs.envelopePassIdx,'poses',{qList}));
if isKey(cache,key);env=cache(key);return;end
env=fr_envelope(C,cs,record.design,kin,qList,record.xiBranch);cache(key)=env;
end
