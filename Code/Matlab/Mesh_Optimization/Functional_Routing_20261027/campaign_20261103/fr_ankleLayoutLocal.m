function [score,fail]=fr_ankleLayoutLocal(cs,d,e)
names={'tibia','tibial lead','foot vessel','sleeve bay','vessel wrap','shaft lead','tendon length'};
margin=[e.boneMargin_m,e.leadMargin_m,e.footMargin_m,e.distalGuideMargin_m, ...
    e.wrapMargin_m,e.shaftCableMargin_m,d.tendon-e.tendonRequired_m];
fail=names(margin < -1e-4);
if any(e.noForce); fail{end+1}='no-force poses'; end
score=100*mean(e.noForce)+1000*sum(max(0,-margin-1e-4));
if isfield(cs,'assemblySizing') && isfield(e,'stopTurn_deg')
    score=score+.01*(e.stopTurn_deg/90)^2;
end
end
