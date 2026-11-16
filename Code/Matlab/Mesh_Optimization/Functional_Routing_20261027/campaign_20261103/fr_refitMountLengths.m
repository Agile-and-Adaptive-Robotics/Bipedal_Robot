function d=fr_refitMountLengths(C,kin,cs,d,branch,G,restTarget)
% Use the BiPulley endpoint residuals to bound the tendon length for a
% handful of BPA rest lengths. The caller re-solves the chosen pair.
if isempty(G)
    route=cs;
    qg=struct('passIdx',1,'qList',{cs.passes{1}.qStructs});
else
    route=G.cs;
    qg=G.qg;
end
if isfield(qg,'resolvedRoute')
    P=qg.resolvedRoute.points; B=qg.resolvedRoute.frames;
else
    [P,B]=fr_locForPass(C,route,d,kin,1);
    qg.resolvedRoute=struct('points',P,'frames',{B});
end
m=fr_mechBuild(C,route,d,kin,branch,qg);
n=numel(cs.movable);
bay=sum(vecnorm(diff(P(1:cs.distalSleeveGuideRow,:)),2,2))-2*cs.fittingLength;
loL=cs.lb(3*n+1)/10; hiL=cs.ub(3*n+1)/10;
rest=unique([linspace(loL,min(hiL,max(loL,bay)),9),min(hiL,max(loL,d.L0))]);
tendonLimit=cs.ub(3*n+2)/10;
if isfield(cs,'assemblySizing'); tendonLimit=max(tendonLimit,d.tendon+.050); end
poses=route.passes{1}.qStructs;
span=fr_spans(C,kin,P,B,poses,'foot');
tendonFloor=max(cs.lb(3*n+2)/10,span.requiredSpan_m);
k=m.ActiveCrossing;
g=m.cell.PulleyGain(k);
count=m.cell.NPulleyBPA(k);
cap=1-C.LminFraction.(['d' cs.diaStr]);
lengths=m.cell.MuscleLength;
delta=m.cell.deltaL(:,:,k);
if nargin>=7
    effectivePath=lengths(:)+(1-1/g)*delta(:);
    minRest=max(loL,range(effectivePath)/cap);
    if isfield(cs,'bpaBendEndRow')
        bendEnd=sum(vecnorm(diff(P(1:cs.bpaBendEndRow,:)),2,2));
        minRest=max(minRest,(bendEnd-cs.fittingLength)/(1-cap));
    end
    maxRest=min(hiL,bay);
    if minRest<=maxRest
        rest=min(maxRest,max(minRest,restTarget));
    end
end
cb=m.cell.cbBracket(:);
cAdd=m.cell.DeflClosure_cAdd(:);
if isempty(cb); cb=zeros(numel(lengths),1); end
if isempty(cAdd); cAdd=zeros(numel(lengths),1); end
base=g*(lengths(:)-m.cell.Xi0-2*cs.fittingLength-m.cell.deltaRoute(:))+(g-1)*delta(:);
stiffnessPerLength=m.kSpr*d.tendon;
best=inf;
for L=rest
    f0=maxBPAforce(L,cs.diaStr,C.pressure);
    a=festo4(cs.diaNum,0,C.pressure)*f0;
    b=festo4(cs.diaNum,1,C.pressure)*f0;
    fcoeff=cAdd-g*cb;
    lower=max((base-g*L+fcoeff*a)/(g+count*a/(g*stiffnessPerLength)));
    upper=min((base-g*L+g*cap*L+fcoeff*b)/(g+count*b/(g*stiffnessPerLength)));
    t=min(tendonLimit,max(tendonFloor,lower+.00025));
    mismatch=max([0,lower-t,t-upper,tendonFloor-t,L-bay]);
    demand=g*t-(base-g*L);
    % Prefer modest peak contraction where the lengths fit.
    proxy=1e4*mismatch+max(0,max(demand))/(g*cap*L);
    if proxy<best; best=proxy; bestL=L; bestT=t; end
end
d.L0=bestL; d.tendon=bestT;
end
