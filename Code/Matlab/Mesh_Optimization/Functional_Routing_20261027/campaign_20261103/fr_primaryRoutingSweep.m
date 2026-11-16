function [tau,active,m,meta]=fr_primaryRoutingSweep(C,kin,cs,d,branch)
% Primary signed torque with force-closure and physical stroke screening.
qg=struct('passIdx',1,'qList',{cs.passes{1}.qStructs});
[m,meta]=fr_mechBuild(C,cs,d,kin,branch,qg);
tau=fr_coordinateTorque(cs,kin,m,1,cs.primaryCrossing);
active=isfinite(m.Ftendon(:)) & m.Ftendon(:)>1e-6 & ...
    ~m.PulleySlack(:) & ~m.PulleyInfeasible(:) & ...
    m.strain_p(:)<=1-m.kmax/m.rest+1e-9;
tau(~active | ~isfinite(tau))=0;
end
