function [tau, rotational, translation] = fr_coordinateTorque(cs, kin, mech, passIdx, crossing)
% Generalized force of the class insertion wrench, in the coordinate units.
% The child-origin translation term is necessary for the rolling knee.
pass=cs.passes{passIdx};coord=pass.sweptCoord;N=numel(pass.qStructs);
tau=zeros(N,1);rotational=tau;translation=tau;h=1e-6;
if crossing==1;fn=pass.T1;else;fn=pass.T2;end
for i=1:N
 q=pass.qStructs{i};qp=q;qm=q;qp.(coord)=q.(coord)+h;qm.(coord)=q.(coord)-h;
 T=fn(kin,q);Tp=fn(kin,qp);Tm=fn(kin,qm);
 A=T(1:3,1:3)'*((Tp(1:3,1:3)-Tm(1:3,1:3))/(2*h));
 omega=[A(3,2)-A(2,3);A(1,3)-A(3,1);A(2,1)-A(1,2)]/2;
 velocity=T(1:3,1:3)'*((Tp(1:3,4)-Tm(1:3,4))/(2*h));
 moment=reshape(mech.cell.Torque_ins(i,:,1,crossing),3,1);
 force=reshape(mech.cell.F_ins(i,:,1,crossing),3,1);
 rotational(i)=dot(moment,omega);translation(i)=dot(force,velocity);
 tau(i)=rotational(i)+translation(i);
end
end
