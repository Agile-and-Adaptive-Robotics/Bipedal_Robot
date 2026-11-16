function kin=fr_robotKneeContext(C,kin)
% Use the saved robot femur-to-t1 motion and register the tibia at home.
file=fullfile(C.repo,'Code','Matlab','Mesh_Optimization','Results', ...
    'Vas_Pam_20mm_Result_cam_20261014_0416.mat');
S=load(file,'ctxCamBest');
x=S.ctxCamBest;
[angles,i]=sort(x.phiD(:));
T=x.T_f_t1(:,:,i);
translation=reshape(T(1:3,4,:),3,[])';
theta=unwrap(atan2(reshape(T(2,1,:),[],1),reshape(T(1,1,:),[],1)));
assert(min(angles)<-119.9 && max(angles)>0,'Unexpected robot knee calibration range');
k=struct('translation',pchip(angles,translation'), ...
    'theta',pchip(angles,theta),'adapter',eye(4));
q=kin.qDefault; q.knee_angle_r=0;
home=fr_chainT(kin,q,'tibia_r','femur_r');
k.adapter=fr_robotKneeT(k,0)\home;
assert(norm(fr_robotKneeT(k,0)-home,'fro')<1e-10);
kin.robotKnee=k;
end
