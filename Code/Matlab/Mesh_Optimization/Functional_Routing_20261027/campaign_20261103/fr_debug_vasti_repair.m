function fr_debug_vasti_repair()
C=fr_setup();kin=fr_kinematics(C.osimStock);cases=fr_cases27(C,kin);cs=cases(strcmp({cases.name},'vas_int'));
disp(cs.seedLoc);disp(cs.movable(2).patch);disp(cs.movable(2).patch.project(0,0,5));
disp(cs.xSeed);disp(cs.lb);disp(cs.ub);d=fr_unpackX(cs,cs.xSeed);
[loc,bodies]=fr_locForPass(C,cs,d,kin,1);disp(loc);disp(bodies);
qg=struct('passIdx',1,'qList',{cs.passes{1}.qStructs});[m,meta]=fr_mechBuild(C,cs,d,kin,'rigid',qg);
disp(m.sleeveLength_m);disp(meta.loc);
q=qg.qList{1};for j=1:size(loc,1);T=fr_bodyT(kin,q,bodies{j});disp((T(1:3,1:3)*loc(j,:)'+T(1:3,4))');end
end
