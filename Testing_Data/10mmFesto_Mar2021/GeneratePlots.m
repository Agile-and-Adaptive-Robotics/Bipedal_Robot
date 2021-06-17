%Read PAM Flexion Data into Matlab from Data Set provided by Connor
Flexor_Pam = readmatrix('Lindie Muscle Force and Torque Data 10 mm.xlsx','Sheet', 'Flexor Pam', 'Range','A4:M103');
FAngle_Expected = Flexor_Pam(:,13);
FTorque_Expected =Flexor_Pam(:,7);

%Read Flexion Data into Matlab from Data Set generated on Test Jig
Flexion = readmatrix('TorqueByAngle_Mar21.xlsx','Sheet', 'Flexion', 'Range','A2:C16');
FAngle_Actual= Flexion(:,2);
FTorque_Actual =Flexion(:,3);

%Plote two flexion data sets on one new figure
figure
plot(FAngle_Expected,FTorque_Expected, FAngle_Actual, FTorque_Actual, '*')
xlabel('Knee Flexion/Extension, Degree')
ylabel('Torque, Nm')
xlim([-120 20])
ylim([-20 0])
legend('Model','Experimental');
[t,s] = title ('10mm PAM - Bicep Femoris Short Head', 'Z-Torque');
%t.Fontsize = 16;
%s.FontAngle = 'italic';


%Read PAM Extension Data into Matlab from Data Set provided by Connor
Extensor_Pam = readmatrix('Lindie Muscle Force and Torque Data 10 mm.xlsx','Sheet', 'Extensor Pam', 'Range','A4:M103');
EAngle_Expected = Extensor_Pam(:,13);
ETorque_Expected = Extensor_Pam(:,7);

%Read Extension Data into Matlab from Data Set generated on Test Jig
Extension = readmatrix('TorqueByAngle_Mar21.xlsx','Sheet', 'Extension', 'Range','A2:C10');
EAngle_Actual = Extension(:,2);
ETorque_Actual = Extension(:,3);

%Plote two Extension data sets on one new figure 
figure
plot(EAngle_Expected, ETorque_Expected, EAngle_Actual, ETorque_Actual,'*')
xlabel('Knee Flexion/Extension, Degree')
ylabel('Torque, Nm')
xlim([-150 20])
ylim([-20 20])
legend('Model','Experimental');
[t,s] = title ('10 mm PAM - Vastus Group', 'Z-Torque');
%t.Fontsize = 16;
%s.FontAngle = 'italic';