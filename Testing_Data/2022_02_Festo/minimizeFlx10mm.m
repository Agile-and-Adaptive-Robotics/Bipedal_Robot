
%Minimization scheme

clear; clc; close all

[a, b, ~] = minimizeFlx(0,Inf,Inf);         %Get current goodness of fit measures with no extra length and infinite bracket stiffness

%% Use solution from optimizer and check validity on biomimetic knee

load minimizeFlxPin10_results.mat sol_actual

g = sol_actual(1,1:3);
[u,v,bpa] = minimizeFlx(g(1),g(2),g(3));           % Now pull bpa structures out
%% Plot torque curves, Optimized and validation 
load ForceStrainForFit.mat z
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation

str = ["Optimization"; "Validation"];
for i = 1:2
    %Find old old torque, then plot everything
    clear Yq Xq Vq Fold Fq Mold
    Yq = bpa(i).strain./((bpa(i).rest-bpa(i).Kmax)/bpa(i).rest);
    Xq = bpa(i).P;
    Vq = zeros(size(bpa(i).unitD,1),1);
    for j = 1:size(bpa(i).unitD, 1)
        if bpa(i).strain(j) >= 0 && Yq(j) <=1
            Vq(j) = interp2(X, Y, z, Xq, Yq(j));
%             Vq(j) = f10(Yq(j),Xq);
        elseif Yq(j)>1
            Vq(j) = 0;
        elseif bpa(i).strain(j) < 0
            Vq(j) = NaN;
        end
    end
    if bpa(i).dBPA == 20
        Vq = 1500/630*Vq;
    end
    Fold = Vq.*bpa(i).unitD;    %Force vector
    Fq = (Fold(:,1).^2+Fold(:,2).^2).^(1/2);
    Mold = -bpa(i).mA.*Fq;
    
    
    figure
    ax = gca;
    hold on
    scatter(bpa(i).Aexp,bpa(i).Mexp,[],'filled','DisplayName','Experiment')
    scatter(bpa(i).A_h,bpa(i).M_h,[],'filled','DisplayName','Hybrid')
    plot(bpa(i).Ak,bpa(i).M_p(:,3),'DisplayName','New predict')
    plot(bpa(i).Ak,bpa(i).M,'DisplayName','Predict original') %"original" is with updated BPA characterization
    plot(bpa(i).Ak,Mold,'DisplayName','Predict old') %"old" is before updated BPA characterization
    hold off
    title(str(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Torque, N\cdotm')
    legend
end

%% Plot muscle length, optimization and validation
for i = 1:2
    figure
    ax = gca;
    Lm = bpa(1).Lmt-2*bpa(i).fitn-bpa(i).ten;      %Original predicted muscle length
    Lm_p = bpa(1).Lmt_p-2*bpa(i).fitn-bpa(i).ten;    %Optimized muscle length (Lmt_p uses sol_actual(1)
    hold on
    scatter(bpa(i).A_h,bpa(i).Lm_h,[],'filled','DisplayName','Measured')
    plot(bpa(i).Ak,Lm_p,'DisplayName','New predict')
    plot(bpa(i).Ak,Lm,'DisplayName','Predict original')
    hold off
    title(str(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end

%% Plot moment arm, optimization and validation
for i = 1:2
    figure
    ax = gca;
    G_p = (bpa(i).mA_p(:,1).^2+bpa(i).mA_p(:,2).^2).^(1/2);      %z-axis moment arm for optimized
    hold on
    scatter(bpa(i).A_h,bpa(i).mA_h,[],'filled','DisplayName','Measured')
    plot(bpa(i).Ak,G_p,'DisplayName','New predict')
    plot(bpa(i).Ak,bpa(i).mA,'DisplayName','Predict original')
    hold off
    title(str(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end
