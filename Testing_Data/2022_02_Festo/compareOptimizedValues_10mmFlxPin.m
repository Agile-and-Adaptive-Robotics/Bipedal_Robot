%compare values of optimized and validation from 10mm Flx_pinned pareto
%front

PAR = [Pareto_front, Pareto_Fvals];                     %Combine Pareto front and fvals from optimizer
[Par2, index] = sortrows(PAR, [5 4 6]);                 %Resort Pareto front values in the order FVU, RMSE, Largest Residual
PAR3 = [index, Par2];                                   %Add the index back into the matrix
PAR4 = unique(PAR3(:,2:7),'rows');                      %Find non-repeating values
X1 = [PAR4(:,1)./100, 10.^PAR4(:,2), 10.^PAR4(:,3)];    %Input variables scaled back to proper numbers
Z1 = PAR4(:,4:6);                                       %Output of optimized   

Lz = size(X1, 1);           %Length of Z1
Z2 = zeros(size(Z1));       %empty Z2

for i = 1:Lz
    [~, Z2(i,:)] = minimizeFlx(X1(i,1),X1(i,2),X1(i,3)); %Find GOF values for all non-repeating pareto front solutions
end


%% Scatter 3D plot for all
sz = 3*ones(Lz, 1);
for i= 1:Lz
    sz(i) = sz(i)+50*i/Lz;  %scale for scatter points to try and see if low fvals correspond to the same pareto front value
end

figure
hold on
scatter3(Z1(:,1),Z1(:,2),Z1(:,3),sz,'filled','DisplayName','Optimized')
scatter3(Z2(:,1),Z2(:,2),Z2(:,3),sz,'filled','DisplayName','Validation')
hold off
xlabel('RMSE')
ylabel('FVU')
zlabel('Residual')
lgd = legend;
view(40,35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';

%% Plot best points 
% Found the best results wrt to optimized and validated data based on vecnorm, FVU, RMSE, and largest residual
%i.e.[M, I] = min(Z(:, ),[],'linear')
BEST = [X1([1 62 103 306],:),Z1([1 62 103 306],:),Z2([1 62 103 306],:)];

%scale points size and color to see if the best of the best is readily
%apparent
sz1 = 10*ones(size(BEST, 1),1);
alpha = size(sz1);
c1 = zeros(size(sz1,1),3);
c2 = zeros(size(sz1,1),3);
for i= 1:length(sz1)
    sz1(i) = sz1(i)+15*(i-1);
    c1(i,:) = [0, 0.54+0.22/3*(i-1), 0.6+0.4/3*(i-1)];
    c2(i,:) = [1-0.1*(i-1), 0, 0];
    alpha(i) = 1-0.15*(i-1); 
end

figure
hold on
scatter3(BEST(:,4),BEST(:,5),BEST(:,6),sz1,c1,'filled','DisplayName','Optimized')
scatter3(BEST(:,7),BEST(:,8),BEST(:,9),sz1,c2,'filled','DisplayName','Validation')
hold off
xlabel('RMSE')
ylabel('FVU')
zlabel('Residual')
lgd = legend;
view(40,35);
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.ZGrid = 'on';

%%
uB = cell(size(BEST,1)+1,1);
uB{1} = u;
vB = cell(size(uB));
vB{1} = v;
bpaB = cell(size(uB));
bpaB{1} = bpa;
for i = 2:length(uB)
    [uB{i},vB{i},bpaB{i}] = minimizeFlx(BEST(i-1,1),BEST(i-1,2),BEST(i-1,3));           % Now pull bpa structures out       
end

%% Plot torque curves, Optimized and validation 
load ForceStrainForFit.mat z
X = linspace(0,620,20); %Pressure for interpolation
X = X(2:20);
Y = linspace(0,1,30);   %Relative strain range for interpolation
str1 = ["Opt1","Best1","Best2","Best3","Best4"];
str2 = ["Optimization"; "Validation"];
for i = 1:2
    %Find old old torque, then plot everything
    clear Yq Xq Vq Fold Fq Mold
    Yq = bpa(i).strain./((bpa(i).rest-bpa(i).Kmax)/bpa(i).rest);
    Xq = bpa(i).P;
    Vq = zeros(size(bpa(i).unitD,1),1);   
        for j = 1:size(bpa(i).unitD, 1)
            if bpa(i).strain(j) >=-0.03 && Yq(j) <=1
                Vq(j) = interp2(X, Y, z, Xq, Yq(j));
%             Vq(j) = f10(Yq(j),Xq);
            elseif Yq(j)>1
            Vq(j) = 0;
            elseif bpa(i).strain(j) < -0.03
            Vq(j) = NaN;
            end
        end

    Fold = Vq.*bpa(i).unitD;    %Force vector
    Fq = (Fold(:,1).^2+Fold(:,2).^2).^(1/2);
    Mold = -bpa(i).mA.*Fq;
    
    figure
    ax = gca;
    hold on
    plot(bpa(i).Ak,bpa(i).M,'.-','DisplayName','Predict original') %"original" is with updated BPA characterization
    plot(bpa(i).Ak,Mold,':','DisplayName','Predict old') %"original" is with updated BPA characterization
    scatter(bpa(i).A_h,bpa(i).M_h,[],'filled','DisplayName','Hybrid')
    scatter(bpa(i).Aexp,bpa(i).Mexp,[],'filled','DisplayName','Experiment')
    for k = 1:length(uB)
        plot(bpaB{k}(i).Ak,bpaB{k}(i).M_p(:,3),'DisplayName',str1{k})
    end
    hold off
    title(str2(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Torque, N\cdotm')
    legend
end

%% Plot muscle length, optimization and validation
for i = 1:2
    figure
    ax = gca;
    Lm = bpa(i).Lmt-2*bpa(i).fitn-bpa(i).ten;      %Original predicted muscle length
    Lm_p = cell(size(uB));
    for kk = 1:length(uB)
     Lm_p{kk} = bpaB{kk}(i).Lmt_p-2*bpaB{kk}(i).fitn-bpaB{kk}(i).ten;    %Optimized muscle length (Lm_p uses optimized values
    end
    hold on
    plot(bpa(i).Ak,Lm,'DisplayName','Predict original')
    scatter(bpa(i).A_h,bpa(i).Lm_h,[],'filled','DisplayName','Measured')
    for k = 1:length(uB)
        plot(bpaB{k}(i).Ak,Lm_p{k},'DisplayName',str1{k})
    end
    hold off
    title(str2(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end

%% Plot moment arm, optimization and validation
for i = 1:2
    figure
    ax = gca;
    G_p = cell(size(uB));
    for kk = 1:length(uB)
        G_p{kk} = (bpaB{kk}(i).mA_p(:,1).^2+bpaB{kk}(i).mA_p(:,2).^2).^(1/2);      %z-axis moment arm for optimized
    end
    hold on
    plot(bpa(i).Ak,bpa(i).mA,'DisplayName','Predict original')
    scatter(bpa(i).A_h,bpa(i).mA_h,[],'filled','DisplayName','Measured')
    for k = 1:length(uB)
        plot(bpaB{k}(i).Ak,G_p{k},'DisplayName',str1{1})
    end
    hold off
    title(str2(i))
    xlabel('\theta_{k}, \circ')
    ylabel('Length, m')
    legend
end
