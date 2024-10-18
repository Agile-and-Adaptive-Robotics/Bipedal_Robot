%compare values of optimized and validation from 10mm Flx_pinned pareto
%front

PAR = [Pareto_front, Pareto_Fvals];
[Par2, index] = sortrows(PAR, [5 4 6]);
PAR3 = [index, Par2];
PAR4 = unique(PAR3(:,2:7),'rows');
X1 = [PAR4(:,1)./100, 10.^PAR4(:,2), 10.^PAR4(:,3)];
Z1 = PAR4(:,4:6);

Lz = size(X1, 1);
Z2 = zeros(size(Z1));
g = zeros(size(Z1));

for i = 1:Lz
    [~, g(i,:)] = minimizeFlx(X1(i,1),X1(i,2),X1(i,3));
    Z2(i,:)= g(i,:);
end

sz = 3*ones(Lz, 1);
for i= 1:Lz
    sz(i) = sz(i)+50*i/Lz;
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