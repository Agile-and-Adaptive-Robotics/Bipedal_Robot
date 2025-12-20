openfig('Fmax_3D.fig')
ax1=gca;
openfig('Fnorm.fig');
ax2=gca;

 
figure;
tc1=tiledlayout(1,2);
ax1.Parent=tc1;
ax1.Layout.Tile=1;
ax2.Parent=tc1;
ax2.Layout.Tile=2;