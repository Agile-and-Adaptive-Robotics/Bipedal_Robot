openfig('MaxForce10.fig')
ax1=gca;
openfig('MaxForce20.fig');
ax2=gca;
openfig('MaxForce10_all.fig');
ax3=gca;
openfig('kmax10.fig')
ax4=gca;
 
figure;
tc1=tiledlayout(2,1);
ax1.Parent=tc1;
ax1.Layout.Tile=1;
ax2.Parent=tc1;
ax2.Layout.Tile=2;

figure;
tc2=tiledlayout(2,1);
ax3.Parent=tc2;
ax3.Layout.Tile=3;
ax4.Parent=tc2;
ax4.Layout.Tile=4;