


set(gcf,'Units','inches','Position',[2,3,10,8.18]); %Golden ratio 1.618

ax = gca; 
ax.FontSize=30;
% ax.Position=[0.1,0.13,0.8,0.8];
ax.FontSize=30;
ax.FontName='Times New Roman';
ax.FontWeight='Normal';
% ax.OuterPosition=[0,0,1,1];

ax.Units='inches';
ax.Position=[1.2,1.2,8,6.4];

xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
ax.XLabel.Units='inches';
ax.XLabel.Position=[4,-0.6,0];

ylabel({'Square error'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
ax.YLabel.Units='inches';
ax.YLabel.Position=[-0.8,3,0];

title('(b)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')

%% 
% ax.OuterPosition=[0,0.0,1,1];
ax.OuterPosition=[0,0.05,0.9,0.95];



ylabel({'Square error'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
title('(a)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')




ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=4;
dB=0;
ax.YTick=logspace(dB,uB,uB-dB+1);
ax.YLim=[10^dB,10^uB];



%% For Profile plot
% ylabel({'x_2'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
% xlabel({'x_1'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');