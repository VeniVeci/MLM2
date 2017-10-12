%FigExport2


figure
nPod=[5:5:60];

% boxplot(SSE_dxdt_U_GlobalSnapS(:,nPod),nPod)
boxplot(SSE_dxdt_U_GPESnapS(:,nPod),nPod)

title(sprintf('L^2 Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i', DrMethod,Num_Train,Num_Test,Num_Snapshot))
set(gca,'yscale','log');


ax = gca; 
ax.FontSize=30;
ax.FontName='Times New Roman';
ax.FontWeight='normal';
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'Relative error'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
title('(b)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')
% 
ax.YLabel.Position=[-1.2,0.01,-1];

% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position')- [0.2 0])

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=3;
dB=-5;
ax.YTick=logspace(dB,uB,uB-dB+1);
ax.YLim=[10^dB,10^uB];

% %Adjust windows size 
% hFig = gcf;
% h.WindowStyle='normal';     %This would release figure from editor mode


ax = gca; 
ax.FontSize=30;
ax.FontName='Times New Roman';
ax.FontWeight='normal';
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'x_2'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
xlabel({'x_1'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');


fig = gcf; 
set(gcf,'PaperPositionMode','auto')
print(fig,'MySavedPlot','-depsc')
