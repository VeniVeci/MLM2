% Fig2Eps_Boxplot_Demo
%
% It uses boxplot generate by GPE directly.
%
% Modifications:
% 12-April-2016, WeiX, first edition 

ax = gca; 
ax.FontSize=30;


ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=6;
dB=1;
tick=logspace(dB,uB,uB-dB+1);

ax.YTick=tick(1:1:6);
ax.YLim=[10^dB,10^uB];

ax = gca; 
ax.XTick=1:2:15;
ax.XTickLabel=1:2:15;


ylabel({'Square error'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');



ax.Position=[0.1,0.13,0.8,0.8];

% Set(ax,yla

ylabel(sprintf('Square error'),'FontSize',24,'FontWeight','bold');
xlabel(sprintf('POD Basis'),'FontSize',24,'FontWeight','bold');
% ylabel(sprintf('Relative error (x10^{-2})'),'FontSize',24,'FontWeight','bold');
% ylabel('10_2','FontSize',24,'FontWeight','bold');
% ylabel('x_{2}');

% set(ax, 'xlim', [0 1])
set(ax, 'ylim', [0.0001 1])

% ax.YTick = [0:0.01:0.03];

hFig = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 800 600];   %Set position and size 



