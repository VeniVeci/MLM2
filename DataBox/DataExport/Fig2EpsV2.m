
% ax = gca; 
% xlabel({'POD basis dimension'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
% ax.XLabel.Units='inches';
% ax.XLabel.Position=[3,-0.6,0];
% 
% ylabel({'Square error'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
% ax.YLabel.Units='inches';
% ax.YLabel.Position=[-0.8,1.8,0];

% xlabel({'x'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
% ylabel({'u(x,t)'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
% title('(c)','FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal')


FONT_SIZE=28;
FONT_NAME='Times New Roman';

ax = gca; 
ax.FontSize=FONT_SIZE;


ax.FontSize=FONT_SIZE;
ax.FontName=FONT_NAME;
ax.FontWeight='Normal';


    ax.XLabel.FontSize=FONT_SIZE; %Or 'italic'
    ax.XLabel.FontName=FONT_NAME;
    ax.XLabel.HorizontalAlignment='center';
    ax.XLabel.VerticalAlignment='Top';
%     ax.XLabel.FontAngle='italic';
    ax.XLabel.FontAngle='normal';

    ax.YLabel.FontSize=FONT_SIZE; %Or 'italic'
    ax.YLabel.FontName=FONT_NAME;
    ax.YLabel.HorizontalAlignment='center';
    ax.YLabel.VerticalAlignment='bottom';    
%     ax.YLabel.FontAngle='italic';
    ax.YLabel.FontAngle='normal';   
    
    ax.Title.FontSize=FONT_SIZE;
    ax.Title.FontName=FONT_NAME;
    ax.Title.FontWeight='Normal';
    
    


ax.XLabel.FontSize=FONT_SIZE; %Or 'italic'
ax.XLabel.HorizontalAlignment='center';
ax.XLabel.VerticalAlignment='Top';
% ax.XLabel.FontAngle='italic';

ax.YLabel.FontSize=FONT_SIZE; %Or 'italic'
ax.YLabel.HorizontalAlignment='center';
ax.YLabel.VerticalAlignment='bottom';
% ax.YLabel.FontAngle='italic';

ax.Title.FontSize=FONT_SIZE;

ax.Units='Inch';
 ax.Position=[1.2, 1.0, 6.25, 4.25];
%  ax.Position=[0.1, 0.16, 0.8, 0.75];

% ax.YLabel.Position=[-1,1];
ax.YLabel.Units='Inch';
ax.YLabel.Position=[-0.8,2];

ax.XLabel.Units='Inch';
ax.XLabel.Position=[1.5,2];


set(gcf,'Units','inches','Position',[2,3,8.042,5.75]); %Golden ratio 1.618
set(gca,'Units','inches','Position',[1.2, 1.0, 6.25, 4.25]);


%% set axis
    ax = gca; 
    ax.YScale='log';
    %Range for log scale
    ax = gca; 
    highBound=0;
    lowBound=-8;
    yTick=logspace(lowBound,highBound,highBound-lowBound+1);
    ax.YTick=yTick(1:2:end);   
    ax.YLim=[10^lowBound,10^highBound];




% print(gcf.name(1:end-4),'-depsc')




%%
% highInch=10;
% widthInch=8;
% 
% FONT_SIZE=27;
% 
% set(gcf,'Units','inches','Position',[2,3,8,6]); %Golden ratio 1.618
% 
% ax = gca; 
% ax.Units='inches';
% ax.Position=[1.2,1.2,6.25,4.25];
% 
% 
% ax.FontSize=FONT_SIZE;
% % ax.Position=[0.1,0.13,0.8,0.8];
% ax.FontSize=FONT_SIZE;
% ax.FontName='Times New Roman';
% ax.FontWeight='Normal';
% 




% ax.OuterPosition=[0,0,1,1];

    % xlabel({'Approximate manifold dimension'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
%     xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
%     ax.XLabel.Units='inches';
%     ax.XLabel.Position=[3,-0.6,0];
% 
%     ylabel({'Square error'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
%     ax.YLabel.Units='inches';
%     ax.YLabel.Position=[-0.8,1.8,0];






% 
% 
% title('(b)','FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal')
% 
% 
% xlabel({'Approximate manifold dimension'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
% % xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% ax.XLabel.Units='inches';
% ax.XLabel.Position=[2.5,-0.6,0];
% 
% ylabel({'Square error'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
% ax.YLabel.Units='inches';
% ax.YLabel.Position=[-0.8,2,0];
% 
% 
% title('(a)','FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal')









% 
% ax = gca; 
% ax.XTick=[1:1:30];
% ax.XLabel=[4:4:28];
% 
% 
% 
% 
% xlabel({'x'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
% ax.XLabel.Units='inches';
% ax.XLabel.Position=[3,-0.6,0];
% 
% ylabel({'u(x,t)'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
% ax.YLabel.Units='inches';
% ax.YLabel.Position=[-0.8,2,0];
% 
% title('(a)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')
% 
% %%
% 
% 
% 
% 
% 
% 
% xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% ax.XLabel.Units='inches';
% ax.XLabel.Position=[3,-0.6,0];
% 
% ylabel({'Square error'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% ax.YLabel.Units='inches';
% ax.YLabel.Position=[-0.8,2,0];
% 
% title('(b)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')
% 
% %% 
% % ax.OuterPosition=[0,0.0,1,1];
% ax.OuterPosition=[0,0.05,0.9,0.95];
% 
% 
% 
% ylabel({'Square error'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% xlabel({'Approximate manifold dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% % xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% title('(a)','FontSize',30,'FontName','Times New Roman','FontWeight','normal')
% 
% 
% 
% 
% ax = gca; 
% ax.YScale='log';
% %Range for log scale
% ax = gca; 
% uB=4;
% dB=0;
% ax.YTick=logspace(dB,uB,uB-dB+1);
% ax.YLim=[10^dB,10^uB];
% 


%% For Profile plot
% ylabel({'x_2'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');
% xlabel({'x_1'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal','FontAngle','italic');