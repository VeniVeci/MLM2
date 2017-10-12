%Figure to .Eps formatting file.
%Just run the script
%
% Weix 5-Oct-2016


%% Setting paper and graph size (Unit:Inch)
xInch=10;
yInch=7;
leftMarginRatio=0.15;   %Left magin ratio
rightMarginRatio=0.10; 
topMarginRatio=0.15;
bottomtMarginRatio=0.15;



FONT_SIZE=30;

graphRatio=(1-leftMarginRatio-rightMarginRatio)*xInch;


widthCenter=xInch*(marginRatio+0.5*graphRatio);
hightCenter=yInch*(marginRatio+0.5*graphRatio);

%Paper size
set(gcf,'Units','inches','Position',[2,3,xInch,yInch]); %Golden ratio 1.618
ax = gca; 
ax.Units='inches';
% ax.Position=[FONT_SIZE/20,FONT_SIZE/20,highInch-2*FONT_SIZE/20,widthInch-2*FONT_SIZE/20];
ax.Position=[xInch*marginRatio,yInch*marginRatio,xInch*graphRatio,yInch*graphRatio];

ax.FontSize=FONT_SIZE;
% ax.Position=[0.1,0.13,0.8,0.8];
ax.FontName='Times New Roman';
ax.FontWeight='Normal';
% ax.OuterPosition=[0,0,1,1];


% 
% title('(b)','FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal')
% 
% 
% xlabel({'Approximate manifold dimension'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
% % xlabel({'POD basis dimension'},'FontSize',30,'FontName','Times New Roman','FontWeight','normal');
% ax.XLabel.Units='inches';
% ax.XLabel.Position=[5,-0.6,0];
% 
% ylabel({'Square error'},'FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal');
% ax.YLabel.Units='inches';
% ax.YLabel.Position=[-0.8,0,0];
% 
% 
% title('(a)','FontSize',FONT_SIZE,'FontName','Times New Roman','FontWeight','normal')
