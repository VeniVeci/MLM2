%Figure to .Eps formatting file.
%Formatting figure without adding anything.
%
% Weix 5-Oct-2016


%% Setting paper and graph size (Unit:Inch)
widthInch=10;
hightInch=7;
leftMarginRatio=0.15;               %Left magin ratio
rightMarginRatio=0.10; 
topMarginRatio=0.10;
bottomtMarginRatio=0.15;

FONT_SIZE=28;
FONT_NAME='Times New Roman';

%% Calculating 
widthGrapgRatio=1-leftMarginRatio-rightMarginRatio;
hightGraphRatio=1-topMarginRatio-bottomtMarginRatio;


%% Setting
% Size & Position
set(gcf,'Units','inches','Position',[2,3,widthInch,hightInch]); %Golden ratio 1.618
set(gca,'Units','inches','Position',[widthInch*leftMarginRatio,hightInch*bottomtMarginRatio,widthInch*widthGrapgRatio,hightInch*hightGraphRatio]);

% Axis Font
ax = gca; 
ax.FontSize=FONT_SIZE;
ax.FontName=FONT_NAME;
ax.FontWeight='Normal';


    ax.XLabel.FontSize=FONT_SIZE; %Or 'italic'
    ax.XLabel.FontName=FONT_NAME;
    ax.XLabel.HorizontalAlignment='center';
    ax.XLabel.VerticalAlignment='Top';
    
    
    % ax.XLabel.FontAngle='italic';

    ax.YLabel.FontSize=FONT_SIZE; %Or 'italic'
    ax.YLabel.FontName=FONT_NAME;
    ax.YLabel.HorizontalAlignment='center';
    ax.YLabel.VerticalAlignment='bottom';
    % ax.YLabel.FontAngle='italic';

    % Title Font
    ax.Title.FontSize=FONT_SIZE;
    ax.Title.FontName=FONT_NAME;
    ax.Title.FontWeight='Normal';

 %% Set xy limit and steps
    ax = gca; 
    ax.YScale='log';
    %Range for log scale
    ax = gca; 
    highBound=0;
    lowBound=-8;
    yTick=logspace(lowBound,highBound,highBound-lowBound+1);
    ax.YTick=yTick(1:1:end);   
    ax.YLim=[10^lowBound,10^highBound];
    
    
%% label position
% ax.YLabel.Units='Inch';
% ax.YLabel.Position=[-0.8,2];
% 
% ax.XLabel.Units='Inch';
% ax.XLabel.Position=[3,-0.5];

ff=gcf;
ff.PaperPositionMode='auto'     %!important! ensure print size = visual size !!!

