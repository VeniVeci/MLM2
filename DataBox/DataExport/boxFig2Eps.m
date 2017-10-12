function [] = boxFig2Eps()
%convert all box plot .fig file to .eps file.
%
% Wei Xing, 22-12-2015, first edition 


%% Fetch all figure in folder
figList= dir('*.fig');
nFig=length(figList);

%--use reference figure
% refFig=dir('Reference.fig');
% if length(refFig)==0; error('No reference figure is found'); end 

% use specify parameters
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

%%
for i = 1:nFig
    iFig=openfig(figList(i).name);

    %% get figure 
    ax = gca; 
    
    %% overall font
    ax.FontSize=FONT_SIZE;
    ax.FontName=FONT_NAME;
    ax.FontWeight='Normal';
    
    %% title font
    ax.Title.FontSize=FONT_SIZE;
    ax.Title.FontName=FONT_NAME;
    ax.Title.FontWeight='Normal';

    %% Axis font
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
    

%     ax.Units='Inch';
%      ax.Position=[1.2, 1.0, 6.25, 4.25];
    %  ax.Position=[0.1, 0.16, 0.8, 0.75];

    % ax.YLabel.Position=[-1,1];
%     ax.YLabel.Units='Inch';
%     ax.YLabel.Position=[-0.8,2];

%     ax.XLabel.Units='Inch';
%     ax.XLabel.Position=[1.5,2];

    
    %% set frame
    figHandle=gcf;
    figHandle.WindowStyle='normal';
    
    set(gca,'Units','inches','Position',[1.4, 1.0, 6.25, 4.25]);
    
    set(gcf,'Units','inches','Position',[2,3,8.042,5.75]); 

    


%% Extra:    
%     % set axis

%     ax = gca; 
%     ax.YScale='log';
%     %Range for log scale
%     ax = gca; 
%     highBound=-2;
%     lowBound=-3;
%     yTick=logspace(lowBound,highBound,highBound-lowBound+1);
%     ax.YTick=yTick(1:1:end);   
%     ax.YLim=[10^lowBound,10^highBound];

    
    %% Export
    NUMBER_EXTENSION_CHAR=4;    %.fig equals 4;
    
    figHandle.PaperPositionMode='auto'     %!important! ensure print size = visual size !!!
%     print(figure(1),'123','-depsc')
    print(figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'-depsc')
    
    close(iFig);
    
    

end




end

