function [] = FigAutoFormate()
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

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

    %% setting formate
   


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
    ax.XLabel.FontAngle='italic';
%     ax.XLabel.FontAngle='normal';

    ax.YLabel.FontSize=FONT_SIZE; %Or 'italic'
    ax.YLabel.FontName=FONT_NAME;
    ax.YLabel.HorizontalAlignment='center';
    ax.YLabel.VerticalAlignment='middle';    
    ax.YLabel.FontAngle='italic';
%     ax.YLabel.FontAngle='normal';   
    
    ax.Title.FontSize=FONT_SIZE;
    ax.Title.FontName=FONT_NAME;
    ax.Title.FontWeight='Normal';
    

    ax.Title.FontSize=FONT_SIZE;

%     ax.Units='Inch';d 
%      ax.Position=[1.2, 1.0, 6.25, 4.25];
    %  ax.Position=[0.1, 0.16, 0.8, 0.75];

    % ax.YLabel.Position=[-1,1];
%     ax.YLabel.Units='Inch';
%     ax.YLabel.Position=[-0.8,2.2];

%     set(gcf,'Units','inches','Position',[2,3,8.042,5.75]); %Golden ratio 1.618
%     set(gca,'Units','inches','Position',[1.2, 1.0, 6.25, 4.25]);
    

    set(gcf,'Units','inches','Position',[2,3,8.042,5.75]); %Golden ratio 1.618
    set(gca,'Units','inches','Position',[1.2, 1.0, 6.25, 4.25]);
    
%     ax = gca; 
%     ax.YLabel.HorizontalAlignment='center';
    
    ax.YLabel.Units='Inch';
    ax.YLabel.Position=[-0.8,2];

    NUMBER_EXTENSION_CHAR=4;
%     saveas(gcf,figList(i).name(1:end-NUMBER_EXTENSION_CHAR))

%     set(gcf,'Renderer','OpenGL')
    print(figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'-depsc')
    
%     print(figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'-deps','-zbuffer')
    
    %-Change renderer.  Normally useed when BUG or large output result.
    % set(gcf,'Renderer','OpenGL')
    % set(gcf,'Renderer','Painters')
    % set(gcf,'Renderer','zbuffer')
%     
%     print(figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'-depsc','-painters')
%     
    close(iFig);
    
    

end




end

