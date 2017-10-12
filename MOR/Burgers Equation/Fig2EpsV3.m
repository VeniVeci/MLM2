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
    
    ax.YLabel.Position=[-1,leftMarginRatio*(bottomtMarginRatio+hightGraphRatio/2)]
    
    
%     ax.Units='normalized';
%     ax.YLabel.Position=[-0.1,0.5];
    

    % Title Font
    ax.Title.FontSize=FONT_SIZE;
    ax.Title.FontName=FONT_NAME;
    ax.Title.FontWeight='Normal';

 %%

%  folder  = 'C:\YourFolder';
% list    = dir(fullfile(folder, '*.m');
% nFile   = length(list);
% success = false(1, nFile);
% for k = 1:nFile
%   file = list(k).name;
%   try
%     run(fullfile(folder, file));
%     success(k) = true;
%   catch
%     fprintf('failed: %s\n', file);
%   end
% end
%  