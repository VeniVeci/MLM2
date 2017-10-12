function [] = FigAutoFormateSize()
%% This script tries to locat figure in current folder and export them as eps formate.
%  For font size 
% 08-01-2017, WeiX, first edition 

%% Fetch all figure in folder
figList= dir('*.fig');
nFig=length(figList);

%--use reference figure

if length(dir('Reference.fig'))==0; error('No reference figure is found'); end 
refFig=openfig('Reference.fig');

%%
for i = 1:nFig
    iFig=openfig(figList(i).name);

    %% setting formate

    
    
    iFig.Units=refFig.Units;
    iFig.Position=refFig.Position;
    
    iFig.CurrentAxes.Units=refFig.CurrentAxes.Units;
    iFig.CurrentAxes.Position=refFig.CurrentAxes.Position;
    iFig.CurrentAxes.FontName=refFig.CurrentAxes.FontName;
    iFig.CurrentAxes.FontSize=refFig.CurrentAxes.FontSize;
    iFig.CurrentAxes.FontWeight=refFig.CurrentAxes.FontWeight;
    iFig.CurrentAxes.FontAngle=refFig.CurrentAxes.FontAngle;
  
    iFig.CurrentAxes.XLabel=refFig.CurrentAxes.XLabel;
    iFig.CurrentAxes.YLabel=refFig.CurrentAxes.YLabel;
    
    iFig.CurrentAxes.Title=refFig.CurrentAxes.Title;
    
    NUMBER_EXTENSION_CHAR=4;
%     saveas(gcf,figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'epsc')
    print(figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'-depsc')
    %-Change renderer.  Normally useed when BUG or large output result.
    % set(gcf,'Renderer','OpenGL')
    % set(gcf,'Renderer','Painters')
    % set(gcf,'Renderer','zbuffer') %This renderer is useful when the
    % generated image is too large.
%     
%     print(figList(i).name(1:end-NUMBER_EXTENSION_CHAR),'-depsc','-painters')
%     
    close(iFig);
    
end
    close(refFig);



end

