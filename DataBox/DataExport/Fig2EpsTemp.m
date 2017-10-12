% Fig2EpsTemp
% A template to follow to output MATLAB figure to eps for Latex. 


%% METHOD I: DIRECT SAVE AS .EPS
%Used with Figure properties and directly use save as eps rather than 'print' command
%Usage:
%1.On figure click Edit/Figure properties to enter editor mode
%2.Set axis position for all using: ax = gca;ax.Position=[0.05,0.13,0.8,0.8]; Shown followed 
%3.Adjust the window to maximize figure area occupation 
%4.Copy and past following command. Adjust value to needs.
%(4.5). Before next step one can use adjust the window size to ensure unity of all figure
%5.Up to the point all figure should have same scale/position/ratio. Then use File/Save as into .eps file. //


%change Axis(Font, whole graph position)
ax = gca; 
ax.FontSize=28;
ax.Position=[0.1,0.13,0.8,0.8];

%YScale
ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=4;
dB=-6;
ax.YTick=logspace(dB,uB,uB-dB+1);
ax.YLim=[10^dB,10^uB];


%Label
ylabel({'SSE'},'FontSize',24,'FontWeight','bold');
xlabel({'Dimensions'},'FontSize',24,'FontWeight','bold');

ylabel({'Relative error'},'FontSize',24,'FontWeight','bold');
xlabel({'Dimensions'},'FontSize',24,'FontWeight','bold');

ylabel({'\chi m^{-1}'},'FontSize',24,'FontWeight','bold');
xlabel({'\xi m^{-1}'},'FontSize',24,'FontWeight','bold');


%Legend (Font, position)
leg = findobj(gcf,'Tag','legend');
leg.FontSize=20;
leg.Position=[0.65,0.68,0.2,0.2];


%color bar range
caxis([-0.007, 0.007])

%Adjust windows size 
hFig = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 800 600];   %Set position and size 



%% METHOD II: PRINT COMMAND (Not Finished)
% This method is said to generate higher quality eps but more difficult to handle.
set(0,'DefaultFigureWindowStyle','normal')

%Set paper size/position
set(gcf,'units','centimeters')
set(gcf,'PaperPositionMode', 'manual')
set(gcf,'papersize',[5,5])
set(gcf,'paperposition',[0,0,5,5])
set(gcf, 'renderer', 'painters');

%Set paper size/position with easier command
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 5];
fig.PaperPositionMode = 'manual';


legend ('FontSize',24)
% xlabel(?Time(s)?,...
% ?FontUnits?,?points?,...
% ?FontWeight?,?normal?,...
% ?FontSize?,7,...
% ?FontName?,?Times?)

%
%Output eps
print -depsc myplot.eps
