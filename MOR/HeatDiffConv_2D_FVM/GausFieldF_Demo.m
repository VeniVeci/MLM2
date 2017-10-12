%% GausFieldF_Demo


Paras.Lx=1;
Paras.Ly=1;

Paras.Num_node_x=50;
Paras.Num_node_y=50;

Options=[];

[T,X,Y] = GausFieldF(Paras);

%% 
figure
surf(X,Y,T);
caxis([min(T(:))-.5*range(T(:)),max(T(:))]);
% axis([-3 3 -3 3 0 .4])
xlabel('x'); ylabel('y'); zlabel('Heat temperature');

ax = gca; 
ax.FontSize=20;
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'\xi '},'FontSize',24,'FontWeight','bold');
xlabel({'\eta '},'FontSize',24,'FontWeight','bold');
zlabel({'contaminant concentration'},'FontSize',24,'FontWeight','bold');

set(ax, 'xlim', [0 1])
set(ax, 'ylim', [0 1])
set(ax, 'zlim', [0 50])

ax.XTick = [0:0.2:1];
ax.YTick = [0:0.2:1];
ax.ZTick = [0:10:50];

%Adjust windows size 
h = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 700 700];   %Set position and size 

%%-------------
figure
pcolor(X,Y,T),shading interp;

colorbar

txt = 'Gaussian 1';
text(0.1,0.3,txt,'Color','Red','FontSize',20)

txt = 'Gaussian 2';
text(0.1,0.9,txt,'Color','Red','FontSize',20)

txt = 'Gaussian 3';
text(0.7,0.9,txt,'Color','Red','FontSize',20)

hold on
scatter(0.8,0.2,100,'filled','d','r')
txt = 'Sensor';
text(0.75,0.3,txt,'Color','Red','FontSize',20)
hold off

hold on
scatter(0.5,0.5,100,'filled','d','r')
txt = 'Sensor';
text(0.45,0.6,txt,'Color','Red','FontSize',20)
hold off


ax = gca; 
ax.FontSize=20;
ax.Position=[0.1,0.13,0.8,0.8];

ylabel({'\xi '},'FontSize',24,'FontWeight','bold');
xlabel({'\eta '},'FontSize',24,'FontWeight','bold');

set(ax, 'xlim', [0 1])
set(ax, 'ylim', [0 1])

ax.XTick = [0:0.2:1];
ax.YTick = [0:0.2:1];

% caxis([0 50])      %color bar
hcb=colorbar;
set(hcb, 'ylim', [0 50])
set(hcb,'YTick',[0:10:50])


%Adjust windows size 
h = gcf;
h.WindowStyle='normal';     %This would release figure from editor mode
h.Position=[500 500 700 700];   %Set position and size 