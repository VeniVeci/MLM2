%% uvp SIMPLE solver

%Simulating the Velosity and Pressure field 2D equation  
...by the Finite Volume Method using SIMPLE algorithm 
...in a square plane in uniform grid in a "merge grid" n
%    
% in y direction
%     
% Modifications:
% 12-May-2014, WeiX, first edition 

%%
clear
%% -----------------------Setup Parameters---------------------------------
%  -----------------------Problem Parameters---------------------------
Lx=1;                   % length of x axis is 1m
Ly=1;                   % length of y axis is 1m
rho=1;                 % flow density
k=0.01;                    % ???????????

u_T=1;            % Vertical & horizontal velocity of Top    lid                 
u_B=0;            % Vertical & horizontal velocity of Bottom lid
v_L=0;            % Vertical & horizontal velocity of Left   lid   
v_R=0;            % Vertical & horizontal velocity of Right  lid   

% -----------------------Solver Parameters---------------------------
grid_Num_x=81;        % number of volumes on x direction
grid_Num_y=81;        % number of volumes on y direction
MiniConvergenceRate=0.01;
MaxIteration=500;

alpha_p=0.1;

%%
Grid = RetanG_Maker(grid_Num_x,grid_Num_y,Lx,Ly,rho,k,u_T,u_B,v_L,v_R);

p_record=[];
for i=1:MaxIteration
    
    [u,Grid] = RetanG_solveru(Grid);
    [v,Grid] = RetanG_solverv(Grid);
    % Update uv field-----------------------------------------
    Grid.u=u;
    Grid.v=v;
    p_corrector = RetanG_solverp(Grid,'penta');
    AveConvRate=sum(abs((p_corrector./Grid.p)))/(Grid.Num_x_p*Grid.Num_y_p);
    Grid = RetanG_correction(Grid,p_corrector,alpha_p);
    
    if AveConvRate<MiniConvergenceRate
        Grid = RetanG_correction(Grid,p_corrector,alpha_p);
        break;
    end
    % recording ------------
    p_abs=sum(abs(Grid.p(:)));
    p_record=[p_record;p_abs];

    
    
end


%% Visualization ----------------------------------------------
 uplot(1:Grid.Num_x_p,1:Grid.Num_y_p)=0.5*(u(1:Grid.Num_x_u-1,1:Grid.Num_y_u)+u(2:Grid.Num_x_u,1:Grid.Num_y_u));
 vplot(1:Grid.Num_x_p,1:Grid.Num_y_p)=0.5*(v(1:Grid.Num_x_v,1:Grid.Num_y_v-1)+v(1:Grid.Num_x_v,2:Grid.Num_y_v));
for x= 1:Grid.Num_x_p
    for y=1:Grid.Num_y_p
        corx(x,y)=Grid.dx*(x-0.5);
        cory(x,y)=Grid.dy*(y-0.5);
    end
end

% Show number of iterations
iteration=i

figure(1);
quiver(corx,cory,uplot,vplot,1,'k-');
hold on
pcolor(corx,cory,Grid.p),shading interp,
hold off

figure(2)
pcolor(corx,cory,Grid.p),shading interp,
title('Pressure (Steady State)'),xlabel('x'),ylabel('y'),colorbar,axis([0,Lx,0,Ly]),axis equal;
hold on 
% streamline(corx',cory',uplot',vplot',rand(25,1),rand(25,1),[0.1,1000]);
streamline(corx',cory',uplot',vplot',rand(grid_Num_x,1),rand(grid_Num_x,1),[0.1,1000]);
hold off


figure(2); 

% corx=Grid.dx*(0.5):Grid.dx:Grid.dx*(Grid.Num_x_p-0.5);
% cory=Grid.dy*(0.5):Grid.dy:Grid.dy*(Grid.Num_x_p-0.5);

% umin = min(uplot(:));
% umax = max(uplot(:));
% vmin = min(vplot(:));
% vmax = max(vplot(:));
% streamline(dx*0.5:dx:,1:Grid.Num_x_p,uplot,vplot,rand(5,1),rand(5,1));
 
 


