%% Us_HeatDC_MOR_SolveFunc_Demo

clear 

Paras.Lx=1;                % length of x axis is 1m
Paras.Ly=1;                % length of y axis is 1m
Paras.k=1;               % Thermal conductivity is 1000W/m/K
Paras.rho=1;                 % flow density

% Paras.T_L=100;               % Temperature of Left    boundary 
% Paras.T_R=100;               % Temperature of Right   boundary 
% Paras.T_B=400;               % Temperature of Bottom  boundary 
% Paras.T_T=400;               % Temperature of Top     boundary 

Paras.T_L=0;               % Temperature of Left    boundary 
Paras.T_R=0;               % Temperature of Right   boundary 
Paras.T_B=0;               % Temperature of Bottom  boundary 
Paras.T_T=0;               % Temperature of Top     boundary 

% u_x=1;                 % Velocity x direction 
% u_y=-1;                 % Velocity y direction 


% -----------------------Solver Parameters---------------------------
Paras.Num_node_x=20;        % number of volumes on x direction
Paras.Num_node_y=20;        % number of volumes on y direction


uParas.ampx=1;
uParas.ampy=1;
uParas.x0=0.5;
uParas.y0=0.5;

Paras.T_0=1000;              % initial temperature

Paras.t_end=0.2;             % total simulation time 
Paras.dt=0.002;                  % time step


approximate_degree=5;


% ----------Generate Initial field-------------

[T,X,Y] = GausFieldF(Paras);
X_int=T(:);

%--------------Main--------------------------

[T_Rec,Time]=Us_HeatDC_SolveFunc(Paras,uParas,X_int);


[U,S,~]=svd(T_Rec);  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
U=U(:,1:approximate_degree);
eigenvalues=diag(S);

[T_Rec_mor,Time_mor]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U);


%-------------------Plot----------------------------------------------
dx=Paras.Lx/Paras.Num_node_x;
dy=Paras.Ly/Paras.Num_node_y;
[Xcor,Ycor] = meshgrid(0+dx/2:dx:Paras.Lx-dx/2,0+dy/2:dy:Paras.Lx-dy/2);    

Index_frame=1;
for i=1:Paras.t_end/Paras.dt+1
    
%     Tx=reshape(P(1,:),Num_node_y,Num_node_x);       % Rearrange coordinate
%     Ty=reshape(P(2,:),Num_node_y,Num_node_x);

    T =reshape(T_Rec(:,i),Paras.Num_node_y,Paras.Num_node_x);
    T2=reshape(T_Rec_mor(:,i),Paras.Num_node_y,Paras.Num_node_x);

    figure(1)
    title('HDM')
    subplot(2,1,1),contour(Xcor,Ycor,T),
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(2,1,2),pcolor(Xcor,Ycor,T),shading interp,
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    
    
    figure(2)
    title('Rom')
    subplot(2,1,1),contour(Xcor,Ycor,T2),
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal
    subplot(2,1,2),pcolor(Xcor,Ycor,T2),shading interp,
    title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,%axis([0,Lx,0,Ly]),%axis equal

%     figure(2)
%     pcolor(Tx,Ty,T),shading interp,
%     title('Temperature (Steady State)'),xlabel('x'),ylabel('y'),colorbar,axis([0,Lx,0,Ly]),axis equal;
% 
%     % figure(3)
%     % mesh(Tx,Ty,T),shading interp,
%     figure(4)
%     surf(Tx,Ty,T),shading interp,
%     figure(5)
%     contour(Tx,Ty,T)


     F(Index_frame)=getframe;
     Index_frame=Index_frame+1;

end
