%% Unsteady state heat diffusion convection 2D problem solved by FVM. Test using various input parameter (generate by sobol sequence) and test.


clear
%% -----------------------Setup Parameters----------------------------------
%  -----------------------Problem Parameters---------------------------
Paras.Lx=10;                % length of x axis is 1m
Paras.Ly=10;                % length of y axis is 1m
Paras.k=10;                 % Thermal conductivity is 1000W/m/K     % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)
Paras.rho=1;                % flow density       % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)
Paras.T_L=0;                % Temperature of Left    boundary 
Paras.T_R=0;                % Temperature of Right   boundary 
Paras.T_B=0;                % Temperature of Bottom  boundary 
Paras.T_T=0;                % Temperature of Top     boundary 
Paras.u_x=2;               % Velocity x direction 
Paras.u_y=-20;              % Velocity y direction 

Paras.T_0=0;                % initial temperature

Paras.t_end=0.5;            % total simulation time 
Paras.dt=0.01;              % time step

% -----------------------Solver Parameters---------------------------
Paras.Num_node_x=20;        % number of volumes on x direction
Paras.Num_node_y=20;        % number of volumes on y direction

% --------------------Complex Boundary condition parameters-------------
Paras.b=100;
Paras.a=10;

% --------------------Sobol Sequence the parameters---------------------
Num_X=200;
rang=[20,20,10;-20,-20,0]; 
[X] = Sobo_Generator(Num_X,rang);
Index_snapshot=1:5:51;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

for i = 9:Num_X
    Paras.u_x=X(i,1);
    Paras.u_y=X(i,2);
    Paras.a=X(i,3);
    
    [T,time]=Us_HeatDC_2D_FVM(Paras);
    T=T(:,Index_snapshot);
    Y(i,:)=T(:)';
    Time(i,1)=time;
end

%   save('ExpData.mat','X','Y','Paras');






