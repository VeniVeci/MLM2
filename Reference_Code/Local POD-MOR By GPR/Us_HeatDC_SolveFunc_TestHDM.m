%% Us_HeatDC_SolveFunc_TestHDM

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

Paras.T_0=1000;              % initial temperature

Paras.t_end=0.2;             % total simulation time 
Paras.dt=0.005;                  % time step

% -----------------------Solver Parameters---------------------------
Paras.Num_node_x=50;        % number of volumes on x direction
Paras.Num_node_y=50;        % number of volumes on y direction


IntParas.ampx=1;
IntParas.ampy=1;
IntParas.x0=0.5;
IntParas.y0=0.7;




% --------------------Sobol Sequence the parameters---------------------
Num_X=500;
rang=[20,20,10;-20,-20,0]; 
[X] = Sobo_Generator(Num_X,rang);
Index_snapshot=1:5:51;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

for i = 1:Num_X
    Paras.u_x=X(i,1);
    Paras.u_y=X(i,2);
    Paras.a=X(i,3);
    
    [T,time]=Us_HeatDC_2D_FVM(Paras);
    [T,Time]=Us_HeatDC_SolveFunc(Paras,IntParas);
    
    T=T(:,Index_snapshot);
    Y(i,:)=T(:)';
    Time(i,1)=time;
    
    
    
end

%   save('ExpData.mat','X','Y','Paras');




