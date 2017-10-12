%% Unsteady state heat diffusion convection 2D problem solved by FVM. Test using various input parameter (generate by sobol sequence) and test.


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

Paras.T_0=1000;              % initial temperature

Paras.t_end=0.1;             % total simulation time 
Paras.dt=0.001;                  % time step

% -----------------------Solver Parameters---------------------------
Paras.Num_node_x=20;        % number of volumes on x direction
Paras.Num_node_y=20;        % number of volumes on y direction

IntParas.ampx=1;
IntParas.ampy=1;
IntParas.x0=0.5;
IntParas.y0=0.7;


% --------------------Sobol Sequence the parameters---------------------
Num_X=5;
% rang=[20,20,10;-20,-20,0]; 
% [X] = Sobo_Generator(Num_X,rang);

rang=[1,1;0,0]; 
[X] = Sobo_Generator(Num_X,rang);

% Index_snapshot=1:5:51;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;



h = waitbar(0,'TestHDM is working very hard for you');
for i = 1:Num_X
    
    IntParas.x0=X(i,1);
    IntParas.y0=X(i,2);
   
    [y,Time]=Us_HeatDC_SolveFunc(Paras,IntParas);
%     [Rec_X,Time]=Us_HeatDC_SolveFunc(Paras,IntParas)
    
    Y_Rec(:,:,i)=y; 
    Time_Rec(i,1)=Time;
        
%     y=y(:,Index_snapshot);
%     Y(i,:)=y(:)';
   
    waitbar(i/Num_X);
%     waitbar(i/Num_X,'TestHDM is working very hard for you');
    
%     h=waitbar(i/Num_X,'TestHDM is working very hard for you');
%     h=waitbar(i/Num_X,sprintf('TestHDM is working very hard for you. %d%% done...',i/Num_X));
    
end
close(h);

%     save('ExpData11.mat','X','Y','Paras','rang');
    save('HeatDC_DBC_ExpData1.mat','X','Y_Rec','Time','Paras','rang');





