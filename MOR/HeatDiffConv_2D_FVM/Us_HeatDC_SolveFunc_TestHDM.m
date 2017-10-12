%% Us_HeatDC_SolveFunc_TestHDM

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

Paras.t_end=0.2;             % total simulation time 
Paras.dt=0.002;                  % time step

% -----------------------Solver Parameters---------------------------
Paras.Num_node_x=20;        % number of volumes on x direction
Paras.Num_node_y=20;        % number of volumes on y direction

uParas.ampx=1;
uParas.ampy=1;
uParas.x0=0.5;
uParas.y0=0.7;

% ----------Generate Initial field-------------

[T,X,Y] = GausFieldF(Paras);
X_int=T(:);

% --------------------Sobol Sequence the parameters---------------------
Num_X=500;
% rang=[20,20,10;-20,-20,0]; 
% [X] = Sobo_Generator(Num_X,rang);

rang=[1,1;0,0]; 
[X] = Sobo_Generator(Num_X,rang);
% Index_snapshot=1:5:51;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;



h = waitbar(0,'TestHDM is working very hard for you');
for i = 1:Num_X
    
    uParas.x0=X(i,1);
    uParas.y0=X(i,2);
   
    [y,Time]=Us_HeatDC_SolveFunc(Paras,uParas,X_int);
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

% Exam generated data
yy=reshape(Y_Rec,[],Num_X); %500 for this dataset
yy=yy';
syy=sum(yy,2);
med=median(syy);
med_abs=abs(med);

alpha=10;
index=find(syy<(alpha*med_abs+med) & syy>(-alpha*med_abs+med));

X=X(index,:);
Y_Rec=Y_Rec(:,:,index);
Time_Rec=Time_Rec(index,:);


%     save('ExpData11.mat','X','Y','Paras','rang');
%     save('HeatDC_DBC_ExpData21.mat','X','Y_Rec','Time_Rec','Paras','uParas','rang','X_int');





