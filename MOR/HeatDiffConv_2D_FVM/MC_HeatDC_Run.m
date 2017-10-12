%% MC_HeatDC_Run
% Monte Carlo process
%
% Modifications:
% WeiX, 9-4-2016, Create.



clear

% ----------------------Sampling-------------------
% rand(10,1)
num_sample=3000;
mu = [0.5,0.5];
sigma = [0.01,0;0,0.01];
% rng default  % For reproducibility
rng(1)
X_sample = mvnrnd(mu,sigma,num_sample);

% figure
% plot(X_sample(:,1),X_sample(:,2),'+')
figure
plotmatrix(X_sample)

%---------------------------------------------------

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

Paras.T_0=0;              % initial temperature

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

[T,Xcor,Ycor] = GausFieldF(Paras);
X_int=T(:);

h = waitbar(0,'TestHDM is working very hard for you');
for i = 1:num_sample
    
    uParas.x0=X_sample(i,1);
    uParas.y0=X_sample(i,2);
   
    [y,Time]=Us_HeatDC_SolveFunc(Paras,uParas,X_int);
%     [Rec_X,Time]=Us_HeatDC_SolveFunc(Paras,IntParas)

    Y_Rec(i,:,1)=y(450,:); 
    Y_Rec(i,:,2)=y(490,:); 
    Y_Rec(i,:,3)=y(1960,:); 
    Y_Rec(i,:,4)=y(1990,:); 
    Y_Rec(i,:,5)=y(1250,:); 
% 
%     Sum_Y_Rec=sum(Y_Rec,2);
%     Sum_Y_Rec=reshape(Sum_Y_Rec,num_sample,[]);
    
%     Y_Rec(:,:,i)=y;
    


%     Y_Rec1(i,:,1)=y(50,:); 
%     Y_Rec2(i,:,2)=y(490,:); 
%     Y_Rec3(i,:,3)=y(1960,:); 
%     Y_Rec4(i,:,4)=y(1990,:); 
%     Y_Rec5(i,:,5)=y(1250,:); 
    
%     Y_Rec(:,:,i)=y; 
    Time_Rec(i,1)=Time;
        
%     y=y(:,Index_snapshot);
%     Y(i,:)=y(:)';
   
    waitbar(i/num_sample);
%     waitbar(i/Num_X,'TestHDM is working very hard for you');
    
%     h=waitbar(i/Num_X,'TestHDM is working very hard for you');
%     h=waitbar(i/Num_X,sprintf('TestHDM is working very hard for you. %d%% done...',i/Num_X));
    
end
close(h);

% for i=1:num_sample
%     Y_Rec2(i,:,1)=Y_Rec(84,:,i); 
%     Y_Rec2(i,:,2)=Y_Rec(96,:,i); 
%     Y_Rec2(i,:,3)=Y_Rec(324,:,i);  
%     Y_Rec2(i,:,4)=Y_Rec(336,:,i);  
%     Y_Rec2(i,:,5)=Y_Rec(200,:,i);   
% end





% save('HeatDC_MC_data21.mat','X_sample','Y_Rec','Time_Rec','Paras','uParas','mu','sigma','X_int');