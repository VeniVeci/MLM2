
%% GPE_MORx3_MC
%  GPE+MOR 3 method for multiple test cases
% Modifications:
% 30-3-2016, WeiX, first edition 

clear

% %------------------Basic Parameters--------------------------------------
% Paras.Lx=10;                % length of x axis is 1m
% Paras.Ly=10;                % length of y axis is 1m
% Paras.k=10;                 % Thermal conductivity is 1000W/m/K     % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)
% Paras.rho=1;                % flow density       % A sensitive parameter effects the stablility of solution(highter rho require more nodes for a grid)
% Paras.T_L=0;                % Temperature of Left    boundary 
% Paras.T_R=0;                % Temperature of Right   boundary 
% Paras.T_B=0;                % Temperature of Bottom  boundary 
% Paras.T_T=0;                % Temperature of Top     boundary 
% Paras.u_x=2;               % Velocity x direction 
% Paras.u_y=-20;              % Velocity y direction 
% 
% Paras.T_0=0;                % initial temperature
% 
% Paras.t_end=0.5;            % total simulation time 
% Paras.dt=0.01;              % time step
% 
% % -----------------------Solver Parameters---------------------------
% Paras.Num_node_x=20;        % number of volumes on x direction
% Paras.Num_node_y=20;        % number of volumes on y direction
% 
% % --------------------Complex Boundary condition parameters-------------
% Paras.b=100;
% Paras.a=10;
% %-----------------------------------------------------------------------
Num_Train=80;   
Num_Test=20;
Test_StartIndex=401;

Num_Snapshot=40;
n_Ubases=10;     %Number of POD bases

dim_new=10;      % For LPCA more than 10 require. The new dimension of reduced emulatior



%----------------Load dataset-------------------------------------------
load('ExpData.mat')
Num_Trian=141;
Index_test=152;%151 is difficulty

% X_star=[10,-10,5];
X_star=X(Index_test,:);
Y_orig=Y(Index_test,:);

Paras.u_x=X(Index_test,1);
Paras.u_y=X(Index_test,2);
Paras.a=X(Index_test,3);

% ---------------HDM----------------------------------------------------
[T_orig,Time_HDM]=Us_HeatDC_2D_FVM(Paras);
Time_HDM;


% ---------------MOR Bases by Record-----------------------------------
T_snaps_orig=reshape(Y_orig',Paras.Num_node_x*Paras.Num_node_y,[]);

approximate_degree=10;
[U_orig,S,~]=svd(T_snaps_orig);  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
U_orig=U_orig(:,1:approximate_degree);
eigenvalues=diag(S);

[T_U_orig,Time_U_orig]=Us_HeatDC_2D_FVM_Mor(Paras,U_orig);
Time_U_orig;




%  FiledPlot (Paras,T_orig,T_U_orig)

% -------------MOR Bases by GPR------------------------------------------
approximate_degree=10;
new_dim=10;     % For LPCA more than 10 require
% -------------------------------LPCA-GPR--------------------------------
% [Y_star,~,time]=GPR_SVD(X(1:Num_Trian,:),Y(1:Num_Trian,:),X_star,new_dim);

% % -------------------------------KPCA-GPR--------------------------------
% options = struct('ker','rbf','arg',50000,'new_dim',new_dim); 
% [Y_star,~,time]=GPR_KPCA(X(1:Num_Trian,:),Y(1:Num_Trian,:),X_star,options);

% -------------------------------LPCA-ANN--------------------------------
[Y_star,time]=ANN_SVD(X(1:Num_Trian,:)',Y(1:Num_Trian,:)',X_star',new_dim,10); %works fine
Y_star=Y_star';

% -------------------------------KPCA-ANN--------------------------------
% options = struct('ker','rbf','arg',50000,'new_dim',new_dim); 
% [Y_star,time]=ANN_KPCA(X(1:Num_Trian,:)',Y(1:Num_Trian,:)',X_star',options,10);
% Y_star=Y_star';

%-------------------------------Rearrange---------------------------------
T_snaps_star=reshape(Y_star',Paras.Num_node_x*Paras.Num_node_y,[]);


% %-----------------Show prediction by GPR and original---------------------
% T_snaps_orig=reshape(Y_orig',Paras.Num_node_x*Paras.Num_node_y,[]);
% FiledPlot (Paras,T_snaps_orig,T_snaps_star);


%--------------------------------MOR ------------------------------------
[U_star,S,~]=svd(T_snaps_star);  % U*S*V'=Rec_X
U_star=U_star(:,1:approximate_degree);
eigenvalues=diag(S);


[T_U_star,Time_U_star]=Us_HeatDC_2D_FVM_Mor(Paras,U_star);
Time_U_star;

% FiledPlot (Paras,T_U_star)
% FiledPlot (Paras,T_U_orig-T_orig)
% FiledPlot (Paras,T_U_orig,T_U_star)


% Make sure the index end means the same frame!!!
figure()
hold on
L1=plot(T_orig(:,end),'k-');
L2=plot(T_snaps_star(:,end),'b--');
L3=plot(T_U_star(:,end),'r-.');
legend([L1,L2,L3],'Original T filed','GPR Predicted T filed','GPR-MOR T filed')
hold off

% FiledPlot (Paras,T_orig,T_U_orig,T_U_star)


%% Additional try % not working

T_Ref=T_snaps_star;

Plotline=[];
for i =1:10
   
    [U_star,S,~]=svd(T_Ref);  % U*S*V'=Rec_X
    U_star=U_star(:,1:approximate_degree);
    eigenvalues=diag(S);

    [T_U_star,Time_U_star]=Us_HeatDC_2D_FVM_Mor(Paras,U_star);
    Time_U_star;
    
       
%     L1=plot(T_Ref(:,end),'k-');
    
    T_Ref=T_U_star;
    Plotline=[Plotline,T_U_star(:,end)];
end

figure()
hold on
L1=plot(Plotline);
L1=plot(T_orig(:,end),'k--');





