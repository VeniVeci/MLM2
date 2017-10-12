%% MC_HeatDC_Run_Part2
%Monte Carlo process using surrogate.
%
% Modifications:
% WeiX, 11-4-2016, Create.


clear
%% ----HDM MC-----

% load('HeatDC_MC_data13.mat')
% 
% [num_sample,~,~]=size(Y_Rec);
% 
% Sum_Y_Rec=sum(Y_Rec,2);
% Sum_Y_Rec=reshape(Sum_Y_Rec,num_sample,[]);
% 
% 
% SS_Y_Rec=sum(Sum_Y_Rec,2);
% med=median(SS_Y_Rec);
% 
% index=find(SS_Y_Rec<5*med & SS_Y_Rec>=0);
% 
% Sum_Y_Rec=Sum_Y_Rec(index,:);

%%---------------Setting parameter---------------------------------------
Num_Train=80;   
% Num_Test=100;
Test_StartIndex=301;

Num_Snapshot=10;

n_Ubases=50;     %Number of POD bases
dim_new=10;      % For LPCA more than 10 require. The new dimension of reduced emulatior

num_sample=3000;
mu = [0.5,0.5];
sigma = [0.01,0;0,0.01];

%% DR method and parameters     
DrMethod='kPCA';
% DrMethod='DiffusionMaps';
% DrMethod='Isomaps';
% DrMethod='PCA';

switch DrMethod
    
    case 'kPCA'
        options.ker='gaussian';   
        options.new_dim=dim_new;
        options.FullRec=0;       
        options.arg=1000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
        options.kAuto=1;
        
    case 'DiffusionMaps'
        options.metric ='euclidean';
        options.kernel ='gaussian'; 
        options.dim_new = dim_new;              
        options.t = 1;                     
        options.FullRec = 0;      
        % Doptions.kpara = 10000;             
        options.kAuto=1;
        options.Ztype=0;    %Output type. With/without 1st component
        
    case 'Isomaps'
        options.dim_new=dim_new;                % New dimension
        options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
        options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
        options.metric='euclidean';             % Method of measurement. Metric

        %Isomap PreImage options
        preoptions.ReCoverNeighborType='k';     % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
        preoptions.ReCoverNeighborPara=10;      % Parameter of neighbor of new point
        
    case 'PCA'
        options=[];
        
    otherwise 
        error('No such DR method')
end


options.DrMethod=DrMethod;
% PreImage options----------------------
preoptions.type='Exp';  %'LSE', 'Dw' OR 'Exp'
preoptions.neighbor=10;
% preoptions.type='LpcaI';
% preoptions.dim_new=10; % Use to stable the result but sacrefy accuracy
% preoptions.InMethod='ANN';


%% ----------------Load dataset-------------------------------------------
% load('HeatDC_DBC_ExpData21.mat')
load('HeatDC_DBC_ExpData13.mat')
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best

% X_star=[10,-10,5];

Paras.t_n=Paras.t_end/Paras.dt;
Index_snapshot=1:(Paras.t_n/Num_Snapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

[~,~,num_Y]=size(Y_Rec);
Y=Y_Rec(:,Index_snapshot,:);
Y=reshape(Y,[],num_Y); %500 for this dataset
Y=Y';
% Time_HDM=Time;

% X_star=X(Test_StartIndex:Test_StartIndex+Num_Test-1,:);
% Y_starorig=Y(Test_StartIndex:Test_StartIndex+Num_Test-1,:);

X=X(1:Num_Train,:);
Y=Y(1:Num_Train,:);

% Paras.Re=X(Index_test,1);
% Paras.u0a=X(Index_test,2);
% Paras.u0b=X(Index_test,3);

%% ---------------MOR Golboal Bases -----------------------------------
Y_Global_snaps=reshape(Y',Paras.Num_node_y*Paras.Num_node_x,[]);
% n_Ubases=5;
% Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point

% [U_GlobalSnapS,~,~]=svd(Y_Global_snaps);  % U*S*V'=Rec_X
[U_GlobalSnapS,~,~]=svds(Y_Global_snaps,n_Ubases); % Accelerate SVD decomposition 

% [U_orig,S,~]=svd(T_orig);
U_GlobalSnapS=U_GlobalSnapS(:,1:n_Ubases);
% eigenvalues_GlobalSnapS=diag(S);

% num_sample=3000;
% mu = [0.5,0.5];
% sigma = [0.01,0;0,0.01];
% rng default  % For reproducibility
rng(2)
X_sample = mvnrnd(mu,sigma,num_sample);

tic;
h = waitbar(0,'Test MOR Bases by Golboal Bases');
for i =1:num_sample
    
    
    uParas.x0=X_sample(i,1);
    uParas.y0=X_sample(i,2);

    [y,Time]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U_GlobalSnapS);
%     [Rec_X,Time]=Us_HeatDC_SolveFunc(Paras,IntParas)

%     Y_UGlobal_Rec(i,:,1)=y(84,:); 
%     Y_UGlobal_Rec(i,:,2)=y(96,:); 
%     Y_UGlobal_Rec(i,:,3)=y(324,:); 
%     Y_UGlobal_Rec(i,:,4)=y(336,:); 
%     Y_UGlobal_Rec(i,:,5)=y(200,:); 
    
    Y_UGlobal_Rec(i,:,1)=y(450,:); 
    Y_UGlobal_Rec(i,:,2)=y(490,:); 
    Y_UGlobal_Rec(i,:,3)=y(1960,:);
    Y_UGlobal_Rec(i,:,4)=y(1990,:);
    Y_UGlobal_Rec(i,:,5)=y(1250,:);
        
    Time_UGlobal_Rec(i,1)=Time;
    
    waitbar(i/num_sample);
end
close(h);
totalTime_Global=toc;
% save('HeatDC_MC_data13.mat','Y_UGlobal_Rec','Time_UGlobal_Rec','-append')

% Data process
% Sum_Y_UGlobal_Rec=sum(Y_UGlobal_Rec,2);
% Sum_Y_UGlobal_Rec=reshape(Sum_Y_UGlobal_Rec,num_sample,[]);
% 
% 
% SS_Y_UGlobal_Rec=sum(Sum_Y_UGlobal_Rec,2);
% med=median(SS_Y_UGlobal_Rec);
% index=find(SS_Y_UGlobal_Rec<5*med & SS_Y_UGlobal_Rec>=0);
% 
% Sum_Y_UGlobal_Rec=Sum_Y_UGlobal_Rec(index,:);


%% -------------MOR Bases by GPE Predicted snapshot-------------------------
% n_Ubases=5;   %number of POD basics
 

% num_sample=3000;
% mu = [0.5,0.5];
% sigma = [0.01,0;0,0.01];
% rng default  % For reproducibility
rng(1)
X_sample = mvnrnd(mu,sigma,num_sample);


[Y_GPE,Time_Emu]=Func_DrGPE(X(1:Num_Train,:),Y(1:Num_Train,:),X_sample,options,preoptions);    

% toc0=toc;
tic;
h = waitbar(0,'Test MOR Bases by GPE snapshot');
for i =1:num_sample
    
    uParas.x0=X_sample(i,1);
    uParas.y0=X_sample(i,2);

    Y_GPESnapS =reshape(Y_GPE(i,:)',Paras.Num_node_y*Paras.Num_node_x,[]);
  
%     [U_GPESnapS,~,~]=svd(Y_GPESnapS);  % U*S*V'=Rec_X
    [U_GPESnapS,~,~]=svds(Y_GPESnapS,n_Ubases); % Accelerate SVD decomposition 
    
    U_GPESnapS=U_GPESnapS(:,1:n_Ubases);
%     eigenvalues_Ustar=diag(S);
    [y,Time]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U_GPESnapS);
    
%     Y_GPESnapS_Rec(i,:,1)=y(84,:); 
%     Y_GPESnapS_Rec(i,:,2)=y(96,:); 
%     Y_GPESnapS_Rec(i,:,3)=y(324,:); 
%     Y_GPESnapS_Rec(i,:,4)=y(333,:); 
%     Y_GPESnapS_Rec(i,:,5)=y(200,:); 

    Y_GPESnapS_Rec(i,:,1)=y(450,:); 
    Y_GPESnapS_Rec(i,:,2)=y(490,:); 
    Y_GPESnapS_Rec(i,:,3)=y(1960,:);
    Y_GPESnapS_Rec(i,:,4)=y(1990,:);
 	Y_GPESnapS_Rec(i,:,5)=y(1250,:);
    
    
    Time_GPESnapS_Rec(i,1)=Time;
    
    waitbar(i/num_sample);
end
close(h);
% Time_ROM_U_GPESnapS=toc-toc0;
totalTime_GPEPOD=toc;

% save('HeatDC_MC_data13.mat','Y_GPESnapS_Rec','Time_GPESnapS_Rec','-append')

% Data process
Sum_Y_GPESnapS_Rec=sum(Y_GPESnapS_Rec,2);
Sum_Y_GPESnapS_Rec=reshape(Sum_Y_GPESnapS_Rec,num_sample,[]);


SS_Y_GPESnapS_Rec=sum(Sum_Y_GPESnapS_Rec,2);
med=median(SS_Y_GPESnapS_Rec);
index=find(SS_Y_GPESnapS_Rec<5*med & SS_Y_GPESnapS_Rec>=0);

Sum_Y_GPESnapS_Rec=Sum_Y_GPESnapS_Rec(index,:);


Para_Dr.DrMethod=DrMethod;
Para_Dr.Num_Train=Num_Train;
Para_Dr.Num_Snapshot=Num_Snapshot;
Para_Dr.n_Ubases=n_Ubases;
Para_Dr.dim_new=dim_new;
% save('HeatDC_MC_data13.mat','Para_Dr','-append')


% save('HeatDC_MC_data21_1.mat','Y_UGlobal_Rec','Time_UGlobal_Rec','Y_GPESnapS_Rec','Time_GPESnapS_Rec','Para_Dr')
%  save('HeatDC_MC_data13_SS100U30.mat','Y_UGlobal_Rec','Time_UGlobal_Rec','Y_GPESnapS_Rec','Time_GPESnapS_Rec','Para_Dr','totalTime_Global','totalTime_GPEPOD')




