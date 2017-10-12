%% POD_Burger1D_Part2_POD_DEIM_M
% S for Single case
%
% Modifications:
% 16-Jun-2016, WeiX, first edition 
clear

%%---------------Setting parameter---------------------------------------
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best

Index_Train=[1:80];
Index_Test=[301:400];

Index_SS=[1:2:201]; % Index of snapshot


n_Ubases=10;     %Number of POD bases
n_DEIM=15;


%% DR method and parameters
dim_new=10;           
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
load('ExpDataT2_2.mat') 
% load('ExpDataT2_1.mat') 
% load('ExpDataT2_2F1.mat') 


X_Trian=X(Index_Train,:);
Y_Trian=Y_Rec(:,Index_SS,Index_Train);

% X_star=[10,-10,5];
X_star=X(Index_Test,:);
Y_orig=Y_Rec(:,:,Index_Test);

Paras.Re=X(Index_Test,1);
Paras.u0a=X(Index_Test,2);
Paras.u0b=X(Index_Test,3);

%% ---------------G-POD Full -----------------------------------
display('Starting Global POD ...')

Y_Global_snaps=reshape(Y_Trian,Paras.n+1,[]);

[U_GlobalSnapS,S,~]=svd(Y_Global_snaps(2:end-1,:));  % U*S*V'=Rec_X

U_GlobalSnapS=U_GlobalSnapS(:,1:n_Ubases);
eigenvalues_GlobalSnapS=diag(S);

for i=1:length(Index_Test)
    
    Paras.Re=X(Index_Test(i),1);
    Paras.u0a=X(Index_Test(i),2);
    Paras.u0b=X(Index_Test(i),3);
    
    [Y_U_GlobalSnapS,~,Time_ROM_U_GlobalSnapS]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_GlobalSnapS);    
    Y_GU(:,:,i)=Y_U_GlobalSnapS;   
    Time_GU(i)=Time_ROM_U_GlobalSnapS;  
    
end

%% --------------- G-POD DEIM Global-----------------------------------
display('Starting Global POD DEIM...')
Y_F=Y_Global_snaps.^2;
[UF,~,~]=svd(Y_F(2:end-1,:));  % F(y)=y.^2 in Burgers' equation
UF=UF(:,1:n_DEIM);
[phi,UG,P] = DEIM(UF);

for i=1:length(Index_Test)
    
    Paras.Re=X(Index_Test(i),1);
    Paras.u0a=X(Index_Test(i),2);
    Paras.u0b=X(Index_Test(i),3);

    [Y_GPOD_DEIM,~,Time_GPOD_DEIM]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U_GlobalSnapS,UF,P);
    Y_GU_DEIM(:,:,i)=Y_GPOD_DEIM;   
    Time_GU_DEIM(i)=Time_GPOD_DEIM;  
    
end


%% --------------- G-POD GP DEIM + GP-POD GPE DEIM  -----------------------------------
display('Starting Global POD GP DEIM...')

% Rearrange
Y_Train_v=reshape(Y_Trian,[],length(Index_Train))';
% Y_Trian2 =reshape(Y_Train_v',Paras.n+1,[],length(Index_Train));  % Recover formular
[Y_GPE,Time_Emu]=Func_DrGPE(X_Trian,Y_Train_v,X_star,options,preoptions);
% --Rearrange--   
Y_GPESnapS=reshape(Y_GPE',Paras.n+1,[],length(Index_Test));

% --MOR FEM--
for i=1:length(Index_Test)
    
    Y_F=Y_GPESnapS(:,:,i).^2;
    [UF,~,~]=svd(Y_F(2:end-1,:),'econ');  % F(y)=y.^2 in Burgers' equation
    UF=UF(:,1:n_DEIM);
    [phi,Uu,P] = DEIM(UF);

    Paras.Re= X_star(i,1);
    Paras.u0a=X_star(i,2);
    Paras.u0b=X_star(i,3);
    
    [Y_GPOD_GPEDEIM,~,Time_GPOD_GPEDEIM]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U_GlobalSnapS,UF,P);    
    Y_GU_GPEDEIM(:,:,i)=Y_GPOD_GPEDEIM;   
    Time_GU_GPEDEIM(i)=Time_GPOD_GPEDEIM;  
    
    
    %---------GP-POD GP DEIM  ------------------
    [U_GPESnapS,S,~]=svd(Y_GPESnapS(2:end-1,:,i));  % U*S*V'=Rec_X
    U_GPESnapS=U_GPESnapS(:,1:n_Ubases);
    
    [Y_GPEPOD_GPEDEIM,~,Time_GPEPOD_GPEDEIM]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U_GPESnapS,UF,P);    
    Y_GPEU_GPEDEIM(:,:,i)=Y_GPEPOD_GPEDEIM;   
    Time_GPEU_GPEDEIM(i)=Time_GPEPOD_GPEDEIM;  
    
end

%% -----------------Error Analysis---------------------------------------
% SE_U_OrigSnapS=(Y_U_OrigSnapS-Y_orig).^2;
SE_GU=(Y_GU-Y_orig).^2;
SE_GU_DEIM=(Y_GU_DEIM-Y_orig).^2;
SE_GU_GPEDEIM=(Y_GU_GPEDEIM-Y_orig).^2;
SE_GPEU_GPEDEIM=(Y_GPEU_GPEDEIM-Y_orig).^2;

SSE_GU=sum(reshape(SE_GU,[],length(Index_Test)));
SSE_GU_DEIM=sum(reshape(SE_GU_DEIM,[],length(Index_Test)));
SSE_GU_GPEDEIM=sum(reshape(SE_GU_GPEDEIM,[],length(Index_Test)));
SSE_GPEU_GPEDEIM=sum(reshape(SE_GPEU_GPEDEIM,[],length(Index_Test)));


RE_GU=SSE_GU./sum(reshape(Y_orig.^2,[],length(Index_Test)));
RE_GU_DEIM=SSE_GU_DEIM./sum(reshape(Y_orig.^2,[],length(Index_Test)));
RE_GU_GPEDEIM=SSE_GU_GPEDEIM./sum(reshape(Y_orig.^2,[],length(Index_Test)));
RE_GPEU_GPEDEIM=SSE_GPEU_GPEDEIM./sum(reshape(Y_orig.^2,[],length(Index_Test)));


MSSE_GU=mean(SSE_GU);
MSSE_GU_DEIM=mean(SSE_GU_DEIM);
MSSE_GU_GPEDEIM=mean(SSE_GU_GPEDEIM);
MSSE_GPEU_GPEDEIM=mean(SSE_GPEU_GPEDEIM);


figure
boxplot([RE_GU',RE_GU_DEIM',RE_GU_GPEDEIM',RE_GPEU_GPEDEIM'],'labels',{'Global POD','Global POD DEIM','GLobal POD GPE-DEIM','GPE POD GPE-DEIM'})
title(sprintf('Rate Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DEIM}=%i, Num_{DimNew}=%i', DrMethod,length(Index_Train),length(Index_Test),length(Index_SS),n_Ubases,n_DEIM,dim_new))
hold off
set(gca,'yscale','log');









