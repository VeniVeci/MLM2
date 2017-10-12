%% POD_Burger1D_Part2_POD_DEIM_S
% S for Single case
%
% Modifications:
% 15-Jun-2016, WeiX, first edition 
clear
%%---------------Setting parameter---------------------------------------
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best

Index_Train=[1:160];
Index_Test=[408];
Index_Test=[438];

Index_SS=[1:2:201]; % Index of snapshot


n_Ubases=10;     %Number of POD bases
n_DEIM=10;


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
load('ExpDataT2_1.mat') 

X_Trian=X(Index_Train,:);
Y_Trian=Y_Rec(:,Index_SS,Index_Train);

% X_star=[10,-10,5];
X_star=X(Index_Test,:);
Y_orig=Y_Rec(:,:,Index_Test);

Paras.Re=X(Index_Test,1);
Paras.u0a=X(Index_Test,2);
Paras.u0b=X(Index_Test,3);



%% ---------------G-POD Full -----------------------------------
% Y_Global_snaps=reshape(Y(1:Num_Trian,:)',Paras.n+1,[]);
Y_Global_snaps=reshape(Y_Trian,Paras.n+1,[]);

% n_Ubases=5;
% Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point
[U_GlobalSnapS,S,~]=svd(Y_Global_snaps(2:end-1,:));  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
U_GlobalSnapS=U_GlobalSnapS(:,1:n_Ubases);
eigenvalues_GlobalSnapS=diag(S);

[Y_GPOD,~,Time_GPOD]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_GlobalSnapS);

%% --------------- G-POD DEIM Global-----------------------------------
Y_F=Y_Global_snaps.^2;
[UF,~,~]=svd(Y_F(2:end-1,:));  % F(y)=y.^2 in Burgers' equation
UF=UF(:,1:n_DEIM);
[phi,Uu,P] = DEIM(UF);

[Y_GPOD_DEIM,~,Time_GPOD_DEIM]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U_GlobalSnapS,UF,P);


%% --------------- G-POD DEIM GP-----------------------------------
Y_Train_v=reshape(Y_Trian,[],length(Index_Train))';

Y_F_Train_v=Y_Train_v.^2;
[Y_F_GPE,Time_Emu]=Func_DrGPE(X_Trian,Y_F_Train_v,X_star,options,preoptions);     

% --Rearrange--
Y_F_GPE=reshape(Y_F_GPE',Paras.n+1,[]);

[UF_GP,~,~]=svd(Y_F_GPE(2:end-1,:));  % F(y)=y.^2 in Burgers' equation
UF_GP=UF_GP(:,1:n_DEIM);
[phi,Uu,P] = DEIM(UF_GP);

[Y_GPOD_GPDEIM,~,Time_GPOD_GPDEIM]=Burger1D_FEM_DBC_MOR_DEIM_SolverF(Paras,U_GlobalSnapS,UF_GP,P);





%% -----------------Error Analysis---------------------------------------
SE_GPOD=(Y_GPOD-Y_orig).^2;
SE_GPOD_DEIM=(Y_GPOD_DEIM-Y_orig).^2;
SE_GPOD_GPDEIM=(Y_GPOD_GPDEIM-Y_orig).^2;


% ----Animation plot------------------------------------------------------
x = 0:1/Paras.n:1;  % coordinate sequence
% Y_Maxi=max([Y_orig(:);Y_U_OrigSnapS(:);Y_U_GPESnapS(:)]);
% Y_Mini=min([Y_orig(:);Y_U_OrigSnapS(:);Y_U_GPESnapS(:)]);
Y_Maxi=max(Y_orig(:));
Y_Mini=min(Y_orig(:));
for i = 1:Paras.t_n+1
    figure(5)   
%     title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))
%     subplot(2,1,1)
    L1=plot(x,Y_orig(:,i),'k-');
    hold on
    L2=plot(x,Y_GPOD(:,i),'b--');
    L3=plot(x,Y_GPOD_DEIM(:,i),'r-.');
    L4=plot(x,Y_GPOD_GPDEIM(:,i),'g--');
   
    axis([0,1,Y_Mini,Y_Maxi]);
    legend([L1,L2,L3,L4],'Full FEM Solution','POD Full','POD DEIM','POD GP+DEIM')
    title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n)))) 
    hold off
    
% %     figure(5)
%     subplot(2,1,2)
%     L1=plot(x(1:i),sum(Y_Full(:,1:i)-Y_Full(:,1:i)),'k-');
%     hold on
%     L2=plot(x(1:i),sum(Y_U_OrigSnapS(:,1:i)-Y_Full(:,1:i)),'b--');
%     L3=plot(x(1:i),sum(Y_U_GPESnapS(:,1:i)-Y_Full(:,1:i)),'r-.');
%     L4=plot(x(1:i),sum(Y_U_GlobalSnapS(:,1:i)-Y_Full(:,1:i)),'g--');
%     L5=plot(x(1:i),sum(Y_U_GPE(:,1:i)-Y_Full(:,1:i)),'y--');
%     
%     axis([0,1,Y_Mini,Y_Maxi]);
%     legend([L1,L2,L3,L4,L5],'Full FEM Solution','Perfect MOR FEM Solution','GPE Snapshot bases MOR FEM Solution','Global bases MOR FEM Solution','GPE bases MOR FEM Solution')
%     title(sprintf('L^2 Error Acumulation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n)))) 
%     hold off

    F(i) = getframe;    
end    



