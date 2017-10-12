%% POD_Burger1D_Part2_POD_M
% M for Multiple case
%
% Modifications:
% 16-Jun-2016, WeiX, first edition 
clear

%%---------------Setting parameter---------------------------------------
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best

Index_Train=[1:120];
Index_Test=[301:500];

Index_SS=[1:4:201]; % Index of snapshot


n_Ubases=10;     %Number of POD bases
% n_DEIM=10;


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

X_Trian=X(Index_Train,:);
Y_Trian=Y_Rec(:,Index_SS,Index_Train);

% X_star=[10,-10,5];
X_star=X(Index_Test,:);
Y_orig=Y_Rec(:,:,Index_Test);

Paras.Re=X(Index_Test,1);
Paras.u0a=X(Index_Test,2);
Paras.u0b=X(Index_Test,3);

%% ---------------MOR Bases by perfect snapshot -----------------------------------
% % Y_orig_snaps=reshape(Y_orig',Paras.n+1,[]);
% Y_orig_snaps=Y_orig(:,Index_SS);
% 
% % n_Ubases=5;
% % Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point
% [U_OrigSnapS,S,~]=svd(Y_orig_snaps(2:end-1,:));  % U*S*V'=Rec_X
% % [U_orig,S,~]=svd(T_orig);
% U_OrigSnapS=U_OrigSnapS(:,1:n_Ubases);
% eigenvalues_OrigSnapS=diag(S);
% 
% [Y_U_OrigSnapS,~,Time_ROM_U_OrigSnapS]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_OrigSnapS);

%% ---------------MOR Global Bases -----------------------------------
display('Starting Global POD...')
% Y_Global_snaps=reshape(Y(1:Num_Trian,:)',Paras.n+1,[]);
Y_Global_snaps=reshape(Y_Trian,Paras.n+1,[]);

% n_Ubases=5;
% Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point
[U_GlobalSnapS,S,~]=svd(Y_Global_snaps(2:end-1,:));  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);
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

%% -------------MOR Bases by GPE Predicted snapshot-------------------------
display('Starting GPES POD...')
% n_Ubases=5;   %number of POD basics
Y_Train_v=reshape(Y_Trian,[],length(Index_Train))';
% Y_Trian2 =reshape(Y_Train_v',Paras.n+1,[],length(Index_Train));

[Y_GPE,Time_Emu]=Func_DrGPE(X_Trian,Y_Train_v,X_star,options,preoptions);
% --Rearrange--   
Y_GPESnapS=reshape(Y_GPE',Paras.n+1,[],length(Index_Test));

% --MOR FEM--
for i=1:length(Index_Test)
    
    [U_GPESnapS,S,~]=svd(Y_GPESnapS(2:end-1,:,i));  % U*S*V'=Rec_X
    U_GPESnapS=U_GPESnapS(:,1:n_Ubases);
%     eigenvalues_Ustar=diag(S);

    Paras.Re= X_star(i,1);
    Paras.u0a=X_star(i,2);
    Paras.u0b=X_star(i,3);
    
    [Y_U_GPESnapS,~,Time_ROM_U_GPESnapS]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_GPESnapS);
    Y_GPESU(:,:,i)=Y_U_GPESnapS;   
    Time_GPESU(i)=Time_ROM_U_GPESnapS;  
end


%% -------------MOR Bases by GPE Prediction-------------------------
display('Starting GPEU POD...')
% -------------Figure local basis-------------------------
% Figure local basis
for i=1:length(Index_Train)

    Yi_SnapS=Y_Trian(:,:,i);
%     Yi=Y(i,:);
%     Yi_SnapS=reshape(Yi',Paras.n+1,[]);
    [U_YiSnapS,~,~]=svd(Yi_SnapS(2:end-1,:),'econ');  % U*S*V'=Rec_X
    U_YiSnapS=U_YiSnapS(:,1:n_Ubases);

    U_Y_Train(i,:)=U_YiSnapS(:)';

end
 % ------------ Main -------------------------
 
[U_GPE,Time_Emu]=Func_DrGPE(X_Trian,U_Y_Train,X_star,options,preoptions);     
% --Rearrange--
U_GPE=reshape(U_GPE',Paras.n-1,[],length(Index_Test));

for i=1:length(Index_Test)
    
    Paras.Re= X_star(i,1);
    Paras.u0a=X_star(i,2);
    Paras.u0b=X_star(i,3);
    
    [Y_U_GPE,~,Time_ROM_U_GPE]=Burger1D_FEM_DBC_MOR_SolverF(Paras,U_GPE(:,:,i));
    Y_GPEU(:,:,i)=Y_U_GPE;   
    Time_GPEU(i)=Time_ROM_U_GPE;  
end



%% -----------------Error Analysis---------------------------------------
% SE_U_OrigSnapS=(Y_U_OrigSnapS-Y_orig).^2;
SE_GPESU=(Y_GPESU-Y_orig).^2;
SE_GU=(Y_GU-Y_orig).^2;
SE_GPEU=(Y_GPEU-Y_orig).^2;

SSE_GPESU=sum(reshape(SE_GPESU,[],length(Index_Test)));
SSE_GU=sum(reshape(SE_GU,[],length(Index_Test)));
SSE_GPEU=sum(reshape(SE_GPEU,[],length(Index_Test)));

RE_GPESU=SSE_GPESU./sum(reshape(Y_orig.^2,[],length(Index_Test)));
RE_GU=SSE_GU./sum(reshape(Y_orig.^2,[],length(Index_Test)));
RE_GPEU=SSE_GPEU./sum(reshape(Y_orig.^2,[],length(Index_Test)));



MSSE_GPESU=mean(SSE_GPESU);
MSSE_GU=mean(SSE_GU);
MSSE_GPEU=mean(SSE_GPEU);

% MSSE_GPESU=sum(SE_GPESU(:))/length(Index_Test);
% MSSE_GU=sum(SE_GU(:))/length(Index_Test);
% MSSE_GPEU=sum(SE_GPEU(:))/length(Index_Test);

% figure
% boxplot([SSE_GPESU',SSE_GU',SSE_GPEU']);
% figure
% boxplot([RE_GPESU',RE_GU',RE_GPEU']);

figure
boxplot([RE_GPESU',RE_GU',RE_GPEU'],'labels',{'GPE POD','Global POD','GPEU POD'})
title(sprintf('Rate Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DimNew}=%i', DrMethod,length(Index_Train),length(Index_Test),length(Index_SS),n_Ubases,dim_new))
hold off
set(gca,'yscale','log');


% % ----Animation plot------------------------------------------------------
% x = 0:1/Paras.n:1;  % coordinate sequence
% Y_Maxi=max([Y_orig(:);Y_U_OrigSnapS(:);Y_U_GPESnapS(:)]);
% Y_Mini=min([Y_orig(:);Y_U_OrigSnapS(:);Y_U_GPESnapS(:)]);
% for i = 1:Paras.t_n+1
%     figure(5)   
% %     title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n))))
% %     subplot(2,1,1)
%     L1=plot(x,Y_orig(:,i),'k-');
%     hold on
%     L2=plot(x,Y_U_OrigSnapS(:,i),'b--');
%     L3=plot(x,Y_U_GPESnapS(:,i),'r-.');
%     L4=plot(x,Y_U_GlobalSnapS(:,i),'g--');
%     L5=plot(x,Y_U_GPE(:,i),'y--');
%     
%     axis([0,1,Y_Mini,Y_Maxi]);
%     legend([L1,L2,L3,L4,L5],'Full FEM Solution','Perfect MOR FEM Solution','GPE Snapshot bases MOR FEM Solution','Global bases MOR FEM Solution','GPE bases MOR FEM Solution')
%     title(sprintf('Animation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n)))) 
%     hold off
%     
% % %     figure(5)
% %     subplot(2,1,2)
% %     L1=plot(x(1:i),sum(Y_Full(:,1:i)-Y_Full(:,1:i)),'k-');
% %     hold on
% %     L2=plot(x(1:i),sum(Y_U_OrigSnapS(:,1:i)-Y_Full(:,1:i)),'b--');
% %     L3=plot(x(1:i),sum(Y_U_GPESnapS(:,1:i)-Y_Full(:,1:i)),'r-.');
% %     L4=plot(x(1:i),sum(Y_U_GlobalSnapS(:,1:i)-Y_Full(:,1:i)),'g--');
% %     L5=plot(x(1:i),sum(Y_U_GPE(:,1:i)-Y_Full(:,1:i)),'y--');
% %     
% %     axis([0,1,Y_Mini,Y_Maxi]);
% %     legend([L1,L2,L3,L4,L5],'Full FEM Solution','Perfect MOR FEM Solution','GPE Snapshot bases MOR FEM Solution','Global bases MOR FEM Solution','GPE bases MOR FEM Solution')
% %     title(sprintf('L^2 Error Acumulation t= %0.3f',((i-1)*(Paras.t_end/Paras.t_n)))) 
% %     hold off
% 
%     F(i) = getframe;    
% end    



