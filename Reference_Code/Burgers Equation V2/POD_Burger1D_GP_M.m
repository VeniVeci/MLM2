%% POD_Burger1D_GP_M
% S for Single case
%
% Modifications:
% 16-Jun-2016, WeiX, first edition 
clear

Index_Train=[1:80];
Index_Test=[201:500];

Index_SS=[1:2:201]; % Index of snapshot

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
Y_orig1=Y_Rec(:,:,Index_Test);
Y_orig=Y_orig1(:,Index_SS,:);


%% ---------------- GP -------------------------------------------

% Rearrange
Y_Train_v=reshape(Y_Trian,[],length(Index_Train))';
% Y_Trian2 =reshape(Y_Train_v',Paras.n+1,[],length(Index_Train));  % Recover formular
[Y_GPE,Time_Emu]=Func_DrGPE(X_Trian,Y_Train_v,X_star,options,preoptions);
% --Rearrange--   
Y_GPESnapS=reshape(Y_GPE',Paras.n+1,[],length(Index_Test));



%% -----------------Error Analysis---------------------------------------
% SE_U_OrigSnapS=(Y_U_OrigSnapS-Y_orig).^2;
SE=(Y_GPESnapS-Y_orig).^2;
SSE=sum(reshape(SE,[],length(Index_Test)));
RE=SSE./sum(reshape(Y_orig.^2,[],length(Index_Test)));

MSSE=mean(SSE);
figure
boxplot(RE');
title(sprintf('Rate Error Boxplot. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i,  Num_{DimNew}=%i', DrMethod,length(Index_Train),length(Index_Test),length(Index_SS),dim_new))
hold off
set(gca,'yscale','log');


%% -----------------Filter---------------------------------------

[Sort_RE,index]=sort(RE');

% index_delete=index(end-20:end);
index_select=index(1:end-20);
index_all=sort(index_select);


% index_all=[1:200,index_select']';
% index_all=sort(index_all);
% Y_Rec=[Y_Rec(:,:,1:200);Y_orig1(:,:,index_all)];

X=[X(1:200,:);X_star(index_all,:)];
Y_Rec=Y_Rec(:,:,1:200);
Y_Rec(:,:,201:480)=Y_orig1(:,:,index_all);

%     save('ExpDataT2_2.mat','X','Y_Rec','Time','T1','Paras','rang');




