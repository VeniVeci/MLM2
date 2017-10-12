%% MOR_IncBase
%  Test ROM wiht increasing basis number
%
% Modifications:
% 12-4-2016, WeiX, first edition 
% 05-01-2017, WeiX, Modify relative error showcase

clear

%%---------------Setting parameter---------------------------------------
Num_Train=80;   
Num_Test=300;
Test_StartIndex=201;

Num_Snapshot=100;

% n_Ubases=50;     
n_Ubases=[5:5:50]; %Number of POD bases

dim_new=10;      % For LPCA more than 10 require. The new dimension of reduced emulatior

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
load('HeatDC_DBC_ExpData13.mat')
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best

% X_star=[10,-10,5];
[~,~,num_Y]=size(Y_Rec);

Paras.t_n=Paras.t_end/Paras.dt;
Index_snapshot=1:(Paras.t_n/Num_Snapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

Y=Y_Rec(:,Index_snapshot,:);
Y=reshape(Y,[],num_Y); %500 for this dataset
Y=Y';
% Time_HDM=Time;

X_star=X(Test_StartIndex:Test_StartIndex+Num_Test-1,:);
Y_starorig=Y(Test_StartIndex:Test_StartIndex+Num_Test-1,:);

X=X(1:Num_Train,:);
Y=Y(1:Num_Train,:);

Y_Rec_star=Y_Rec(:,:,Test_StartIndex:end);

%% ---------------MOR Golboal Bases -----------------------------------
Y_Global_snaps=reshape(Y',Paras.Num_node_y*Paras.Num_node_x,[]);
% n_Ubases=5;
% Y_orig_snaps=Y_orig_snaps(2:end-1,:); %take out boundary point
[U_GlobalSnapS,~,~]=svd(Y_Global_snaps);  % U*S*V'=Rec_X
% [U_orig,S,~]=svd(T_orig);

toc0=toc;
h = waitbar(0,'Test MOR Bases by Golboal Bases');
for j=n_Ubases
    
   U=U_GlobalSnapS(:,1:j);
    
   for i =1:Num_Test    
        %     uParas.ampx=1;
        %     uParas.ampy=1;
        uParas.x0=X_star(i,1);
        uParas.y0=X_star(i,2);

        [Y_U,~]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U);

        Y_orig=Y_Rec_star(:,:,i);        
        SquErr=(Y_U-Y_orig).^2;        

        SSE_G(i,j)=sum(SquErr(:)); 
        RSSE_G(i,j)=sum(SquErr(:))/sum(Y_orig(:));           
        RE_G(i,j)= mean(sqrt(sum(SquErr,1))./sum(Y_orig,1));
        
       
   end
     waitbar(j/Num_Test);
end
close(h);

%% -------------MOR Bases by GPE Predicted snapshot-------------------------
% n_Ubases=5;   %number of POD basics
[Y_GPE,Time_Emu]=Func_DrGPE(X(1:Num_Train,:),Y(1:Num_Train,:),X_star,options,preoptions);     

h = waitbar(0,'Test MOR Bases by GPE snapshot');
for i =1:Num_Test
    
%     uParas.ampx=1;
%     uParas.ampy=1;
    uParas.x0=X_star(i,1);
    uParas.y0=X_star(i,2);

    Y_GPESnapS =reshape(Y_GPE(i,:)',Paras.Num_node_y*Paras.Num_node_x,[]);
    [U_GPESnapS,~,~]=svd(Y_GPESnapS);  % U*S*V'=Rec_X
    
    for j = n_Ubases
        U=U_GPESnapS(:,1:j);
    %     eigenvalues_Ustar=diag(S);
        [Y_U,~]=Us_HeatDC_MOR_SolveFunc(Paras,uParas,X_int,U);
        
        Y_orig=Y_Rec_star(:,:,i);        
        SquErr=(Y_U-Y_orig).^2;        
        
        SSE_K(i,j)=sum(SquErr(:));    
        RSSE_K(i,j)=sum(SquErr(:))/sum(Y_orig(:));     
        RE_K(i,j)= mean(sqrt(sum(SquErr,1))./sum(Y_orig,1));

    end
    waitbar(i/Num_Test);
end
close(h);

% Time_ROM_U_GPESnapS=toc-toc0;

%%
nPOD=n_Ubases;
figure
boxplot(RE_G(:,nPOD))
figure
boxplot(RE_K(:,nPOD))

figure
boxplot(SSE_G(:,nPOD))
figure
boxplot(SSE_K(:,nPOD))

figure
boxplot(RSSE_G(:,nPOD))
figure
boxplot(RSSE_K(:,nPOD))

%%
% nPOD=5:5:50;
nPOD=n_Ubases;

figure
boxplot(RE_G(:,nPOD))
% title(sprintf('SSE Global basis. Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i',Num_Train,Num_Test,Num_Snapshot,n_Ubases))
hold off

xlabel('POD basis dimension');
ylabel('Square error');
title('(b)');

set(gca,'yscale','log');

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=-1;
dB=-6;

ax.YLim=[10^dB,10^uB];
ytick=logspace(dB,uB,uB-dB+1);
ax.YTick=ytick;
ax.XTick=nPOD;
ax.XTick=1:length(nPOD);
ax.XTickLabel=nPOD;

%% 
figure
boxplot(RE_K(:,nPOD))
% title(sprintf('SSE %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{Basis}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,n_Ubases,dim_new))
hold off

xlabel('POD basis dimension');
ylabel('Square error');
title('(a)');
set(gca,'yscale','log');

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=-1;
dB=-6;

ax.YLim=[10^dB,10^uB];
ytick=logspace(dB,uB,uB-dB+1);
ax.YTick=ytick;

ax.XTick=1:2:15;
ax.XTickLabel=1:2:15;

ax.XTick=nPOD;
ax.XTick=1:length(nPOD);
ax.XTickLabel=nPOD;
