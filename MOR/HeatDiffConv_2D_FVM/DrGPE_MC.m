% DrGPE_MC
% Dimension reduction Gausian process emulator on data Multiple cases.
% Used to find the proper data set for demonstration (parameters in proper range).
% 
% Modifications:
% WeiX, 5-4-2016, Create.



clear
% close all 


%% ----------------Load dataset-------------------------------------------
Num_Train=80;   
Num_Test=300;
Test_StartIndex=201;
Num_Snapshot=10;

load('HeatDC_DBC_ExpData13.mat') %7 is ok for HH, HH fails in 4.
% Num_Trian=40;   % 80 is best for 'ExpData3,4.mat'
% Index_test=404; % 400,403 is tracky. mostly fine. in ExpData4 404 best

% X_star=[10,-10,5];
Paras.t_n=Paras.t_end/Paras.dt;
Index_snapshot=1:(Paras.t_n/Num_Snapshot):Paras.t_n+1;     %Shoube be adjusted accroding to Paras.dt & Paras.t_end;

Y=Y_Rec(:,Index_snapshot,:);
Y=reshape(Y,[],500); %500 for this dataset
Y=Y';


% Time_HDM=Time_Rec;

X_star=X(Test_StartIndex:Test_StartIndex+Num_Test-1,:);
Y_starorig=Y(Test_StartIndex:Test_StartIndex+Num_Test-1,:);

X=X(1:Num_Train,:);
Y=Y(1:Num_Train,:);


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
        options.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
        options.kAuto=1;
        
    case 'DiffusionMaps'
        options.metric ='euclidean';
        options.kernel ='gaussian'; 
        options.dim_new = dim_new;              
        options.t = 1;                     
        options.FullRec = 0;      
        options.kpara = 10;             
%         options.kAuto=0;
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

% PreImage options----------------------
preoptions.type='Exp';  %'LSE', 'Dw' OR 'Exp'
% preoptions.type='LpcaI';
% preoptions.dim_new=5; % Use to stable the result but sacrifice accuracy
preoptions.neighbor=10;


%% GPE structure and parameters

% %Structure 1 
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
% likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);


% %Structure 2 %Isotropic covariance function (simple)
% % covfunc = @covSEiso; 
% % hyp.cov = [0; 0];
% 
% likfunc = @likGauss; 
% sn = 0.1;
% hyp.lik = log(sn);
% 
% meanfunc=[];
% hyp.mean=[];


%Structure 3
Dim_X=1; % only a temporary value as reminder. Would be detect later by the script.
covfunc = {@covSum,{@covSEard,@covNoise}}; 
hyp.cov = [zeros(Dim_X+1,1);0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];


InfMethod=@infExact;
% InfMethod=@infMCMC;


%% Main
%Initialize

%% Dimension reduction
switch DrMethod

    case 'kPCA'
        [Z,model] = Kpca(Y,options);    
    case 'DiffusionMaps'
        [Z,model] = DiffusionMaps(Y,options);   
    case 'Isomaps'
        [Z,model] = Isomaps(Y,options);
    case 'PCA'
        [Z,model] = Lpca(Y,dim_new,options);
    otherwise 
        error('No such DR method')
end



[np_train,Dim_X]=size(X);

%% Univariate GPE
%Assumened uncorrelated MISO GP on Z   
[num_Z,dim_Z]=size(Z);     
for i=1:dim_Z
    
    %Clean memory of hyperparameter
%     [num_X,dim_X]=size(X);
    hyp.cov = [zeros(Dim_X+1,1);0];
    sn = 0.1;
    hyp.lik = log(sn);
    hyp.mean=[];

    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));
    exp(hyp.lik);
    nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));        
    [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
end
    
%% Pre-image
%Clean memory
SSE=zeros(Num_Test,dim_new);  
    
h = waitbar(0,'Running');
for k=1:dim_new 
    switch DrMethod    
        case 'kPCA'
            Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);
        case 'DiffusionMaps'   
            switch options.Ztype
                case 0
                    Y_star = DiffusionMaps_PreImage(m(:,1:k),model,preoptions);
                case 1
                    Y_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
            end            
            
        case 'Isomaps'
            Y_star= Isomaps_PreImage(m(:,1:k),model,preoptions);            
        case 'PCA'
            Y_star=Lpca_PreImage(m(:,1:k),model);
            
        otherwise 
            error('No such DR method')
    end

    SquErr=(Y_starorig-Y_star).^2;        
    MSS_orig=mean(Y_starorig.^2,2);

    SSE(:,k)=sum(SquErr,2);
    RSSE(:,k)=mean(SquErr,2)./MSS_orig;
    
    
%     Y=Y_Rec(:,Index_snapshot,:);
%     Y=reshape(Y,[],500); %500 for this dataset
    
    %Relative error
    for i=1:Num_Test
        iY=Y_star(i,:);
        iFieldY=reshape(iY,2500,length(Index_snapshot));
        iFieldYOrig=Y_Rec(:,Index_snapshot,Test_StartIndex+i-1);
        RE(i,k)= mean (sqrt(sum((iFieldY-iFieldYOrig).^2,1) ./ sum(iFieldYOrig.^2,1)) ) ;
        RE2(i,k)= mean ((sum((iFieldY-iFieldYOrig).^2,1) ./ sum(iFieldYOrig.^2,1)) ) ;
        
    end
    
waitbar(k/dim_new);
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

end
close(h);
%% Result presentation

% %Ratio of variance.mean of coefficients
% VarRelaError=abs(s./m);
% sumVarRelaError(i,:)=sum(VarRelaError);     
% figure(10+i)
% boxplot(VarRelaError)
% title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
% set(gca,'yscale','log');
% xlabel('Index of principal direction')
% ylabel('Predicted variance ratio')

% %Predictions of coefficients
% for i = 1:dim_new
%     %Extract
%     mu=m(:,i);
%     sigma=s(:,i);
%     %Sort
%     [mu,index]=sort(mu);
%     sigma=sigma(index);
%     %Plot
%     figure(i)
%     errorbar(mu,sigma,'rx')
% end
     
%% Relatice square sum error     
%Figure RSSE
figure
boxplot(RE)
title(sprintf('RE. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,dim_new))
set(gca,'yscale','log');
xlabel('Approximate manifold dimension')
ylabel('Relative error')

title('(b)')


figure
boxplot(RE2)
title(sprintf('RE. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,dim_new))
set(gca,'yscale','log');
xlabel('Approximate manifold dimension')
ylabel('Relative error')
title('(a)')

ax = gca; 
ax.YScale='log';
%Range for log scale
ax = gca; 
uB=1;
dB=-5;

ax.YLim=[10^dB,10^uB];
ytick=logspace(dB,uB,uB-dB+1);

ax.YTick=ytick(1:2:end);




%% Analysis for representative case 
mColumn=15;
data=RE(:,mColumn);

meanData=mean(data);
minData=min(data);
maxData=max(data);
quantileData=quantile(data,[.25 .5 .75]);


CompareValue=1.5*quantile(data,[.75]); %Upper whisker
%  CompareValue=meanData;
% CompareValue=maxData;


sortData=(data-CompareValue).^2;
[sortData,sortData_Index]=sort(sortData);
SortX_star=X_star(sortData_Index,:);

RE=data(sortData_Index(1))
xStar=X_star(sortData_Index(1),:)






% figure
% boxplot(RSSE)
% title(sprintf('RSSE. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,dim_new))
% set(gca,'yscale','log');
% xlabel('Approximate manifold dimension')
% ylabel('Relative error')
% 
% title('(a)')

% ax = gca; 
% ax.YScale='log';
% %Range for log scale
% ax = gca; 
% uB=4;
% dB=0;
% 
% ax.YLim=[10^dB,10^uB];
% ytick=logspace(dB,uB,uB-dB+1);
% 
% ax.YTick=ytick(1:2:end);


% MRSSE=mean(RSSE);
% 
% figure
% boxplot(SSE)
% title(sprintf('SSE. %s, Num_{Trian}=%0i, Num_{Test}=%0i, Num_{Snapshot}=%0i, Num_{DimNew}=%i', DrMethod,Num_Train,Num_Test,Num_Snapshot,dim_new))
%  set(gca,'yscale','log');
% xlabel('Approximate manifold dimension')
% ylabel('Square error')
% 
% title('(b)')

%Figure SSE
% figure
% boxplot(SSE)
% title(sprintf('SSE -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Square sum error')

% ax = gca; 
% ax.YScale='log';
% %Range for log scale
% ax = gca; 
% uB=6;
% dB=0;
% 
% ax.YLim=[10^dB,10^uB];
% ytick=logspace(dB,uB,uB-dB+1);
