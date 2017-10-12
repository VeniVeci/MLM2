% DrGPE_SS
% Dimension reduction Gausian process on snapshots. This script is used to
% find the proper data set for demonstration (parameters in proper range).
% 
% Modifications:
% WeiX, 21-3-2016, Create.




clear
% close all 


%% ----------------Load dataset-------------------------------------------
load('ExpData8.mat') %7 is ok for HH,
Num_Train=150;   
Num_Test=200;
Test_StartIndex=301;

% X_star=[10,-10,5];
X_Data=X;
Y_Data=Y;

X=X_Data(1:Num_Train,:);
Y=Y_Data(1:Num_Train,:);

X_star=X_Data(Test_StartIndex:Test_StartIndex+Num_Test-1,:);
Y_starorig=Y_Data(Test_StartIndex:Test_StartIndex+Num_Test-1,:);


%% DR method and parameters
dim_new=10;           
DrMethod='kPCA';
% DrMethod='DiffusionMaps';
DrMethod='Isomaps';
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

%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

end

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
boxplot(RSSE)
title(sprintf('RSSE -%s -%d Training points', DrMethod, Num_Train))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Square sum error')

MRSSE=mean(RSSE);

%Figure SSE
% figure
% boxplot(SSE)
% title(sprintf('SSE -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Square sum error')

