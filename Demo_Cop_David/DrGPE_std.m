% DrGPE_std
% Dimension reduction GPE with standarization
%
% 
% Modifications:
% WeiX, 7-7-2016, Create.

clear
% close all 



%% Dataset Parameter

num_train=40;
num_test=100;

d_X=10;         %Number of coefficients used for input

% Test_split=701;
% num_train_ref=[50:50:150]'; % must be in INCREASING order

%% Generate dataset
load('X_design');
load('Y_design');
load('X_test');
load('Y_test');

% X_design=X_design';
% Y_design=Y_design';
% X_test=X_test';
% Y_test=Y_test';

X=X_design(1:d_X,1:num_train,:)';
Y=Y_design(:,1:num_train)';

X_star=X_test(1:d_X,1:num_test)';
Y_starorig=Y_test(:,1:num_test)';


%% PCA
[Z1,model1] = Lpca(Y,40,[]);



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
        options.arg0=1;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
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
% preoptions.type='Dw';
% preoptions.para = 1;
% preoptions.type='LpcaI';
% preoptions.dim_new=5; % Use to stable the result but sacrefy accuracy
% preoptions.neighbor=num_train/10;
preoptions.neighbor=5;


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

% Structure 4 Polynomial
% 
% % cp = {@covPoly ,3}; c = 2; hypp = log([c;sf]);
% 
% covfunc={@covPoly ,3};
% hyp.cov = log([2;2]);
% 
% likfunc = @likGauss; 
% sn = 0.1;
% hyp.lik = log(sn);
% 
% meanfunc=[];
% hyp.mean=[];

%% Main
%Initialize

%% Dimension reduction
switch DrMethod

    case 'kPCA'
        [Z,model] = Kpca(Z1,options);    
    case 'DiffusionMaps'
        [Z,model] = DiffusionMaps(Z1,options);   
    case 'Isomaps'
        [Z,model] = Isomaps(Z1,options);
    case 'PCA'
        [Z,model] = Lpca(Z1,dim_new,options);
    otherwise 
        error('No such DR method')
end



[np_train,Dim_X]=size(X);

%% Univariate GPE
%Assumened uncorrelated MISO GP on Z   

[num_Z,dim_Z]=size(Z);     
for i=1:dim_Z
    
    % Clean memory of hyperparameter
    [num_X,dim_X]=size(X);

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
SSE=zeros(num_test,dim_new);  
    
for k=1:dim_new 
    switch DrMethod    
        case 'kPCA'
            Z1_star = Kpca_PreImage(m(:,1:k),model,preoptions);
        case 'DiffusionMaps'   
            switch options.Ztype
                case 0
                    Z1_star = DiffusionMaps_PreImage(m(:,1:k),model,preoptions);
                case 1
                    Z1_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
            end            
            
        case 'Isomaps'
            Z1_star= Isomaps_PreImage(m(:,1:k),model,preoptions);            
        case 'PCA'
            Z1_star=Lpca_PreImage(m(:,1:k),model);
            
        otherwise 
            error('No such DR method')
    end
    
     Y_star=Lpca_PreImage(Z1_star,model1);
    
    SquErr=(Y_starorig-Y_star).^2;        
    MSS_orig=mean(Y_starorig.^2,2);

    SSE(:,k)=sum(SquErr,2);
    RE(:,k)=mean(SquErr,2)./MSS_orig;

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
boxplot(RE)
title(sprintf('Error rate -%s -%d Training points d_x=%d', DrMethod, num_train,d_X))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Error rate')

MRSSE=mean(RE);


% figure
% boxplot(SSE)
% title(sprintf('SSE -%s -%d Training points d_x=%d', DrMethod, num_train,d_X))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Error rate')
% 




% figure
% plot([MRE_U',MRE_V',MRE_P'],'-o');
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Error rate')
% legend('U','V','P');
% title(sprintf('Rate Error -%s -%d Training points', DrMethod, num_train))

% figure
% boxplot(RE_P)
% title(sprintf('Error rate P -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Error rate')

%Figure SSE
% figure
% boxplot(SSE)
% title(sprintf('SSE -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Square sum error')


%% Detail plot
% index=1;
% figure
% field=reshape(Y_star(index,:),[25,25]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf('%s GPE Predicted field -%d Training points x1=%3.2f x2= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))
% 
% figure
% field=reshape(Y_starorig(index,:),[25,25]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf(' Real field x1=%3.2f x2= %3.3f ',X_star(index,1),X_star(index,2)))

