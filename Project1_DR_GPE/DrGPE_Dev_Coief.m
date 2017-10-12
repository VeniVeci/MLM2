% DrGPE_HQ_Train
% A HQ terminal for dimension reduction GPE

%% Initialize
clear
close all 


%% Dataset Parameter
index_dataset=5;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR
% num_train=200;                 

% num_train_ref=[50:50:150]'; % must be in INCREASING order

num_train=200;
num_test=300;
dim_new=10;           


%% DR method and parameters
DrMethod='kPCA';
% DrMethod='DiffusionMaps';

switch DrMethod
    
    case 'kPCA'
        options.ker='gaussian';
        options.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
        % Koptions.arg=500;    % 500 For index_dataset=8
        options.new_dim=dim_new;
        options.FullRec=0;       
    case 'DiffusionMaps'
        options.metric ='euclidean';
        options.kernel ='gaussian'; 
        % Doptions.kpara = 10000;             
        options.kAuto=1;
        options.dim_new = dim_new;              
        options.t = 1;                     
        options.FullRec = 0;      
    case 'Isomaps'
        
    otherwise 
        error('No such DR method')
end

% PreImage options----------------------
preoptions.type='Exp';  %'LSE' OR 'Dw'
% preoptions.para=2;
% preoptions.neighbor=10;



%% GPE structure and parameters
% % Structure 1
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
% likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);


% %Structure 2
% %Isotropic covariance function
% % covfunc = @covSEiso; 
% % hyp.cov = [0; 0];
% 
% covfunc =@covSEard;
% hyp.cov = [];       % will be updated later 
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

% Z_Rec=zeros(num_train_ref(end),dim_new+1,length(num_train_ref));    % +1 just for diffusion maps
% h1 = waitbar(0,'Loop1');
% 
% for i=1:length(num_train_ref)
%     
%     waitbar(i/length(num_train_ref),h1)
    
    % Generate dataset
    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    % Dimension reduction
    switch DrMethod
    
        case 'kPCA'
            % Auto select kernel parameter
            Distance =pdist2(Y,Y,'euclidean');
            options.arg=sum(Distance(:).^2)/(num_train^2);   
            options.arg=sqrt(options.arg/2);
            [Z,model] = Kpca(Y,options);    

        case 'DiffusionMaps'
            [Z,model] = DiffusionMaps(Y,options);
            
        case 'Isomaps'
            
        otherwise 
        error('No such DR method')
    end
    
    
    [num_Z,dim_Z]=size(Z);    
    % Assumened uncorrelated MISO GP on Z
    
    %Rescale
    [Z,ZKey] = DataPP(Z);
    
    for i=1:dim_Z
        [num_X,dim_X]=size(X);
        hyp.cov = [zeros(Dim_X+1,1);0];
        hyp.lik = log(sn);
        hyp.mean=[];
  
        %GP
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));
        exp(hyp.lik)
        nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));        
        [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
    end
    
    %Plot rescale coefficients
     for dim = 1:dim_Z
        mu=m(:,dim);
        sigma=s(:,dim);

        [mu,index]=sort(mu);
        sigma=sigma(index);
        %         figure(i)
        figure(dim)
        title('rescale coiefficients');
        errorbar(mu,sigma,'rx')
     end

    %Rescale recover
    [m] = DataPP_Rocv(m,ZKey);
    [s] = DataPP_Rocv(s,ZKey);
    
    %Clean memory
    SSE=zeros(num_test,dim_new);  
    
    for k=1:dim_new    
        switch DrMethod    
            case 'kPCA'
                Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);
            case 'DiffusionMaps'
                Y_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
            case 'Isomaps'

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
    
%      VarRelaError=abs(s./m);
%      sumVarRelaError(i,:)=sum(VarRelaError);
     
%     figure(10+i)
%     boxplot(VarRelaError)
%     title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
%     set(gca,'yscale','log');
%     xlabel('Index of principal direction')
%     ylabel('Predicted variance ratio')


%Plot coefficients

 for dim = 1:dim_new
    mu=m(:,dim);
    sigma=s(:,dim);

    [mu,index]=sort(mu);
    sigma=sigma(index);
    %         figure(i)
    figure(dim+20)
    title('recover rescale coiefficients');
    errorbar(mu,sigma,'rx')
 end
     
     
figure
boxplot(RSSE)
title(sprintf('Boxplot of Square sum error for different subspace dimension -%s -%d Training points', DrMethod, num_train))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Square sum error')

 MRSSE(i,:)=mean(RSSE);

