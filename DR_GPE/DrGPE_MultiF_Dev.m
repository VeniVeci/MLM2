% DrGPE_MultiF_Dev
% Dimension reduction Gausian process multi-field emulation_Developer version
% Offer an Analysis through data with given data and Method/parameters
% 
% 
% Modifications:
% WeiX, 14-1-2016, Create.

clear
% close all 


%% Dataset Parameter
% index_dataset=5;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR
            % 9: Battery
          
% num_train_ref=[50:50:150]'; % must be in INCREASING order
num_train=120;
num_test=300;

%% Generate dataset

%Multi-field dataset
Test_split=201;
Name_field='exp1';
[X_Data,U,V,P,option_Data] = Exp2UVPv(Name_field);

%Test data
% X_star=X(Test_split:Test_split+num_test-1,:);
% U_starorig=U(Test_split:Test_split+num_test-1,:);
% V_starorig=V(Test_split:Test_split+num_test-1,:);
% P_starorig=P(Test_split:Test_split+num_test-1,:);

% num_train_max=num_train_ref(end);
%Train data

% U=U(1:num_train_max,:);
% V=V(1:num_train_max,:);
% P=P(1:num_train_max,:);

[~,Dim_U]=size(U);
[~,Dim_V]=size(V);
[~,Dim_P]=size(P);

Y_Data=[U,V,P];  

Y=Y_Data(1:num_train,:);
X=X_Data(1:num_train,:);

X_star=X_Data(Test_split:Test_split+num_test-1,:);
Y_starorig=Y_Data(Test_split:Test_split+num_test-1,:);


%% DR method and parameters
dim_new=10;           
% DrMethod='kPCA';
% DrMethod='DiffusionMaps';
% DrMethod='Isomaps';
DrMethod='PCA';

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
preoptions.neighbor=10;
% preoptions.type='LpcaI';
% preoptions.dim_new=5; % Use to stable the result but sacrefy accuracy




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
SSE=zeros(num_test,dim_new);  
    
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

    %Result Analysis&Record
    
    SquErr=(Y_starorig-Y_star).^2;        
%     MSS_orig=mean(Y_starorig.^2,2);
    SSE(:,k)=sum(SquErr,2);
%     RSSE(:,k)=mean(SquErr,2)./MSS_orig;

%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

        %Seperate field           
        U_star=Y_star(:,1:Dim_U);
        V_star=Y_star(:,1+Dim_U:Dim_U+Dim_V);
        P_star=Y_star(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
        
        U_starorig=Y_starorig(:,1:Dim_U);
        V_starorig=Y_starorig(:,1+Dim_U:Dim_U+Dim_V);
        P_starorig=Y_starorig(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
             
        %Error record
%         SquErr=(Y_starorig-Y_star).^2;
        
        SquEr_U=(U_starorig-U_star).^2;
        SquEr_V=(V_starorig-V_star).^2;
        SquEr_P=(P_starorig-P_star).^2;
        
        MSS_Uorig=mean(U_starorig.^2,2);    %Mean Square Sum of U original field
        MSS_Vorig=mean(V_starorig.^2,2);
        MSS_Porig=mean(P_starorig.^2,2);
        
%         mean_Uorig=mean(U_starorig,2);
%         mean_Vorig=mean(V_starorig,2);
%         mean_Porig=mean(P_starorig,2);
        
%         RSSE_U(:,k)=sum(SquEr_U,2)./MSS_Uorig;          %Relativeu Square Sum Error
%         RSSE_V(:,k)=sum(SquEr_V,2)./MSS_Vorig;
%         RSSE_P(:,k)=sum(SquEr_P,2)./MSS_Porig;        
%         FSRSSE(:,k)=RSSE_U(:,k)+RSSE_V(:,k)+RSSE_P(:,k); %Field Sum Relative Square Sum Error
        
        
        RSSE_U(:,k)=mean(SquEr_U,2)./MSS_Uorig;        %Relative Mean Square Sum Error
        RSSE_V(:,k)=mean(SquEr_V,2)./MSS_Vorig;
        RSSE_P(:,k)=mean(SquEr_P,2)./MSS_Porig;        
        MRSSE(:,k)=(RSSE_U(:,k)+RSSE_V(:,k)+RSSE_P(:,k))./3; %Field Sum Relative Mean Square Sum Error

end

%% Result presentation

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
boxplot(MRSSE)
title(sprintf('MRSSE -%s -%s -%d Training points -multiF', DrMethod, preoptions.type, num_train))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Square sum error')
MMRSSE=mean(MRSSE);

% figure
% boxplot(SSE)
% title(sprintf('SSE -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
ylabel('Square sum error')


%% Detail plot


