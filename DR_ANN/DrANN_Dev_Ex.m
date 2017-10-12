% Dimension reduction Neural Network Emulator_Developer version_Expanded 
% Include most detial of reduced based ANN with many functions.
%
% Instructions:
% 
% Modifications:
% WeiX, lost date,    first edition 
% WeiX, 5-1-2016, Minor Update

clear
% close all 


%% Dataset Parameter
index_dataset=9;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR
          
% num_train_ref=[50:50:150]'; % must be in INCREASING order
num_train=200;
num_test=300;


%% DR method and parameters
dim_new=10;           
DrMethod='kPCA';
% DrMethod='DiffusionMaps';
% DrMethod='Isomaps';


switch DrMethod
    
    case 'kPCA'
        options.ker='gaussian';   
        options.new_dim=dim_new;
        options.FullRec=0;       
        options.arg=100;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
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

    otherwise 
        error('No such DR method')
end

% PreImage options----------------------
preoptions.type='Exp';  %'LSE', 'Dw' OR 'Exp'
% preoptions.para=2;
% preoptions.neighbor=10;


%% Main
%Initialize

%% Generate dataset
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

%% Dimension reduction
switch DrMethod

    case 'kPCA'
        [Z,model] = Kpca(Y,options);    
    case 'DiffusionMaps'
        [Z,model] = DiffusionMaps(Y,options);   
    case 'Isomaps'
        [Z,model] = Isomaps(Y,options);
    otherwise 
        error('No such DR method')
end



%%ANN structure
MultiVar=1;     %MultiVariate ANN or uncorrelated structure

[num_Z,dim_Z]=size(Z);

switch MultiVar
    case 0
        %% Uncorrelated univariate ANN
        for i=1:dim_Z    
            net=feedforwardnet(10); % One hidden layer with nn nodes; for more layers, 
            % use [nn1 nn2 nn3 ... nnJ] for J layers with nnj nodes in the jth layer 
            net = init(net); % Reinitialize weights
            net.divideParam.trainRatio=0.9; % Fraction of data used for training (cross-validation)
            net.divideParam.valRatio=(1-net.divideParam.trainRatio);% /2; % Fraction of data used for validation
            net.divideParam.testRatio=0;% (1-net.divideParam.trainRatio)/2; % Fraction of data used for testing
            % [net,tr] = trainlm(net,X,Y); % Feedforward with Levenberg-Marquardt backpropagation
            [net,tr] = trainbr(net,X',Z(:,i)'); % Bayesian Regulization

            for j=1:num_test        
                Z_star(i,j)=net(X_star(j,:)');    
            end
        end
        Z_star=Z_star';
        
    case 1
        %% Multi-variate ANN
        net=feedforwardnet(10); % One hidden layer with nn nodes; for more layers, 
        % use [nn1 nn2 nn3 ... nnJ] for J layers with nnj nodes in the jth layer 
        net = init(net); % Reinitialize weights
        net.divideParam.trainRatio=0.9; % Fraction of data used for training (cross-validation)
        net.divideParam.valRatio=(1-net.divideParam.trainRatio);% /2; % Fraction of data used for validation
        net.divideParam.testRatio=0;% (1-net.divideParam.trainRatio)/2; % Fraction of data used for testing
        % [net,tr] = trainlm(net,X,Y); % Feedforward with Levenberg-Marquardt backpropagation
        [net,tr] = trainbr(net,X',Z'); % Bayesian Regulization
        
        for i=1:num_test        
            Z_star(:,i)=net(X_star(i,:)');    
        end
        Z_star=Z_star';
        
    otherwise 
        error('Error ANN MultiVar structure')
        
end


%% Pre-image
%Clean memory
SSE=zeros(num_test,dim_new);  
    
for k=1:dim_new    
    switch DrMethod    
        case 'kPCA'
            Y_star = Kpca_PreImage(Z_star(:,1:k),model,preoptions);
        case 'DiffusionMaps'   
            switch options.Ztype
                case 0
                    Y_star = DiffusionMaps_PreImage(Z_star(:,1:k),model,preoptions);
                case 1
                    Y_star = DiffusionMaps_PreImage(Z_star(:,1:k+1),model,preoptions);
            end            
            
        case 'Isomaps'
            Y_star= Isomaps_PreImage(Z_star(:,1:k),model,preoptions);
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
%      
%Relatice square sum error     
figure
boxplot(RSSE)
title(sprintf('Boxplot of RSSE for different subspace dimension -%s -%d Training points', DrMethod, num_train))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Square sum error')

MRSSE=mean(RSSE);

