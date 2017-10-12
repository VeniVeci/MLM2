% GP_Anal3
%
% 
% Modifications:
% WeiX, 28-6-2016, Create.

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
num_train=300;
num_test=100;
Test_split=901;

%% Generate dataset

%Uni-field dataset
% [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
% [np_train,Dim_X]=size(X);
% [np_test, Dim_Y]=size(Y_starorig);

% load('exp1');
% load('Exp_v2_3');
% load('Exp_v2_13in');
% load('Exp_v2_13in_4');
load('Exp_v2_5in_v2');
[n_data,~]=size(X);

U=reshape(RecU,[],n_data)';
V=reshape(RecV,[],n_data)';
P=reshape(RecP,[],n_data)';

U_starorig=U(Test_split:Test_split+num_test-1,:);
V_starorig=V(Test_split:Test_split+num_test-1,:);
P_starorig=P(Test_split:Test_split+num_test-1,:);
X_star=X(Test_split:Test_split+num_test-1,:);

U=U(1:num_train,:);
V=V(1:num_train,:);
P=P(1:num_train,:);
X=X(1:num_train,:);
% X_star=X_Data(Test_split:Test_split+num_test-1,:);


%% GPE structure and parameters

Dim_X=1; % only a temporary value as reminder. Would be detect later by the script.
covfunc = {@covSum,{@covSEard,@covNoise}}; 
hyp.cov = [zeros(Dim_X+1,1);0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];


InfMethod=@infExact;



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
        options.arg=100;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
        options.kAuto=0;
        
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
        options.neighborPara=round(num_train/10);                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
        options.metric='euclidean';             % Method of measurement. Metric

        %Isomap PreImage options
        preoptions.ReCoverNeighborType='k';     % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
        preoptions.ReCoverNeighborPara=round(num_train/10);      % Parameter of neighbor of new point
        
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
preoptions.neighbor=round(num_train/10);
% preoptions.neighbor=10;

%% P Field
    Y=P;
    Y_starorig=P_starorig;

    % Dimension reduction
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

    % Univariate GPE
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
    

    % Pre-image
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

        SquErr=(Y_starorig-Y_star).^2;        
        MSS_orig=mean(Y_starorig.^2,2);

        SSE(:,k)=sum(SquErr,2);
        RE(:,k)=mean(SquErr,2)./MSS_orig;

    %         SE_D(:,i)=SquErr_D(:);
    %         SE_K(:,i)=SquErr_K(:);

    end
    
    RE_P=RE;
    P_star=Y_star;

%% U Field
    Y=U;
    Y_starorig=U_starorig;

    % Dimension reduction
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

    % Univariate GPE
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
    

    % Pre-image
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

        SquErr=(Y_starorig-Y_star).^2;        
        MSS_orig=mean(Y_starorig.^2,2);

        SSE(:,k)=sum(SquErr,2);
        RE(:,k)=mean(SquErr,2)./MSS_orig;

    %         SE_D(:,i)=SquErr_D(:);
    %         SE_K(:,i)=SquErr_K(:);

    end
    
    RE_U=RE;
    U_star=Y_star;    
    
    
 %% V Field
    Y=V;
    Y_starorig=V_starorig;

    % Dimension reduction
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

    % Univariate GPE
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
    

    % Pre-image
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

        SquErr=(Y_starorig-Y_star).^2;        
        MSS_orig=mean(Y_starorig.^2,2);

        SSE(:,k)=sum(SquErr,2);
        RE(:,k)=mean(SquErr,2)./MSS_orig;

    %         SE_D(:,i)=SquErr_D(:);
    %         SE_K(:,i)=SquErr_K(:);

    end
    
    RE_V=RE;
    V_star=Y_star;       
    
 %% Error analysis   
MRE_P=mean(RE_P);
MRE_U=mean(RE_U);
MRE_V=mean(RE_V);


    
    
    
%% Plot   
figure
boxplot(RE_P)
title(sprintf('P Error rate -%s -%d Training points', DrMethod, num_train))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Error rate')

% figure
% boxplot(RE_U)
% title(sprintf('U Error rate -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Error rate')
% 
% figure
% boxplot(RE_V)
% title(sprintf('V Error rate -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Error rate')


figure
plot([MRE_U',MRE_V',MRE_P'],'-o');
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Error rate')
legend('U','V','P');
title(sprintf('Rate Error -%s -%d Training points', DrMethod, num_train))

    
% index=8;
index=8;
field_U=reshape(U_starorig(index,:),nx-1,ny);
field_V=reshape(V_starorig(index,:),nx,ny-1);
field_P=reshape(P_starorig(index,:),nx,ny);
% mini_P=min(field_P(:));
% field_P=field_P-mini_P
figure
DisplayUVP_v2( field_U,field_V,field_P,Re,tf,dt,lx,ly,nx,ny,BC)
title(sprintf('Real field Index=%d',Test_split+index-1));

field_U=reshape(U_star(index,:),nx-1,ny);
field_V=reshape(V_star(index,:),nx,ny-1);
field_P=reshape(P_star(index,:),nx,ny);
% field_P=field_P-mini_P
figure
DisplayUVP_v2( field_U,field_V,field_P,Re,tf,dt,lx,ly,nx,ny,BC )
title(sprintf('%s GPE Predicted field -%d Training points Index=%d ', DrMethod,num_train,Test_split+index-1));



