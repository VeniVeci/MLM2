% GP_Anal
%
% 
% Modifications:
% WeiX, 28-6-2016, Create.

clear
% close all 

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
            % 9: Battery
          
% num_train_ref=[50:50:150]'; % must be in INCREASING order
num_train=100;
num_test=100;
Test_split=201;


%% Generate dataset

%Uni-field dataset
% [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
% [np_train,Dim_X]=size(X);
% [np_test, Dim_Y]=size(Y_starorig);

% load('exp1');
% load('Exp_v2_3');
load('Exp_v2_5in_v2');
% load('Exp_v2_9in_2');
[n_data,~]=size(X);


U=reshape(RecU,[],n_data)';
V=reshape(RecV,[],n_data)';
P=reshape(RecP,[],n_data)';

Y_Data=[U,V,P];  
X_Data=X;



Y=Y_Data(1:num_train,:);
X=X_Data(1:num_train,:);

X_star=X_Data(Test_split:Test_split+num_test-1,:);
Y_starorig=Y_Data(Test_split:Test_split+num_test-1,:);

[~,Dim_U]=size(U);
[~,Dim_V]=size(V);
[~,Dim_P]=size(P);



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
% preoptions.dim_new=5; % Use to stable the result but sacrefy accuracy
preoptions.neighbor=25;



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

    U_star=Y_star(:,1:Dim_U);
    V_star=Y_star(:,1+Dim_U:Dim_U+Dim_V);
    P_star=Y_star(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);

    U_starorig=Y_starorig(:,1:Dim_U);
    V_starorig=Y_starorig(:,1+Dim_U:Dim_U+Dim_V);
    P_starorig=Y_starorig(:,1+Dim_U+Dim_V:Dim_U+Dim_V+Dim_P);
    
    SquEr_U=(U_starorig-U_star).^2;
    SquEr_V=(V_starorig-V_star).^2;
    SquEr_P=(P_starorig-P_star).^2;

    RE_U(:,k)=sum(SquEr_U,2)./sum(U_starorig.^2,2);
    RE_V(:,k)=sum(SquEr_V,2)./sum(V_starorig.^2,2);
    RE_P(:,k)=sum(SquEr_P,2)./sum(P_starorig.^2,2);
    
end
    
    MRE=RE_U+RE_V+RE_P/3;

    MRE_U=mean(RE_U);
    MRE_V=mean(RE_V);
    MRE_P=mean(RE_P);
    
    



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
% figure
% boxplot(MRE)
% title(sprintf('Mean Error rate -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Error rate')
% MRSSE=mean(RSSE);

figure
plot([MRE_U',MRE_V',MRE_P'],'-o');
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Error rate')
legend('U','V','P');
title(sprintf('Rate Error -%s -%d Training points', DrMethod, num_train))

figure
boxplot(RE_P)
title(sprintf('Error rate P -%s -%d Training points', DrMethod, num_train))
set(gca,'yscale','log');
xlabel('Dimension of manifold')
ylabel('Error rate')

%Figure SSE
% figure
% boxplot(SSE)
% title(sprintf('SSE -%s -%d Training points', DrMethod, num_train))
% set(gca,'yscale','log');
% xlabel('Dimension of manifold')
% ylabel('Square sum error')


%% Detail plot
% index=24;
% figure
% field=reshape(Y_star(index,:),[100,100]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf('%s GPE Predicted field -%d Training points x1=%3.2f x2= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))
% 
% figure
% field=reshape(Y_starorig(index,:),[100,100]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf(' Real field x1=%3.2f x2= %3.3f ',X_star(index,1),X_star(index,2)))

% index=24;
% figure
% field=reshape(P_star(index,:),[50,50]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf('%s GPE Predicted field -%d Training points x1=%3.2f x2= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))
% 
% figure
% field=reshape(P_starorig(index,:),[50,50]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf(' Real field x1=%3.2f x2= %3.3f ',X_star(index,1),X_star(index,2)))

BC.uN=0;     BC.vN=0;
BC.uS=0;     BC.vS=0;
BC.uW=0;     BC.vW=0;
BC.uE=0;     BC.vE=0;


index=3;
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




% figure
% field=reshape(P_star(index,:),[nx-1,ny]);
% surf(field)
% contourf(field,8)
% colorbar
% title(sprintf('%s GPE Predicted field -%d Training points x1=%3.2f x2= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))
% 
% 
% DisplayUV( U,V,P,Re,tf,dt,lx,ly,nx,ny )
