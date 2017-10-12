%% gpTest
% Using DR GP only to test data

% Modifications:
% 4-2-2017, WeiX, first edition 

clear

%%---------------Setting parameter---------------------------------------
Num_Train=120;   
Num_Test=80;
Test_StartIndex=201;

Num_Snapshot=50;

dim_new=10;      % For LPCA more than 10 require. The new dimension of reduced emulatior

%% DR method and parameters
     
DrMethod='kPCA';
% DrMethod='DiffusionMaps';
DrMethod='Isomaps';
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
preoptions.neighbor=5;
% preoptions.type='LpcaI';
% preoptions.dim_new=10; % Use to stable the result but sacrefy accuracy
% preoptions.InMethod='ANN';

%%
%% DR method and parameters Auto Selection



%% ----------------Load dataset-------------------------------------------
load('HeatDC_DBC_ExpData13.mat')


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


X_star=X;
Y_starorig=Y;

%% -------------MOR Bases by GPE Predicted snapshot-------------------------
% % n_Ubases=5;   %number of POD basics
% [Y_GPE,Time_Emu]=Func_DrGPE(X,Y,X,options,preoptions);     
% error_train=Y_GPE-Y;
% MSE_train =sum(error_train.^2,2);
% 
% 
% [Y_GPE,Time_Emu]=Func_DrGPE(X,Y,X_star,options,preoptions);     
% error_test=Y_GPE-Y_starorig;
% MSE_test  =sum(error_test.^2,2);



%%




%% GPE structure and parameters

%Structure 3
Dim_X=1; % only a temporary value as reminder. Would be detect later by the code.
covfunc = {@covSum,{@covSEard,@covNoise}}; 
hyp.cov = [zeros(Dim_X+1,1);0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];

InfMethod=@infExact;
% InfMethod=@infMCMC;

%% Dimension reduction
switch options.DrMethod

    case 'kPCA'
        [Z,model] = Kpca(Y,options);    
    case 'DiffusionMaps'
        [Z,model] = DiffusionMaps(Y,options);   
    case 'Isomaps'
        [Z,model] = Isomaps(Y,options);
    case 'PCA'
        dim_new=options.dim_new;
        [Z,model] = Lpca(Y,dim_new,options);
    otherwise 
        error('No such DR method')
end


[np_train,Dim_X]=size(X);
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

k=dim_Z;
switch options.DrMethod    
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

%% Error

error_test=Y_star-Y_starorig;
MSE_test  =sum(error_test.^2,2);


RE_test=MSE_test./sum(Y_starorig.^2,2);

RE=(Y_star-Y_starorig).^2 ./ Y_starorig.^2;
RE=mean(RE,2);





