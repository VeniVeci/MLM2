% Dimension reduction Gausian process emulation_Train
% Trian GPE and save trained parameter and result.
%
% Instructions:
% 
% Modifications:
% WeiX, 5-1-2016, Create

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
          
num_train_ref=[40:40:80]'; % must be in INCREASING order
% num_train=80;
num_test=100;


%% DR method and parameters
dim_new=10;           
% DrMethod='kPCA';
% DrMethod='DiffusionMaps';
DrMethod='Isomaps';


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
% Initialize
Z_Rec=zeros(num_train_ref(end),dim_new,length(num_train_ref));    % +1 just for diffusion maps
h1 = waitbar(0,'Loop1');

for i=1:length(num_train_ref)
    
    %Process meter
    waitbar(i/length(num_train_ref),h1)
    
    %% Generate dataset
    num_train=num_train_ref(i);
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
    
    
    [num_Z,dim_Z]=size(Z);    
    Z_Rec(1:num_Z,1:dim_Z,i)=Z;
    model_Rec(i)=model;
    
    % Assumened uncorrelated MISO GP on Z
    h2 = waitbar(0,'Loop2');
    for j=1:dim_Z
        
        %Initialize GPE 
        %Clean memory. !IMPORTANT! Change it with structure. 
        
        [num_X,dim_X]=size(X);
        hyp.cov = [zeros(dim_X+1,1);0];
        sn = 0.1;
        hyp.lik = log(sn);
        hyp.mean=[];

        hyp = minimize(hyp, @gp, -100, InfMethod, meanfunc, covfunc, likfunc, X, Z(:,j));
%         exp(hyp.lik);
%         nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i))
%         hypArray(i,j)=hyp;
        
        hyp_Rec(i,j)=hyp;
  
        waitbar(j/dim_Z,h2,'Loop2')
    
    end
    close(h2) 
    
end

close(h1) 

% str=func2str(InfMethod);
% filename = sprintf('Data%d%s%s%s',index_dataset,DrMethod,str);
% save(filename)


