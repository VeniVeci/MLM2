% DrGPE_HQ_Train_NS
% A HQ terminal for dimension reduction GPE _Special design for Navier
% stock equation for the special format of dataset( e.g. using U,V,P as Y
% which are multi-output.)

%% Initialize
clear
close all 


%% Dataset Parameter
Name_field='exp1';

[X,U,V,P,option_Data] = Exp2UVPv(Name_field);

Y_Data=P;
X_Data=X;

clearvars -except Y_Data X_Data

num_train_ref=[40:40:200]'; % must be in INCREASING order
% num_train_ref=[100:100:400]'; % must be in INCREASING order

num_test=100;
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
%  preoptions.neighbor=10;



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
% Initialize
Z_Rec=zeros(num_train_ref(end),dim_new+1,length(num_train_ref));    % +1 just for diffusion maps
h1 = waitbar(0,'Loop1');

for i=1:length(num_train_ref)
    
    waitbar(i/length(num_train_ref),h1)
    
    % Generate dataset
    num_train=num_train_ref(i);
    
    Y=Y_Data(1:num_train,:);
    X=X_Data(1:num_train,:);
    
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
    Z_Rec(1:num_Z,1:dim_Z,i)=Z;
    model_Rec(i)=model;
    
    % Assumened uncorrelated MISO GP on Z
    h2 = waitbar(0,'Loop2');
    for j=1:dim_Z
        
        %Initialize GPE 
        %Clean memory. !IMPORTANT! Change it with structure. 
        
%         [num_X,dim_X]=size(X);
%         hyp.cov = ones(dim_X+1,1);
%         sn = 0.1;
%         hyp.lik = log(sn);
%         hyp.mean=[];     
        
        [num_X,dim_X]=size(X);
        hyp.cov = [zeros(dim_X+1,1);0];
        hyp.lik = log(sn);
        hyp.mean=[];

        
        
        hyp = minimize(hyp, @gp, -100, InfMethod, meanfunc, covfunc, likfunc, X, Z(:,j));
%             hyp = minimize(hyp, @gp, -100, @infMCMC,  meanfunc, covfunc, likfunc, X, Z(:,i));
%             exp(hyp.lik)
%             nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i))
%             hypArray(i,j)=hyp;
        hyp_Rec(i,j)=hyp;

%       waitbar(((i-1)*dim_new+j)/(length(num_train_ref)*dim_new),'Training System working very hard')
    
        waitbar(j/dim_Z,h2,'Loop2')
    
    end
    close(h2) 
    
end

close(h1) 



% str=func2str(InfMethod);
% filename = sprintf('Data%d%s%s%s',index_dataset,DrMethod,str);
% save(filename)


