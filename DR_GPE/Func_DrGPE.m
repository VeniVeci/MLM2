function [Y_star,time]=Func_DrGPE_Dev(X,Y,X_star,options,preoptions)
% Func_DrGPE_Dev
% Dimension reduction Gausian process emulation function version
% Include most detial of reduced based GPE
%
% Instructions:
% 
% Modifications:
% WeiX, 17-3-2016, First Edition

start_time = cputime;

if nargin < 4, options = []; end

if ~isfield(options,'DrMethod'), options.DrMethod ='PCA'; end
if ~isfield(options,'dim_new'), options.dim_new =10; end

if ~isfield(preoptions,'type'), preoptions.type='Exp';   end
if ~isfield(preoptions,'neighbor'), preoptions.neighbor=10;   end


%% DR method and parameters Auto Selection

% switch DrMethod
%     
%     case 'kPCA'
%         options.ker='gaussian';   
%         options.new_dim=dim_new;
%         options.FullRec=0;       
%         options.arg=1000;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
%         options.kAuto=1;
%         
%     case 'DiffusionMaps'
%         options.metric ='euclidean';
%         options.kernel ='gaussian'; 
%         options.dim_new = dim_new;              
%         options.t = 1;                     
%         options.FullRec = 0;      
%         % Doptions.kpara = 10000;             
%         options.kAuto=1;
%         options.Ztype=0;    %Output type. With/without 1st component
%         
%     case 'Isomaps'
%         options.dim_new=dim_new;                % New dimension
%         options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
%         options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
%         options.metric='euclidean';             % Method of measurement. Metric
% 
%         %Isomap PreImage options
%         preoptions.ReCoverNeighborType='k';     % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
%         preoptions.ReCoverNeighborPara=10;      % Parameter of neighbor of new point
%         
%     case 'PCA'
%         % options=[];
%         
%     otherwise 
%         error('No such DR method')
% end

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


%% Main
%Initialize

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


%% Univariate GPE
%Assumened uncorrelated MISO GP on Z   
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

time = cputime - start_time;
