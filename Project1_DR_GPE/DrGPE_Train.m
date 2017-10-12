% DrGPE_Train
% Prototype I

clear
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

num_train_ref=[40:40:200]'; % must be in INCREASING order

num_test=100;

%%
dim_new=10;           % temp value.

% kernel PCA options----------------------
Koptions.ker='gaussian';
Koptions.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
Koptions.new_dim=dim_new;
Koptions.FullRec=0;

% PreImage options----------------------
preoptions.type='Exp';  %'LSE' OR 'Dw'
% preoptions.para=2;    % These paramether onlys used for 'LSE' and 'Dw'
% preoptions.neighbor=10;


%     %Structure 1
%     % meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
%     % covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
%     % likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

%     %Structure 2
%     covfunc = @covSEiso; 
%     hyp.cov = [0; 0];
% 
%     likfunc = @likGauss; 
%     sn = 0.1;
%     hyp.lik = log(sn);
% 
%     meanfunc=[];
%     hyp.mean=[];

    %Structure 3
    Dim_X=1; % only a temporary value as reminder. Would be detect later by the script.
    covfunc = {@covSum,{@covSEard,@covNoise}}; 
    hyp.cov = [zeros(Dim_X+1,1);0];

    likfunc = @likGauss; 
    sn = 0.1;
    hyp.lik = log(sn);

    meanfunc=[];
    hyp.mean=[];


Z_Rec=zeros(num_train_ref(end),dim_new,length(num_train_ref));
% Z_Rec=[];

h1 = waitbar(0,'Loop1');
for i=1:length(num_train_ref)
    
    waitbar(i/length(num_train_ref),h1)
    
    % Generate dataset
    num_train=num_train_ref(i);
    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    Distance =pdist2(Y,Y,'euclidean');
    Koptions.arg=sum(Distance(:).^2)/(num_train^2);   
    Koptions.arg=sqrt(Koptions.arg/2);
    
    [Z,model] = Kpca(Y,Koptions);    
    
%     Z_Rec=[Z_Rec;Z];
    Z_Rec(1:num_train,1:dim_new,i)=Z;
    model_Rec(i)=model;

    h2 = waitbar(0,'Loop2');
    for j=1:dim_new

        %Clean memory and initilize hyperparameter !IMPORTANT!
%         hyp.cov = [0; 0];sn = 0.1;hyp.lik = log(sn);hyp.mean=[]; 
        hyp.cov = [zeros(Dim_X+1,1);0];
        sn = 0.1;
        hyp.lik = log(sn);
        hyp.mean=[];
            
            
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,j));
%       hyp = minimize(hyp, @gp, -100, @infMCMC,  meanfunc, covfunc, likfunc, X, Z(:,i));
%             exp(hyp.lik)
%             nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i))
%             hypArray(i,j)=hyp;
        hyp_Rec(i,j)=hyp;


    %         %MCMC
    %         hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i));
    %         exp(hyp.lik)
    %         nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i))
    %         [m(:,i) s(:,i)] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
     
%     waitbar(((i-1)*dim_new+j)/(length(num_train_ref)*dim_new),'Training System working very hard')

        
        waitbar(j/dim_new,h2,'Loop2')
    
    end
    close(h2) 
    
end

close(h1) 





