% GPE_Train_Terminal_v1
% The script is used to trina GP given specific dataset and parameters,
% e.g., number of training, test, preserved order,

clear
close all

%% Dataset Parameter
index_dataset=8;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR

num_train=40;                 
num_test=10;     % Wont be used 

% num_train_start=10;    
% num_train_step =10;
% num_train_end  =100;
num_train_ref=[200:50:300]'; % must be in INCREASING order

dim_new_end=10;             % dim_new must start from 1 with step 1. thus no refference array is used
% dim_new_ref=[1:10]';        % must be in INCREASING order. 

%% GP STRUCTURE
% %Structure 1
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
% covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
% likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

%Structure 2
covfunc = @covSEiso; 
hyp.cov = [0; 0];

likfunc = @likGauss; 
sn = 0.1;
hyp.lik = log(sn);

meanfunc=[];
hyp.mean=[];








%------------------------------------------------------------------------------------
%Defining the input
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train_ref(end),num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);

%% PCA


for j=1:length(num_train_ref)
       
    [Z,model] = Lpca(Y(1:num_train_ref(j),:),dim_new_end);
    
        for i=1:dim_new_end; 
            hyp.cov = [0; 0];sn = 0.1;hyp.lik = log(sn);hyp.mean=[]; % Clean memory. !IMPORTANT!
            hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i));
            exp(hyp.lik)
            nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i))
            
    

    %         %MCMC
    %         hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i));
    %         exp(hyp.lik)
    %         nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i))
    %         [m(:,i) s(:,i)] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
    
    
            % Test and plot
            [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X(1:num_train_ref(j),:), Z(1:num_train_ref(j),i), X_star);
            mup=m(:,i)+2*sqrt(s(:,i));
            mdown=m(:,i)-2*sqrt(s(:,i));
            figure(i)
            plot(m(:,i),'-');
            hold on
            plot(mup,'--');
            plot(mdown,'--');
            hold off

%             xScale=[1:num_test]';
%             f = [m(:,i)+2*sqrt(s(:,i)); flipdim(m(:,i)-2*sqrt(s(:,i)),1)];
%             fill([xScale; flipdim(xScale,1)], f, [7 7 7]/8)
               
            hypArray(i,j)=hyp;
%             likPlane(i,j)=exp(hyp.lik);
%             nlml2lane(i,j)=nlml2;
            
        end


    
end

% figure
% plot(likPlane');
% 
% figure
% plot(nlml2lane');


