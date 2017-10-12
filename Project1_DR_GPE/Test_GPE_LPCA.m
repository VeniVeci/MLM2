% Test_GPE_LPCA

clear
close all

%Assign parameter

%% Dataset Parameter
index_dataset=4;
            % 1: 5050PloeFlow_Filtered
            % 2: 5050PloeFlow
            % 3: 100100PloeFlow
            % 4: original_TinMeltingFront100by100
            % 5: originalCon
            % 6: originalSuperConductingWire50by50
            % 7: originalSuperConductingWire100by100
            % 8: CSTR
num_train=40;                 
num_test=50;

dim_new=10;


%Assign memory for recoder
SSE=[];

%% PCA-GPR test

    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %% PCA
    
    [Z,model] = Lpca(Y,dim_new);
    
    
    %% GP
    % %Structure 1
%     meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
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

    % TRAIN & TEST MCMC
    % hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, x, y);
    % exp(hyp.lik)
    % nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, x, y)
    % [m s] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, x, y, z);


    % TRAIN & TEST EXACT
    for i=1:dim_new
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));
        exp(hyp.lik)
        nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i))
        [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
     
%         %MCMC
%         hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i));
%         exp(hyp.lik)
%         nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i))
%         [m(:,i) s(:,i)] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
        
        figure(i)
        mup=m(:,i)+2*sqrt(s(:,i));
        mdown=m(:,i)-2*sqrt(s(:,i));
        figure(i)
        plot(m(:,i),'-');
        hold on
        plot(mup,'--');
        plot(mdown,'--');
        
        xScale=[1:num_test]';
        f = [m(:,i)+2*sqrt(s(:,i)); flipdim(m(:,i)-2*sqrt(s(:,i)),1)];
        fill([xScale; flipdim(xScale,1)], f, [7 7 7]/8)

        hold on; plot(xScale, m(:,i)); 
        
        hold off
        
        legend('mean','+2xSigma','-2xSigma')
        
        
    end
    
    [Y_star]=Lpca_PreImage(m,model);
    

%     for i=1:dim_new
%         figure(i)
%         mup=m(:,i)+2*sqrt(s(:,i));
%         mdown=m(:,i)-2*sqrt(s(:,i));
%         figure(i)
%         plot(m(:,i),'-');
%         hold on
%         plot(mup,'--');
%         plot(mdown,'--');
%         hold off
% 
%     end
     
    figure(dim_new+1)
    SquErr=(Y_starorig-Y_star).^2;
    boxplot(SquErr');
    
  
    figure(dim_new+2)
    ReSquErr=SquErr./abs(Y_starorig);
    boxplot(ReSquErr');
%     axis([0 1])
    ylim([0 0.1])
    
    index=[10:10:50];
    figure
    plot(Y_star(index,:)','--b')
    hold on 
    plot(Y_starorig(index,:)','-k')
    hold off
    
    
    
%     index=1;
%     figure
%     field=reshape(Y_star(index,:),[100,100]);
%     surf(field)
%     title('Prediced')
%     figure
%     field=reshape(Y_starorig(index,:),[100,100]);
%     surf(field)
%     title('Original')
    
    
    %% ADDITION
    Z_starorig=Y_starorig*model.Key';
    
%   %%   

%     exp(hyp.lik)
%     nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Y)
%     [m s] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Y, Z);
%     
%     
%     
%     [X_star]=Lpca_PreImage(Z,model);
%     
%    
%     [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);
%     options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
%     [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
%      
% 
%     SquErr=(Y_starorig-Y_star_kpca).^2;
%     means=mean(Y_starorig,2);
%     SSErr=sum(SquErr,2);
%     RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
%     RateErr=abs(Y_starorig-Y_star_kpca)./Y_starorig;
%     RateErr=mean(RateErr,2);
% 
%     RecSSErr_kpca(:,i)   =SSErr;
%     RecRateSsErr_kpca(:,i)=RatSsErr;
%     RecRateErr_kpca(:,i)=RateErr;
%     RecTime_kpca(:,i)=t_kpca;
% 
%     %SVD part
%     SquErr=(Y_starorig-Y_star_svd).^2;
%     means=mean(Y_starorig,2);
%     SSErr=sum(SquErr,2);
%     RatSsErr=sqrt(SSErr)./(means*Dim_Y); %=sqrt(SSErr)./sum((Y_starorig,2))
%     RateErr=abs(Y_starorig-Y_star_svd)./Y_starorig;
%     RateErr=mean(RateErr,2);
% 
%     RecSSErr_svd(:,i)   =SSErr;
%     RecRateSsErr_svd(:,i)=RatSsErr;
%     RecRateErr_svd(:,i)=RateErr;
%     RecTime_svd(:,i)=t_svd;
% 
% 
%     
% 
%     
%     
%     
% 
% figure(1)
% boxplot(RecRateSsErr_kpca(:,1:1:i-1),{1:1:i-1});
% title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));
% figure(2)
% boxplot(RecRateSsErr_svd(:,1:1:i-1),{1:1:i-1});
% title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));
% 
% RecRateSsErr_comb=[];
% for j=1:i-1
%     RecRateSsErr_comb=[RecRateSsErr_comb,RecRateSsErr_kpca(:,j),RecRateSsErr_svd(:,j)]; 
% end
% figure(3)
% boxplot(RecRateSsErr_comb);
% title(sprintf('Square Sum Error Rate of Each pixel of KLPCA-GPR'));
% 
% 
% 
% 
% figure(4)
% boxplot(RecRateErr_kpca(:,1:1:i-1),{1:1:i-1});
% title(sprintf('Square Sum Error Rate of Each pixel of KPCA-GPR'));
% figure(5)
% boxplot(RecRateErr_svd(:,1:1:i-1),{1:1:i-1});
% title(sprintf('Square Sum Error Rate of Each pixel of LPCA-GPR'));
% 
% RecRateErr_comb=[];
% for j=1:i-1
%     RecRateErr_comb=[RecRateErr_comb,RecRateErr_kpca(:,j),RecRateErr_svd(:,j)]; 
% end
% figure(6)
% boxplot(RecRateErr_comb);
% title(sprintf('Square Sum Error Rate of Each pixel of KLPCA-GPR'));
% 

