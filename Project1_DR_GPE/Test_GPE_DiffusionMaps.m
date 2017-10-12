% Test_GPE_Diffusion maps

clear
close all

%Assign parameter

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
num_train=120;                 
num_test=100;

dim_new=10;

options.metric ='euclidean';
options.kernel ='gaussian'; 
% options.kpara = 10000;             
options.kAuto=1;

options.dim_new = dim_new;              
options.t = 1;                     
options.FullRec = 0;      

% Diffusion PreImage options
preoptions.type='Dw';  %'LSE' OR 'Dw'
preoptions.para=2;
preoptions.neighbor=10;
% preoptions.neighbor=10;

%Assign memory for recoder
SSE=[];

%% PCA-GPR test

    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %% Diffusion maps
    
    [Z,model] = DiffusionMaps(Y,options);
    
    
    
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
    for i=1:dim_new+1
%         hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));
%         exp(hyp.lik)
%         nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i))
%         [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
     
        %MCMC
        cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
%         cov = @covSEiso; hyp.cov = [0; 2]; % FAIL
        lik =  {'likGauss'};                 % likLogistic or likErf
        sn = 0.1;
        hyp.lik = log(sn);
        mn  = {'meanZero'};      hyp.mean = [];
        % 4) prediction using MCMC sampling
        % set MCMC parameters, see some more details in inf/infMCMC.m
        % We have two samplers implemented, namely
        %  hmc - Hybrid Monte Carlo, and
        %  ess - Elliptical Slice Sampling.
        par.sampler = 'hmc'; 
        % par.sampler = 'ess'; 
         par.Nsample = 400;
        % par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
        tic
        [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,X,Z(:,i),par);
        % [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y);
        toc
        [m(:,i),s(:,i),fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,X,posts,X_star);
        
        % FORCE Special Correction Treatment
        m(:,1)=Z(1,1);
        
%         %Plot the predictive mean and variance
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
%         legend('mean','+2xSigma','-2xSigma')



%         figure(i)
%         mup=m(:,i)+2*sqrt(s(:,i));
%         mdown=m(:,i)-2*sqrt(s(:,i));
%         figure(i)
%         plot(m(:,i),'-');
%         hold on
%         plot(mup,'--');
%         plot(mdown,'--');
%         
%         xScale=[1:num_test]';
%         f = [m(:,i)+2*sqrt(s(:,i)); flipdim(m(:,i)-2*sqrt(s(:,i)),1)];
%         fill([xScale; flipdim(xScale,1)], f, [7 7 7]/8)
% 
%         hold on; plot(xScale, m(:,i)); 
%         
%         hold off
%         
%         legend('mean','+2xSigma','-2xSigma')
        
        
    end
    
    for i=1:dim_new
        Y_star = DiffusionMaps_PreImage(m(:,1:i+1),model,preoptions);
        SquErr=(Y_starorig-Y_star).^2;       
        SSE(:,i)=sum(SquErr');

    end

    
    figure
    boxplot(SquErr')
    title('Square error for all grid point for each test point -DiffusionMaps')
    set(gca,'yscale','log');

    figure
    boxplot(SSE)
    title('Square sum error for different subspace dimension -DiffusionMaps')
    set(gca,'yscale','log');

    mean_SE=mean(SSE);
    figure
    semilogy(mean_SE,'d');
    title('Mean of Square sum error for different subspace dimension -DiffusionMaps')
    legend('DiffusionMaps')
    
    
   
% Y_star = DiffusionMaps_PreImage(m,model,preoptions);

%     Y_star=real(Y_star); 
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
     
    figure
    SquErr=(Y_starorig-Y_star).^2;
    boxplot(SquErr');
    
  
    figure
    ReSquErr=SquErr./abs(Y_starorig);
    boxplot(ReSquErr');
%     axis([0 1])
    ylim([0 0.3])

    
%     index=[1,11,21];
%     figure
%     plot(Y_star(index,:)','--b')
%     hold on 
%     plot(Y_starorig(index,:)','-k')
%     hold off
% 
%     index=1;
%     figure
%     field=reshape(Y_star(index,:),[100,100]);
%     surf(field)
%     title('Prediced')
%     figure
%     field=reshape(Y_starorig(index,:),[100,100]);
%     surf(field)
%     title('Original')   
%     
    
    
    
    
    
    
    
    
%     [Y_star]=Lpca_PreImage(m,model);
%     
% 
% %     for i=1:dim_new
% %         figure(i)
% %         mup=m(:,i)+2*sqrt(s(:,i));
% %         mdown=m(:,i)-2*sqrt(s(:,i));
% %         figure(i)
% %         plot(m(:,i),'-');
% %         hold on
% %         plot(mup,'--');
% %         plot(mdown,'--');
% %         hold off
% % 
% %     end
%      
%     figure(dim_new+1)
%     SquErr=(Y_starorig-Y_star).^2;
%     boxplot(SquErr');
%     
%   
%     figure(dim_new+2)
%     ReSquErr=SquErr./abs(Y_starorig);
%     boxplot(ReSquErr');
% %     axis([0 1])
%     ylim([0 0.1])
%     
%     index=[10:10:50];
%     figure
%     plot(Y_star(index,:)','--b')
%     hold on 
%     plot(Y_starorig(index,:)','-k')
%     hold off
%     
%     
%     
% %     index=1;
% %     figure
% %     field=reshape(Y_star(index,:),[100,100]);
% %     surf(field)
% %     title('Prediced')
% %     figure
% %     field=reshape(Y_starorig(index,:),[100,100]);
% %     surf(field)
% %     title('Original')
%     
%     
%     %% ADDITION
%     Z_starorig=Y_starorig*model.Key';
%     
