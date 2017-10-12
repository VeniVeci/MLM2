% Test_GPE_kD
% Test_GPE_kPCA & Diffusion maps

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
num_train=100;                 
num_test=100;

%--------------------------------------------------------------------------
dim_new=10;           % temp value.

% kernel PCA options----------------------
Koptions.ker='gaussian';
Koptions.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
% Koptions.arg=500;    % 500 For index_dataset=8
Koptions.new_dim=dim_new;
Koptions.FullRec=0;


% Diffusion Maps options----------------------
Doptions.metric ='euclidean';
Doptions.kernel ='gaussian'; 
% Doptions.kpara = 10000;             
Doptions.kAuto=1;

Doptions.dim_new = dim_new;              
Doptions.t = 1;                     
Doptions.FullRec = 0;      

% PreImage options----------------------
preoptions.type='Dw';  %'LSE' OR 'Dw'
preoptions.para=2;
preoptions.neighbor=10;
% preoptions.neighbor=10;


%%  test

    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %% Diffusion maps
    
     Distance =pdist2(Y,Y,'euclidean');
     Koptions.arg=sum(Distance(:).^2)/(num_train^2);   
     Koptions.arg=sqrt(Koptions.arg/2);
    
    [Z_K,model_K] = Kpca(Y,Koptions);    
    [Z_D,model_D] = DiffusionMaps(Y,Doptions);
    
    
    
%% GP
% % Structure 1
% meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
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

for i=1:dim_new
    %MLE
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z_K(:,i));
    exp(hyp.lik)
    nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_K(:,i))
    [m_K(:,i) s_K(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_K(:,i), X_star);


%     %MCMC
%     cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
%     lik =  {'likGauss'};                 % likLogistic or likErf
%     sn = 0.1;
%     hyp.lik = log(sn);
%     mn  = {'meanZero'};      hyp.mean = [];
%     % 4) prediction using MCMC sampling
%     % set MCMC parameters, see some more details in inf/infMCMC.m
%     % We have two samplers implemented, namely
%     %  hmc - Hybrid Monte Carlo, and
%     %  ess - Elliptical Slice Sampling.
%     par.sampler = 'hmc'; 
%     % par.sampler = 'ess'; 
%     % par.Nsample = 20;
%     % par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
%     tic
%     [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,X,Z_K(:,i),par);
%     % [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y);
%     toc
%     [m_K(:,i),s_K(:,i),fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,X,posts,X_star);


%     figure(i)
%     mup=m(:,i)+2*sqrt(s(:,i));
%     mdown=m(:,i)-2*sqrt(s(:,i));
%     figure(i)
%     plot(m(:,i),'-');
%     hold on
%     plot(mup,'--');
%     plot(mdown,'--');
%     hold off
% 
%     legend('mean','+2xSigma','-2xSigma')

end
        
for i=1:dim_new+1
        %MLE
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z_D(:,i));
        exp(hyp.lik)
        nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_D(:,i))
        [m_D(:,i) s_D(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_D(:,i), X_star);
     
%         %MCMC
%         cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
%         lik =  {'likGauss'};                 % likLogistic or likErf
%         sn = 0.1;
%         hyp.lik = log(sn);
%         mn  = {'meanZero'};      hyp.mean = [];
%         % 4) prediction using MCMC sampling
%         % set MCMC parameters, see some more details in inf/infMCMC.m
%         % We have two samplers implemented, namely
%         %  hmc - Hybrid Monte Carlo, and
%         %  ess - Elliptical Slice Sampling.
%         par.sampler = 'hmc'; 
%         % par.sampler = 'ess'; 
%         % par.Nsample = 20;
%         % par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
%         tic
%         [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,X,Z_D(:,i),par);
%         % [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y);
%         toc
%         [m_D(:,i),s_D(:,i),fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,X,posts,X_star);
        

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
        
end
 


% SquErr_D=(Y_starorig-Y_D_star).^2;
% SquErr_K=(Y_starorig-Y_K_star).^2;
% 
% ReSquErr_D=SquErr_D./abs(Y_starorig);
% ReSquErr_K=SquErr_K./abs(Y_starorig);
% 
% figure(1)
% boxplot([SquErr_D(:),SquErr_K(:)])
% 
% figure(2)
% boxplot([ReSquErr_D(:),ReSquErr_K(:)])
% 
% SSE_D=sum(SquErr_D(:))
% SSE_K=sum(SquErr_K(:))


%% Plot the trend with number of PC
for i=1:dim_new
    
        Y_D_star = DiffusionMaps_PreImage(m_D(:,1:i+1),model_D,preoptions);
        Y_K_star = Kpca_PreImage(m_K(:,1:i),model_K,preoptions);
    
        SquErr_D=(Y_starorig-Y_D_star).^2;
        SquErr_K=(Y_starorig-Y_K_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);
        
        SSE_D(:,i)=sum(SquErr_D');
        SSE_K(:,i)=sum(SquErr_K');
end



figure(1)
boxplot(SquErr_D')
title('Square error for all grid point for each test point -DiffusionMaps')
set(gca,'yscale','log');

figure(2)
boxplot(SquErr_K')
title('Square error for all grid point for each test point -kPCA')
set(gca,'yscale','log');

figure(3)
boxplot(SSE_D)
title('Square sum error for different subspace dimension -DiffusionMaps')
set(gca,'yscale','log');

figure(4)
boxplot(SSE_K)
title('Square sum error for different subspace dimension -kPCA')
set(gca,'yscale','log');


mean_SE_D=mean(SSE_D);
mean_SE_K=mean(SSE_K);

figure(5)
semilogy(mean_SE_D,'d-');
hold on 
semilogy(mean_SE_K),'k-';
hold off 
title('Mean of Square sum error for different subspace dimension -DiffusionMaps')
legend('DiffusionMaps','kPCA')





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
     
%     figure
%     SquErr=(Y_starorig-Y_star).^2;
%     boxplot(SquErr');
%     
%   
%     figure
%     ReSquErr=SquErr./abs(Y_starorig);
%     boxplot(ReSquErr');
% %     axis([0 1])
%     ylim([0 0.3])
% 
%     
%     index=[1,11,21];
%     figure
%     plot(Y_star(index,:)','--b')
%     hold on 
%     plot(Y_starorig(index,:)','-k')
%     hold off
% 
    index=1;
    figure
    field=reshape(Y_D_star(index,:),[100,100]);
    surf(field)
    title('Prediced')
    figure
    field=reshape(Y_starorig(index,:),[100,100]);
    surf(field)
    title('Original')   
    
    