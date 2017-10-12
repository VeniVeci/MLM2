% Test_GPE_kDI
% Test_GPE_kPCA & Diffusion maps & Isomap

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
num_train=50;                 
num_test=100;

%--------------------------------------------------------------------------
dim_new=20;           % temp value.

% kernel PCA options----------------------
Koptions.ker='gaussian';
Koptions.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
% Koptions.arg=500;    % 500 For index_dataset=8
Koptions.new_dim=dim_new;
Koptions.FullRec=0;


% Diffusion Maps options--------------------------------------
Doptions.metric ='euclidean';
Doptions.kernel ='gaussian'; 
% Doptions.kpara = 10000;             
Doptions.kAuto=1;
Doptions.dim_new = dim_new;              
Doptions.t = 1;                     
Doptions.FullRec = 0;      

% Isomaps options------------------------------------------------------
Ioptions.dim_new=dim_new;                      % New dimension
Ioptions.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
Ioptions.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
Ioptions.metric='euclidean';             % Method of measurement. Metric

% PreImage options------------------------------------------------------
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
%     Distance =pdist2(Y,Y,'euclidean');
%     Koptions.arg=sum(Distance(:))/(num_train^2);   
    
    [Z_K,model_K] = Kpca2(Y,Koptions);    
    [Z_D,model_D] = DiffusionMaps(Y,Doptions);
    [Z_I,model_I] = Isomaps(Y,Ioptions);
    
    
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
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z_K(:,i));
    exp(hyp.lik)
    nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_K(:,i))
    [m_K(:,i) s_K(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_K(:,i), X_star);

%     %MCMC
%     hyp = minimize(hyp, @gp, -100, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i));
%     exp(hyp.lik)
%     nlml2 = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i))
%     [m(:,i) s(:,i)] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
% 
% 
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
    
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z_I(:,i));
    exp(hyp.lik)
    nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_I(:,i))
    [m_I(:,i) s_I(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_I(:,i), X_star);


end
        
for i=1:dim_new+1
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z_D(:,i));
        exp(hyp.lik)
        nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_D(:,i))
        [m_D(:,i) s_D(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z_D(:,i), X_star);
            
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
        Y_I_star = Isomaps_PreImage(m_I(:,1:i),model_I,preoptions);
        
        SquErr_D=(Y_starorig-Y_D_star).^2;
        SquErr_K=(Y_starorig-Y_K_star).^2;
        SquErr_I=(Y_starorig-Y_I_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);
        
        SSE_D(:,i)=sum(SquErr_D');
        SSE_K(:,i)=sum(SquErr_K');
        SSE_I(:,i)=sum(SquErr_I');
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
boxplot(SquErr_I')
title('Square error for all grid point for each test point -Isomap')
set(gca,'yscale','log');



figure(4)
boxplot(SSE_D)
title('Square sum error for different subspace dimension -DiffusionMaps')
set(gca,'yscale','log');

figure(5)
boxplot(SSE_K)
title('Square sum error for different subspace dimension -kPCA')
set(gca,'yscale','log');

figure(6)
boxplot(SSE_I)
title('Square sum error for different subspace dimension -Isomaps')
set(gca,'yscale','log');




mean_SE_D=mean(SSE_D);
mean_SE_K=mean(SSE_K);
mean_SE_I=mean(SSE_I);


figure(7)
semilogy(mean_SE_D,'d');
hold on 
semilogy(mean_SE_K),'t';
semilogy(mean_SE_I,'o');
hold off 
title('Mean of Square sum error for different subspace dimension -DiffusionMaps')
legend('DiffusionMaps','kPCA','Isomaps')





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
    
    