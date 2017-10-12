% Test_GPE_kPCA

clear
% close all

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
num_test=300;

%--------------------------------------------------------------------------
dim_new=10;           % temp value.

options.ker='gaussian';
options.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
options.new_dim=dim_new;
options.FullRec=0;

preoptions.type='Exp';  %'LSE' OR 'Dw'
% preoptions.para=2;
% preoptions.neighbor=10;


% %Kpca options
% kpca.options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
% %    options = struct('ker','rbf','arg',50,'new_dim',new_dim);                  % other options
% %    options = struct('ker','poly','arg',[2,0],'new_dim',new_dim);              % other options
% %    options = struct('ker','sigmoid','arg',[0.000001,0],'new_dim',new_dim);    % other options
% % isomap.options = struct('dim_new',new_dim,'neighbor',10,'d2p_method','Dw', 'd2p_Dwpara',10,'d2p_points',10);      
% 
% %--------------------------------------------------------------------------
% % Isomap options
% isomap.options.dim_new=new_dim;                      % New dimension
% isomap.options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
% isomap.options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
% isomap.options.metric='euclidean';             % Method of measurement. Metric
% 
% % Isomap PreImage options
% isomap.Reoptions.ReCoverNeighborType='k';% Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
% isomap.Reoptions.ReCoverNeighborPara=10;       % Parameter of neighbor of new point
% isomap.Reoptions.Recoverd2pType='Dw';          % Type of distance to coordinate method. Distance weight/Least square estimate
% isomap.Reoptions.Recoverd2pPara=1;             % Parameter of distance to coordinate recover method
% 

%% Main program
% -----------------Defining the input ------------------------------------
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);


    %% kPCA
    Distance =pdist2(Y,Y,'euclidean');
    options.arg=sum(Distance(:).^2)/(num_train^2);   
    options.arg=sqrt(options.arg/2);

    [Z,model] = Kpca(Y,options);
    
    
    %% GP
%     %Structure 1
%     meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
%     covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
%     likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

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
    covfunc = {@covSum,{@covSEard,@covNoise}}; 
    hyp.cov = [zeros(Dim_X+1,1);0];
%     hyp.cov = hyp.cov-inf;
    
    likfunc = @likGauss; 
    sn = 0.1;
    hyp.lik = log(sn);

    meanfunc=[];
    hyp.mean=[];
  
   
    % TRAIN & TEST EXACT
    for i=1:dim_new
        
        % MLE method
        hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i));
        exp(hyp.lik)
        nlml2 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i))        
        [m(:,i) s(:,i)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);
        
        %wipe memory
%         hyp.cov = [zeros(Dim_X+1,1);-Inf];
        hyp.cov = [zeros(Dim_X+1,1);0];
        hyp.lik = log(sn);
        hyp.mean=[];
        
%         %MCMC
%         cov = {@covSum,{@covSEiso,@covNoise}}; hyp.cov = [0; 2; -Inf];       % ell,sf,sn
%         %cov = @covSEiso; hyp.cov = [0; 2]; % FAIL
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
%          par.Nsample = 400;
%         % par.sampler = 'ess'; % par.Nais = 5; par.Nsample = 200; par.Nskip = 100;
%         tic
%         [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,X,Z(:,i),par);
%         % [posts,nlZs,dnlZs] = infMCMC(hyp,mn,cov,lik,x,y);
%         toc
%         [m(:,i),s(:,i),fmus,fs2s,junk,posts] = gp(hyp,@infMCMC,mn,cov,lik,X,posts,X_star);
        
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
        
        
    end
    
    for i=1:dim_new
    
        Y_star = Kpca_PreImage(m(:,1:i),model,preoptions);   
        SquErr=(Y_starorig-Y_star).^2;       
        SSE(:,i)=sum(SquErr');

    end

    


    figure
    boxplot(SquErr')
    title('Square error for all grid point for each test point -kPCA')
    set(gca,'yscale','log');

    figure
    boxplot(SSE)
    title('Square sum error for different subspace dimension -kPCA')
    set(gca,'yscale','log');

    mean_SE=mean(SSE);


    figure
    semilogy(mean_SE),'k-';
    title('Mean of Square sum error for different subspace dimension -kPCA')
    legend('kPCA')
    
    
    
    
%  %Each details ploting   
%     
%     
%     Y_star = Kpca_PreImage(m,model,preoptions);
% 
% %     Y_star=real(Y_star); 
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
