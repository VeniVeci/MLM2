% Test_GPE_ISsomaps

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
num_train=60;                 
num_test=100;

%--------------------------------------------------------------------------
dim_new=20;           % temp value.

options.dim_new=dim_new;                      % New dimension
options.neighborType='k';               % Type of neighbor.Choice:1)'k';Choice:2)'epsilon'
options.neighborPara=10;                % parameter for choosing the neighbor. number of neighbor for "k" type and radius for 'epsilon'
options.metric='euclidean';             % Method of measurement. Metric

% Isomap PreImage options
preoptions.neighborType='k';    % Type of neighbor of new point. Choice:1)'k';Choice:2)'epsilon'
preoptions.neighborPara=10;     % Parameter of neighbor of new point
preoptions.type='Dw';           % Type of distance to coordinate method. Distance weight/Least square estimate
preoptions.para=2;              % Parameter of distance to coordinate recover method


%% Main program
% -----------------Defining the input ------------------------------------
[X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
[np_train,Dim_X]=size(X);
[np_test, Dim_Y]=size(Y_starorig);


    %% ISOMAPS    
%     [Z,model] = Kpca(Y,options);
    [Z,model] = Isomaps(Y,options);
    
    %% GP
    % %Structure 1
%     meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1;1;1];
%     covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]); % could be working EXTREMELY well
%     likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

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
        hold off     
        legend('mean','+2xSigma','-2xSigma')
end
        
        
for i=1:dim_new
    
        Y_star = Isomaps_PreImage(m(:,1:i),model,preoptions);        
        SquErr=(Y_starorig-Y_star).^2;
        SSE(:,i)=sum(SquErr');

        
end
        
    
%     Y_star = Isomaps_PreImage(m,model,preoptions);

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
    
  
    figure
    ReSquErr=SquErr./abs(Y_starorig);
    boxplot(ReSquErr');
%     axis([0 1])
    ylim([0 0.3])

    figure
    boxplot(SSE);    
    
    
    mean_SE=mean(SSE);
    figure
    semilogy(mean_SE,'d-');
    title('Mean of Square sum error for different subspace dimension')

    
    
    index=[10:10:50];
    figure
    plot(Y_star(index,:)','--b')
    hold on 
    plot(Y_starorig(index,:)','-k')
    hold off
    


% different way of sorting the index    
%     [~,Index] = sort(sum(Y,2));       % not good
%     [~,Index] = sort(Y(:,1));
%     [~,Index] = sort(Y(:,end));
    [~,Index] = sort(Z(:,1));           % Just OK
%     [~,Index] = sort(max(Y,[],2));    % not good
    
    SufrY=Y(Index,:);
    
    figure
    mesh(SufrY)
    
    figure
    plot(SufrY','-k')
    

    
    
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
    
    

