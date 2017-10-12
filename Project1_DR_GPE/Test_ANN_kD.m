% Test_ANN_kD
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
num_train=40;                 
num_test=50;

%--------------------------------------------------------------------------
dim_new=20;           % temp value.


Koptions.ker='gaussian';
Koptions.arg=10;   %10 wont work well,WHY? model.K is too normal distributed which leads to slow eigenvalue decay!!! So options.arg tends to be large.
% Koptions.arg=500;    % 500 For index_dataset=8
Koptions.new_dim=dim_new;
Koptions.FullRec=0;

Doptions.metric ='euclidean';
Doptions.kernel ='gaussian'; 
% Doptions.kpara = 10000;             
Doptions.kAuto=1;

Doptions.dim_new = dim_new;              
Doptions.t = 1;                     
Doptions.FullRec = 0;      

% Diffusion PreImage options
preoptions.type='Dw';  %'LSE' OR 'Dw'
preoptions.para=2;
preoptions.neighbor=10;
% preoptions.neighbor=10;



SSE=[];

%%  test

    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);
    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %% Diffusion maps    
    [Z_D,model_D] = DiffusionMaps(Y,Doptions);
    [Z_K,model_K] = Kpca(Y,Koptions);
    
%     Z_K=real(Z_K);   
    
    
%% GP
kfold=1;
for i=1:dim_new
    
    m_D= ANN_prediction_crossvalidation(X',Z_D',X_star',kfold);   
    m_K= ANN_prediction_crossvalidation(X',Z_K',X_star',kfold);

end
m_D=m_D';
m_K=m_K';


%% Plot the trend with number of PC
for i=1:dim_new
    
        Y_D_star = DiffusionMaps_PreImage(m_D(:,1:i+1),model_D,preoptions);
        Y_K_star = Kpca_PreImage(m_K(:,1:i),model_K,preoptions);
        
        
        Y_K_star=real(Y_K_star); % need further investigate
    
        SquErr_D=(Y_starorig-Y_D_star).^2;
        SquErr_K=(Y_starorig-Y_K_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);
        
        SSE_D(:,i)=sum(SquErr_D');
        SSE_K(:,i)=sum(SquErr_K');
        
end

% SSE_K=real(SSE_K);


figure(1)
boxplot(SquErr_D')
title('Square error for all grid point -DiffusionMaps')

figure(2)
boxplot(SquErr_K')
title('Square error for all grid point -kPCA')

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
semilogy(mean_SE_D,'d');
hold on 
semilogy(mean_SE_K,'*');
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
    
    