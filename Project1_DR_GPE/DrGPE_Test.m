% DrGPE_Test
% Prototype I


close all

for i=1:length(num_train_ref)
    
    % Generate dataset
    num_train=num_train_ref(i);
    Z=Z_Rec(1:num_train,1:dim_new,i);
    model=model_Rec(i);
    
    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);

    
    for j=1:dim_new
            hyp=hyp_Rec(i,j);
            [m(:,j) s(:,j)] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, X, Z(:,j), X_star);
         
    %         %MCMC
    %         [m(:,i) s(:,i)] = gp(hyp, @infMCMC, meanfunc, covfunc, likfunc, X, Z(:,i), X_star);     
    end
    
    for k=1:dim_new

        Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);
    
        SquErr=(Y_starorig-Y_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

        SSE(:,k)=sum(SquErr');

    end
        
    
    figure
    boxplot(SSE)
    title('Square sum error for different subspace dimension -kPCA')
    set(gca,'yscale','log');

    mean_SE(i,:)=mean(SSE);
    SSE=[];
    m=[];
    s=[];
   
end

    figure
    semilogy(mean_SE');
    title('Mean of Square sum error for different subspace dimension')
%     legend('DiffusionMaps','kPCA')
    
