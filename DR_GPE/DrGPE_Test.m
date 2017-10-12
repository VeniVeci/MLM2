% DrGPE_Test
% Dimension reduction Gausian process emulation_Test
% Test trained GPE of DrGPE_Train
%
% Instructions:
% 
% Modifications:
% WeiX, 5-1-2016, Create

%% Initialize
% clear
% close all 

% load('Data5kPCAinfExact.mat')
preoptions.type='Exp';  %'Exp' method is more stable and non-parametric
num_test=300;


mean_SE=[];

for i=1:length(num_train_ref)
      
    dim_extract=1;
    % Collect dataset
    num_train=num_train_ref(i);
    Z=Z_Rec(1:num_train,:,i);
    model=model_Rec(i);
    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);

    %Solving Pre-image with increasing dimensions
    
    [num_Z,dim_Z]=size(Z);    
    %Clean memory
    m=zeros(num_test,dim_Z);
    s=zeros(num_test,dim_Z);
      
    for j=1:dim_Z
            hyp=hyp_Rec(i,j);
            [m(:,j) s(:,j)] = gp(hyp, InfMethod, meanfunc, covfunc, likfunc, X, Z(:,j), X_star);  
    end
    
    
    %Clean memory
    SSE=zeros(num_test,dim_new);
    
    for k=1:dim_new    
        switch DrMethod    
            case 'kPCA'
                Y_star = Kpca_PreImage(m(:,1:k),model,preoptions);
            case 'DiffusionMaps'   
                switch options.Ztype
                    case 0
                        Y_star = DiffusionMaps_PreImage(m(:,1:k),model,preoptions);
                    case 1
                        Y_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
                end            

            case 'Isomaps'
                Y_star= Isomaps_PreImage(m(:,1:k),model,preoptions);
            otherwise 
                error('No such DR method')
        end

        SquErr=(Y_starorig-Y_star).^2;        
        MSS_orig=mean(Y_starorig.^2,2);

        SSE(:,k)=sum(SquErr,2);
        RSSE(:,k)=mean(SquErr,2)./MSS_orig;

    %         SE_D(:,i)=SquErr_D(:);
    %         SE_K(:,i)=SquErr_K(:);

    end
      
    % %Ratio of variance.mean of coefficients
    % VarRelaError=abs(s./m);
    % sumVarRelaError(i,:)=sum(VarRelaError);     
    % figure(10+i)
    % boxplot(VarRelaError)
    % title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
    % set(gca,'yscale','log');
    % xlabel('Index of principal direction')
    % ylabel('Predicted variance ratio')
           
    figure(i)
    boxplot(SSE)
    title(sprintf('Boxplot of Square sum error for different subspace dimension -%s -%d Training points', DrMethod, num_train_ref(i)))
    set(gca,'yscale','log');
    xlabel('Dimension of manifold')
    ylabel('Square sum error')

     mean_SE(i,:)=mean(SSE);
%     SSE=[];
%     m=[];
%     s=[];

end

%% Result presentation
markers = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
[num_temp,dim_temp]=size(mean_SE);
figure
hold on
for i=1:num_temp
    semilogy(mean_SE(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('Mean of Square sum error for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimension of manifold')
ylabel('Square sum error mean value')





    
