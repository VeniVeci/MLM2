% DrGPE_HQ_Test_NS
% A HQ terminal for dimension reduction GPE _Special design for Navier
% stock equation for the special format of dataset( e.g. using U,V,P as Y
% which are multi-output.
%% Initialize
% clear
% close all 

% load('Data5kPCAinfExact.mat')
preoptions.type='Exp';  %'Exp' method is more stable and non-parametric
preoptions.neighbor=10;

num_test=300;
Test_split=201;


mean_SE=[];

for i=1:length(num_train_ref)
    
   
    % Special design for special structure of Diffusion maps.
    if strcmp(DrMethod,'DiffusionMaps' )
        dim_extract=dim_new+1;
    else 
        dim_extract=dim_new;
    end
    
    % Collect dataset
    num_train=num_train_ref(i);
    Z=Z_Rec(1:num_train,1:dim_extract,i);
    model=model_Rec(i);
    
    Y=Y_Data(1:num_train,:);
    X=X_Data(1:num_train,:);
    
    X_star=X_Data(Test_split:Test_split+num_test-1,:);
    Y_starorig=Y_Data(Test_split:Test_split+num_test-1,:);
    
    % Solving Pre-image with increasing dimensions
    
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
                Y_star = DiffusionMaps_PreImage(m(:,1:k+1),model,preoptions);
            case 'Isomaps'

            otherwise 
                error('No such DR method')
        end

        SquErr=(Y_starorig-Y_star).^2;
        
%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);

        SSE(:,k)=sum(SquErr');
        
    end
    
     VarRelaError=abs(s./m);
     sumVarRelaError(i,:)=sum(VarRelaError);
%     figure(10+i)
%     boxplot(VarRelaError)
%     title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
%     set(gca,'yscale','log');
%     xlabel('Index of principal direction')
%     ylabel('Predicted variance ratio')
    
           
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


index=[17:20];
figure
hold on
plot(Y_starorig(index,:)','b-')
plot(Y_star(index,:)','k--')
hold off


    
