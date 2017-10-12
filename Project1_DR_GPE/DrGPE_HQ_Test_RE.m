% DrGPE_HQ_Test_RE
% A HQ terminal for dimension reduction GPE show relative error

%% Initialize
% clear
% close all 

% load('Data5kPCAinfExact.mat')
preoptions.type='Exp';  %'Exp' method is more stable and non-parametric
num_test=300;


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
    [X,Y,X_star,Y_starorig]=Dataset_Get(num_train,num_test,index_dataset);

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
        MSS_orig=mean(Y_starorig.^2,2);
        
        SSE(:,k)=sum(SquErr,2);
        RSSE(:,k)=mean(SquErr,2)./MSS_orig;

%         SE_D(:,i)=SquErr_D(:);
%         SE_K(:,i)=SquErr_K(:);
    
    end
    
%      VarRelaError=abs(s./m);
%      sumVarRelaError(i,:)=sum(VarRelaError);
     
%     figure(10+i)
%     boxplot(VarRelaError)
%     title(sprintf('Boxplot of variance ratio for each principal direction  -%s', DrMethod))
%     set(gca,'yscale','log');
%     xlabel('Index of principal direction')
%     ylabel('Predicted variance ratio')
    
           
    figure(i)
    boxplot(RSSE)
    title(sprintf('Boxplot of Square sum error for different subspace dimension -%s -%d Training points', DrMethod, num_train_ref(i)))
    set(gca,'yscale','log');
    xlabel('Dimension of manifold')
    ylabel('Square sum error')

     MRSSE(i,:)=mean(RSSE);
%     SSE=[];
%     m=[];
%     s=[];

end


markers = ['+','o','*','.','x','s','d','^','v','>','<','p','h'];
[num_temp,dim_temp]=size(MRSSE);
figure
hold on
for i=1:num_temp
    semilogy(MRSSE(i,:)',markers(i));
    entries{i}=sprintf('%d Training points',num_train_ref(i));
end
hold off
set(gca,'yscale','log');
title(sprintf('Relative error (MRSSE=Mean SSE/Morig^2) for different subspace dimension -%s', DrMethod))
legend(entries)
xlabel('Dimensions')
ylabel('Relative error')



%Detail plot

    index=26;
    figure
    field=reshape(Y_star(index,:),[50,200]);
    surf(field)
    contourf(field,8)
    colorbar
    title(sprintf('%s GPE Predicted field -%d Training points x1=%3.2f x2= %3.3f ', DrMethod,num_train,X_star(index,1),X_star(index,2)))
    
    figure
    field=reshape(Y_starorig(index,:),[50,200]);
    surf(field)
    contourf(field,8)
    colorbar
    title(sprintf(' Real field x1=%3.2f x2= %3.3f ',X_star(index,1),X_star(index,2)))

    