% DrGPE_HQ_Test_Coefficients
% A HQ terminal for dimension reduction GPE

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
    
%     for j=1:dim_Z
%         mu=m(:,j);
%         sigma=s(:,j);
%               
%         [mu,index]=sort(mu);
%         sigma=sigma(index);
%         
%         figure(j)
%         errorbar(mu,sigma,'rx')
%     end
    
        dim=6;
        mu=m(:,dim);
        sigma=s(:,dim);
              
        [mu,index]=sort(mu);
        sigma=sigma(index);
%         figure(i)
        figure
        errorbar(mu,sigma,'rx')

    

end






    
