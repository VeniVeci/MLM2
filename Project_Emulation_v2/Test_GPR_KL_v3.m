% Test_GPR_KL_v3
%
% This script is written to track down the fail prombel when PCA_GPR is
% applied to the data set: "5050PoleFlow"
% 20-Nov-2013 

clear

%Assign parameter

new_dim=3;
num_test=100;

num_train_star=10;
num_train_step=10;
num_train_end=200;

%Assign memory for recoder
SSE=[];

%% PCA-GPR test
%----------------------
i=1;

for num_train = num_train_star:num_train_step:num_train_end
    
    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %------------------------------------------------------------------------------------------
    % Testing------------------------------
    
    [Y_star_svd,Yvar_star_svd,t_svd]=GPR_SVD(X,Y,X_star,new_dim);
    options = struct('ker','rbf','arg',500,'new_dim',new_dim); 
    [Y_star_kpca,Yvar_star_kpca,t_kpca]=GPR_KPCA(X,Y,X_star,options); %5000000
     
    
    %-------------------------------------
     
    SE =(Y_star_svd-Y_starorig).^2;
    SSE(i,1)=sum(SE(:));
    
    SE =(Y_star_kpca-Y_starorig).^2;
    SSE(i,2)=sum(SE(:));

    i=i+1;
    
end


    
    