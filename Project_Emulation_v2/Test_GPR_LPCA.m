% Test_GPR_LPCA
%
% This script is written to track down the fail prombel when PCA_GPR is
% applied to the data set: "5050PoleFlow"
% 20-Nov-2013 

clear

%Assign parameter
new_dim_start=4;    
new_dim_step =1;
new_dim_end  =5;

num_train=200;                  %Defining the input parameters
num_test=100;

num_train_star=10;
num_train_step=10;
num_train_end=100;

%Assign memory for recoder
SSE=[];

% %% PCA test
% 
% new_dim=5;
% for i = 1:new_dim
%     RecZ{i}=[];
% end
% i=1;
% 
% for num_train = num_train_star:num_train_step:num_train_end
%     
%     %------------------------------------------------------------------------------------
%     %Defining the input
% 
%     [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);
% 
%     [np_train,Dim_X]=size(X);
%     [np_test, Dim_Y]=size(Y_starorig);
% 
%     %------------------------------------------------------------------------------------------
%     % Testing------------------------------
%     
%     [Z,Key] = DataReduc_SVD(Y,new_dim);
%     
%     for j=1:new_dim
%          Z_temp=nan*ones(num_train_end,1);
%          Z_temp(1:num_train,1)=Z(:,j);
%          RecZ{j}=[RecZ{j},Z_temp];             %RecZ{j} means Z for jth principal component
%     end
%     
% end
% 

%% PCA-GPR test
%----------------------
% parameters
new_dim=2;

for i = 1:new_dim
    RecZstar{i}=[];
    RecZ{i}=[];
end

i=1;

for num_train = num_train_star:num_train_step:num_train_end
    
    %------------------------------------------------------------------------------------
    %Defining the input

    [X,Y,X_star,Y_starorig]=DataGenerate(num_train,num_test);

    [np_train,Dim_X]=size(X);
    [np_test, Dim_Y]=size(Y_starorig);

    %------------------------------------------------------------------------------------------
    % Testing------------------------------
    
    [Z,Key] = DataReduc_SVD(Y,new_dim);
    [Z_star,~]= GPR_prediction(X,Z,X_star);
    
    %-------------------------------------
    % Recording the process
    for j=1:new_dim
         Z_temp=nan*ones(num_train_end,1);
         Z_temp(1:num_train,1)=Z(:,j);
         RecZ{j}=[RecZ{j},Z_temp];             %RecZ{j} means Z for jth principal component
    end
      
    for j=1:new_dim                
         Zstar_temp=Z_star(:,j); % Assign memory
         RecZstar{j}=[RecZstar{j},Zstar_temp];  
    end
    %-------------------------------------
    
    Y_star=Z_star*Key;
    
    SE =(Y_star-Y_starorig).^2;
    SSE(i,1)=sum(SE(:));
    i=i+1;
    
end


% SSE=[SSE,SSE_SVD];

% SSE=[SSE_KPCA,SSE_SVD];
    
    
    